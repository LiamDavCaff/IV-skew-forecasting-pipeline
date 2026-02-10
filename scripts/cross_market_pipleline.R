# ===========================================================
# MULTI-INDEX IV SKEW OOS FORECASTING WITH COMMON START DATE
# ------------------------------------------------------------
# - Uses full monthly builder (all IV composites, macro, wings)
# - Computes IV coverage (80–120%)
# - Finds common usable start month
# - Builds monthly panels for each index
# - Runs OOS R² + Clark–West for chosen IV skews
# - Produces heatmap: indices × horizons × window lengths,
#   faceted by predictor 
# - Multiple-testing: BH-adjust CW p-values within each (benchmark, window_type, horizon)
# - Significance = (q_BH <= FDR_LEVEL) & (R2 > 0)
# ===========================================================

# ---- 0) Packages ---------------------------------------------------------

pkgs <- c(
  "here",
  "dplyr","tidyr","ggplot2","lubridate","zoo","slider","scales","stringr",
  "cowplot","scico","ggtext","lmtest","sandwich","broom","np",
  "knitr","kableExtra","tibble","purrr","readr","roll"
)

to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages(
  invisible(lapply(pkgs, library, character.only = TRUE))
)

# ---- 0.1) Project root + output folders -----------------------------------------------------

ROOT <- here::here()

OUT_FIG_XM <- file.path(ROOT, "outputs", "cross_market", "figures")
OUT_TAB_XM <- file.path(ROOT, "outputs", "cross_market", "appendix")

dir.create(OUT_FIG_XM, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB_XM, recursive = TRUE, showWarnings = FALSE)

# ---- 1) load data -----------------------------------------------------------

INDEX_KEYS <- c("SPX","FTSE","DAX","ASX","ESTOX","KOSPI","NASDAQ","NIFTY","NIKKEI","SMI")

# Default naming convention: data/bbg_<lowercase key>_data.rds
INDEX_PATHS <- setNames(
  file.path("data", paste0("bbg_", tolower(INDEX_KEYS), "_data.rds")),
  INDEX_KEYS
)

# Fail fast if missing
missing <- INDEX_PATHS[!file.exists(INDEX_PATHS)]
if (length(missing) > 0) {
  stop(
    paste0(
      "Missing cross-market data files in /data:\n",
      paste(names(missing), "->", missing, collapse = "\n")
    ),
    call. = FALSE
  )
}

# read all into memory (used by coverage step)
INDEX_LIST <- purrr::map(INDEX_PATHS, readRDS)

# ---- 2) Helpers ----------------------------------------------------------

# get last available observation of month (end of month)
last_non_na <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) NA_real_ else tail(y, 1)
}

# Snap EOM macro prints back ≤ maxgap days (handles weekend EOM), then LOCF (no peek)
snap_back_then_ffill <- function(x, date, maxgap = 3) {
  x_nocb <- zoo::na.locf(x, fromLast = TRUE, maxgap = maxgap, na.rm = FALSE)
  zoo::na.locf(x_nocb, na.rm = FALSE)
}

# Rolling mean using cumsums, then lag by 1 to avoid peek
roll_mean_lag1 <- function(y, k) {
  y2  <- replace(y, !is.finite(y), 0)
  cnt <- as.integer(is.finite(y))
  cs  <- c(0, cumsum(y2))
  cc  <- c(0, cumsum(cnt))
  m   <- (cs[-1] - dplyr::lag(cs, k, default = 0)[-1]) /
    (cc[-1] - dplyr::lag(cc, k, default = 0)[-1])
  dplyr::lag(m, 1)
}

# Expanding mean using only info through t-1 (no peek)
exp_mean_lag1 <- function(y) {
  y2  <- ifelse(is.finite(y), y, 0)
  cnt <- as.integer(is.finite(y))
  cs  <- cumsum(y2)
  cc  <- cumsum(cnt)
  m   <- ifelse(cc > 0, cs / cc, NA_real_)
  dplyr::lag(m, 1)
}

# Indices for rolling or expanding estimation window
get_window_indices <- function(i, window_months, window_type = c("rolling", "expanding")) {
  window_type <- match.arg(window_type)
  if (window_type == "rolling") {
    start <- i - window_months
    if (start < 1) return(NULL)
    return(start:(i - 1))
  } else {
    if (i <= window_months) return(NULL)
    return(1:(i - 1))
  }
}

# Multiple testing bias: Used to correct clark-west p values in OOS
bh_adjust <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], method = "BH")
  out
}


# ---- 3) Prepare Monthly Panel -----------------------------------------------------------

prepare_monthly_for_forecasting <- function(df_daily,
                                            date_col  = "date",
                                            price_col = "total_price", # TRI (total return inde( price + dividends))
                                            spot_col  = "price",       # spot index
                                            rf_col    = "rfr_1m",      # annualised %, 1m
                                            vix_col   = "vix",
                                            m80_col   = "80%",
                                            m90_col   = "90%",
                                            m100_col  = "100%",
                                            m110_col  = "110%",
                                            m120_col  = "120%",
                                            eps_col   = "eps_ttm",
                                            bvps_col  = "book",
                                            t3m_col   = "tbl_3m",
                                            lty_col   = "lty",
                                            infl_col  = "inflation_yoy",
                                            ts_col    = "term_spread",
                                            var1m_col = "stock_var_1m",
                                            rf_basis  = 360) {
  
  eps   <- .Machine$double.eps # constant to avod log(0)
  m_vec <- c(0.8, 0.9, 1.0, 1.1, 1.2) #moneyness grid for the quadratic fit in log-moneyness space
  
  #daily prep: align dayes, build daily log price and daily rf pieces
  Dd <- df_daily %>%
    dplyr::arrange(.data[[date_col]]) %>%
    dplyr::mutate(
      
      #standardise date and pull key price series
      date        = as.Date(.data[[date_col]]),
      price       = .data[[price_col]],     # total return index (TRI)
      spot_price  = .data[[spot_col]],      # spot index for valuation ratios
      logp_d      = log(price),
      
      # Risk free rate converted to decimals and LOCF fill
      rf_ann      = zoo::na.locf(.data[[rf_col]] / 100, na.rm = FALSE),
      
      # next calandar day (needed for day count fraction)
      next_date   = dplyr::lead(date),
      
      # rf basis day count fraction from date to next_date
      dcf         = as.numeric(next_date - date) / rf_basis,
      
      # onde day log risk free piece = log(1+rf *dcf)
      rf_piece    = dplyr::if_else(is.na(dcf), NA_real_, log1p(rf_ann * dcf)),
      
      #month start tags the current days month; month end tags month where accrual ends
      month_start = floor_date(date, "month"),
      month_end   = floor_date(next_date, "month")
    )
  
  # Monthly RF: sum daily rf peices whose accrual ends in month m
  rf_by_endmonth <- Dd %>%
    dplyr::group_by(month = month_end) %>%
    dplyr::summarise(
      rf_m = {
        x <- stats::na.omit(rf_piece)
        if (length(x) == 0) NA_real_ else sum(x)
      },
      .groups = "drop"
    )
  
  # Daily feature engineering
  M <- Dd %>%
    dplyr::mutate(
      # carry forward any slow moving fundamentals
      eps_ttm     = zoo::na.locf(.data[[eps_col]], na.rm = FALSE),
      book        = zoo::na.locf(.data[[bvps_col]], na.rm = FALSE),
      t3m_yield   = zoo::na.locf(.data[[t3m_col]],   na.rm = FALSE),
      lty_yield   = zoo::na.locf(.data[[lty_col]],   na.rm = FALSE),
      term_spread = zoo::na.locf(.data[[ts_col]],    na.rm = FALSE),
      
      #Inflation: snap back if stale and then fill forward with gap constraint
      infl_step   = snap_back_then_ffill(.data[[infl_col]], date, maxgap = 3),
      
      # Daily valuation ratios
      log_ep_d    = log(pmax(eps_ttm, eps) / pmax(spot_price, eps)),
      log_bm_d    = log(pmax(book,   eps) / pmax(spot_price, eps)),
      
      # Daily IV levels converted to decimals
      iv_80_d     = .data[[m80_col]]  / 100,
      iv_90_d     = .data[[m90_col]]  / 100,
      iv_100_d    = .data[[m100_col]] / 100,
      iv_110_d    = .data[[m110_col]] / 100,
      iv_120_d    = .data[[m120_col]] / 100,
      vix_ann_d   = .data[[vix_col]]  / 100,
      
      # Daily realised variance and volatility
      var_1m_d    = .data[[var1m_col]],
      rv_1m_d     = sqrt(var_1m_d)
    ) %>%
    
    # collapse daily data to month end
    dplyr::group_by(month = month_start) %>%
    dplyr::summarise(
      iv_80         = last_non_na(iv_80_d),
      iv_90         = last_non_na(iv_90_d),
      iv_100        = last_non_na(iv_100_d),
      iv_110        = last_non_na(iv_110_d),
      iv_120        = last_non_na(iv_120_d),
      vix_ann       = last_non_na(vix_ann_d),
      var_1m_mth    = last_non_na(var_1m_d),
      rv_1m_mth     = last_non_na(rv_1m_d),
      
      # Monthly price level and monthly log price
      price         = last_non_na(price),
      rfr_1m        = last_non_na(.data[[rf_col]]),
      logp          = log(price),
      
      # month end valuation ratios/macro variables
      log_ep        = last_non_na(log_ep_d),
      log_bm        = last_non_na(log_bm_d),
      t3m_yield     = last_non_na(t3m_yield),
      lty_yield     = last_non_na(lty_yield),
      term_spread   = last_non_na(term_spread),
      inflation_yoy = last_non_na(infl_step),
      .groups       = "drop"
    ) %>%
    
    # add monthly rd returns
    dplyr::arrange(month) %>%
    dplyr::left_join(rf_by_endmonth, by = "month")
  
  # Monthly skew measures 
  M <- M %>%
    dplyr::mutate(
      
      # Skew differences
      iv_skew_80_100    = iv_80 - iv_100,
      iv_skew_80_110    = iv_80 - iv_110,
      iv_skew_80_120    = iv_80 - iv_120,
      iv_skew_90_100    = iv_90 - iv_100,
      iv_skew_90_110    = iv_90 - iv_110,
      iv_skew_90_120    = iv_90 - iv_120,
      
      #Skew ratios
      skew_ratio_80_100 = log(pmax(iv_80,  eps) / pmax(iv_100, eps)),
      skew_ratio_90_100 = log(pmax(iv_90,  eps) / pmax(iv_100, eps)),
      skew_ratio_80_110 = log(pmax(iv_80,  eps) / pmax(iv_110, eps)),
      skew_ratio_90_110 = log(pmax(iv_90,  eps) / pmax(iv_110, eps)),
      skew_ratio_80_120 = log(pmax(iv_80,  eps) / pmax(iv_120, eps)),
      skew_ratio_90_120 = log(pmax(iv_90,  eps) / pmax(iv_120, eps)),
      
      # Wing slope/curvature in levels across moneyness (0.8 -> 1.2 is width 0.4)
      wing_slope  = (iv_120 - iv_80) / 0.4,
      wing_curve  = (iv_120 - 2*iv_100 + iv_80) / (0.2^2),
      
      # Convert annualised VIX to 1-month implied vol/variance (21 trading days)
      iv_imp_1m   = vix_ann * sqrt(21/252),
      var_imp_1m  = vix_ann^2 * (21/252),
      
      # Volatility/variance risk premia: implied minus realised
      vol_prem_1m = iv_imp_1m  - rv_1m_mth,
      var_prem_1m = var_imp_1m - var_1m_mth
    )
  
  # Quadratic fit in log moneyness space
  # Store b as slope and 2c as curavture
  mk <- log(m_vec)
  
  # For each month. fit quadratic to log IV across moneyness pts 
  coefs_log <- purrr::pmap(list(M$iv_80, M$iv_90, M$iv_100, M$iv_110, M$iv_120), ~{
    v <- log(pmax(c(..1, ..2, ..3, ..4, ..5), eps))
    if (all(is.finite(v))) lm(v ~ mk + I(mk^2))$coefficients
    else c("(Intercept)" = NA_real_, "mk" = NA_real_, "I(mk^2)" = NA_real_)
  })
  
  #Extract b and c terms across months
  b_log <- purrr::map_dbl(coefs_log, ~ .x[["mk"]])
  c_log <- purrr::map_dbl(coefs_log, ~ .x[["I(mk^2)"]])
  
  # Save monthly slope/curvature
  M$log_slope_quad <- b_log
  M$log_curve_quad <- 2 * c_log
  
  M
}


# ---- 4) Coverage Check & Common Start Date --------------------------------------

INDEX_LIST <- purrr::map(INDEX_PATHS, readRDS)

iv_cols <- c("80%", "90%", "100%", "110%", "120%")

# Use to check when we have full IV across moneyness pts 
iv_coverage_daily <- purrr::map_dfr(names(INDEX_LIST), function(idx) {
  df <- INDEX_LIST[[idx]] %>%
    dplyr::mutate(
      date     = as.Date(date),
      have_all = dplyr::if_all(dplyr::all_of(iv_cols), ~ is.finite(.x))
    )
  
  any_start <- min(df$date)
  any_end   <- max(df$date)
  full_df   <- df %>% dplyr::filter(have_all)
  
  if (nrow(full_df) == 0) {
    return(
      tibble::tibble(
        index             = idx,
        any_start         = any_start,
        any_end           = any_end,
        first_full_80_120 = as.Date(NA),
        last_full_80_120  = as.Date(NA),
        n_days_full       = 0L,
        coverage_share    = 0
      )
    )
  }
  
  tibble::tibble(
    index             = idx,
    any_start         = any_start,
    any_end           = any_end,
    first_full_80_120 = min(full_df$date),
    last_full_80_120  = max(full_df$date),
    n_days_full       = nrow(full_df),
    coverage_share    = nrow(full_df) / nrow(df)
  )
})

latest_first_full  <- max(iv_coverage_daily$first_full_80_120, na.rm = TRUE)
COMMON_START_MONTH <- floor_date(latest_first_full, "month")
message(">>> COMMON_START_MONTH = ", COMMON_START_MONTH)

# Write to latex

TabC3 <- iv_coverage_daily %>%
  mutate(
    any_start         = format(any_start, "%Y-%m-%d"),
    any_end           = format(any_end, "%Y-%m-%d"),
    first_full_80_120 = format(first_full_80_120, "%Y-%m-%d"),
    last_full_80_120  = format(last_full_80_120, "%Y-%m-%d"),
    coverage_share    = round(coverage_share, 3)
  ) %>%
  knitr::kable(
    format   = "latex",
    booktabs = TRUE,
    caption  = "Daily IV coverage by index (finite IV at 80--120\\% moneyness).",
    label    = "tab:iv_coverage_daily",
    align    = "lcccccr"
  ) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"))

kableExtra::save_kable(
  TabC3,
  file = file.path(OUT_TAB_XM, "TabC1_iv_coverage_daily.tex")
)


# ---- 5) Add Excess Returns ------------------------------------------------------

add_excess_returns_monthly <- function(M, horizons_months = c(1, 3, 6, 12)) {
  out    <- M
  rf_fwd <- dplyr::lead(out$rf_m, 1)
  for (h in horizons_months) {
    cum_rf <- slider::slide_index_dbl(
      rf_fwd, out$month, sum,
      .before = 0, .after = h - 1, .complete = TRUE
    )
    out[[paste0("excess_", h, "m")]] <- dplyr::lead(out$logp, h) - out$logp - cum_rf
  }
  out
}


# ---- 6) Function to Create Monthly Panel for Index ------------------------------

build_monthly_for_index <- function(df_daily,
                                    common_start_month = COMMON_START_MONTH,
                                    horizons_months    = HORIZONS_MONTHS) {
  df_daily %>%
    prepare_monthly_for_forecasting() %>%
    add_excess_returns_monthly(horizons_months) %>%
    dplyr::filter(month >= common_start_month)
}

# ---- 7) OOS forecasting for 1 index ---------------------------------------------

get_ols_forecasts_univariate <- function(y, x, h_months, window_months, window_type) {
  n      <- length(y) #number of observations
  
  #vectors to store forecasts and beta at each time
  f_mod  <- rep(NA_real_, n) 
  beta_t <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    
    # get estimation window incides for forecast origin
    idx <- get_window_indices(i, window_months, window_type)
    if (is.null(idx)) next
    
    # get estimation sample for y and x values for window
    yi <- y[idx]; xi <- x[idx]
    ok <- is.finite(yi) & is.finite(xi) # only keep rows where both y and x have data
    if (sum(ok) < window_months) next # require a full window of observations before estimating
    
    #fit the predictive regression y~x
    fit <- lm(yi[ok] ~ xi[ok])
    
    # extract coeficients cf[1] = intercept and cf[2] =slope
    cf  <- coef(fit)
    if (length(cf) < 2L) next
    
    # store slope estimates
    beta_t[i] <- cf[2]
    if (is.finite(x[i])) {
      f_mod[i] <- cf[1] + cf[2] * x[i]
    }
  }
  list(f_mod = f_mod, beta_t = beta_t)
}


oos_stats_one_index <- function(M,
                                horizons_months,
                                window_months,
                                pred,
                                benchmark,
                                window_type,
                                min_points) {
  
  # Loop over horizons h and return a named vector of OOS statistics for each horizon
  vapply(
    horizons_months,
    function(h) {
      
      # y variable: h-month ahead excess return series (e.g. "excess_12m")
      y <- M[[paste0("excess_", h, "m")]]
      
      # Univariate predictor series
      x <- M[[pred]]
      
      # Generate rolling/expanding univariate OLS forecasts:
      # f_mod = model forecast at each time t
      # b_t   = estimated slope coefficient (beta) at each time t
      ols   <- get_ols_forecasts_univariate(y, x, h, window_months, window_type)
      f_mod <- ols$f_mod
      b_t   <- ols$beta_t
      
      # Construct benchmark forecast:
      # rolling: rolling historical mean (lagged)
      # expanding: expanding historical mean (lagged)
      f_bm <- if (benchmark == "rolling") {
        roll_mean_lag1(y, window_months)
      } else {
        exp_mean_lag1(y)
      }
      
      # Keep dates where model forecast, benchmark forecast, and realised y are all available
      ok <- which(is.finite(f_mod) & is.finite(f_bm) & is.finite(y))
      
      #Minimum number of OOS observations to compute stable statistics
      if (length(ok) < min_points) {
        return(c(R2_raw=NA_real_, cw_stat=NA_real_, cw_p=NA_real_, beta_avg=NA_real_))
      }
      
      # Realised returns and forecast errors on the valid OOS dates
      y_ok <- y[ok]
      e_m  <- y_ok - f_mod[ok]  # model error
      e_b  <- y_ok - f_bm[ok]   # benchmark error
      
      # Campbell-thompson style OOS R^2
      R2_raw <- 100 * (1 - sum(e_m^2) / sum(e_b^2))
      
      # Clark–West adjusted loss differential for nested models
      adj_loss <- e_b^2 - (e_m^2 - (f_bm[ok] - f_mod[ok])^2)
      
      # Estimate mean(adj_loss) via intercept-only regression
      cw_fit   <- lm(adj_loss ~ 1)
      
      # HAC variance of the intercept using Newey–West with lag = h 
      V        <- sandwich::NeweyWest(cw_fit, lag = h, prewhite = FALSE)
      
      # Clark–West test statistic and one-sided p-value
      cw_stat  <- coef(cw_fit)[1] / sqrt(V[1, 1])
      cw_p     <- pnorm(cw_stat, lower.tail = FALSE)
      
      # Average beta over the OOS-valid evaluation period (useful for interpreting sign/magnitude)
      beta_avg <- mean(b_t[ok], na.rm = TRUE)
      
      # Return stats for this horizon
      c(R2_raw=R2_raw, cw_stat=cw_stat, cw_p=cw_p, beta_avg=beta_avg)
    },
    
    # vapply template: ensures numeric output with consistent names
    FUN.VALUE = c(R2_raw=NA_real_, cw_stat=NA_real_, cw_p=NA_real_, beta_avg=NA_real_)
  )
}



# ---- 8) Run Forecasts Across All Indices ----------------------------------------

run_oos_across_indices <- function(index_files, predictor,
                                   windows_months, horizons_months,
                                   benchmark, window_type,
                                   min_points) {
  
  purrr::imap_dfr(index_files, function(path, idx) {
    message(">>> Building & forecasting: ", idx)
    
    df_daily <- readRDS(path)
    M        <- build_monthly_for_index(df_daily, COMMON_START_MONTH, horizons_months)
    
    tidyr::crossing(window_months = windows_months,
                    horizon       = horizons_months) %>%
      purrr::pmap_dfr(function(window_months, horizon) {
        
        o <- oos_stats_one_index(
          M               = M,
          horizons_months = horizon,
          window_months   = window_months,
          pred            = predictor,
          benchmark       = benchmark,
          window_type     = window_type,
          min_points      = min_points
        )
        
        tibble::tibble(
          index         = idx,
          predictor     = predictor,
          window_months = window_months,
          horizon       = paste0(horizon, "m"),
          benchmark     = benchmark,
          window_type   = window_type,
          R2_oos_raw    = as.numeric(o["R2_raw", ]),
          cw_p_raw      = as.numeric(o["cw_p", ]),     # keep raw p-values for BH
          beta_avg      = as.numeric(o["beta_avg", ])
        )
      })
  })
}
 
# ---- 9) Create Heatmap ----------------------------------------------------------

plot_heatmap_indices_both <- function(results_both,
                                      predictors_cross,
                                      alpha_txt = 0.85,
                                      fdr_level = 0.10) {
  
  pretty_pred <- function(x) stringr::str_replace_all(x, "_", " ")
  
  df <- results_both %>%
    dplyr::mutate(
      horizon_num   = readr::parse_number(horizon),
      horizon_fac   = factor(horizon_num, levels = sort(unique(horizon_num))),
      window_months = factor(window_months, levels = sort(unique(window_months))),
      predictor     = factor(predictor, levels = predictors_cross),
      signif_flag   = is.finite(cw_q_h) & (cw_q_h <= fdr_level) & (R2_oos_raw > 0),
      R2_cap        = pmax(R2_oos_raw, -50)
    )
  
  idx_order <- df %>%
    dplyr::filter(horizon_num == 12) %>%
    dplyr::group_by(index) %>%
    dplyr::summarise(avg_R2_12 = mean(R2_oos_raw, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(avg_R2_12)) %>%
    dplyr::pull(index)
  
  df$index <- factor(df$index, levels = idx_order)
  
  ggplot(df, aes(x = horizon_fac, y = index, fill = R2_cap)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(
      aes(
        label    = ifelse(is.finite(R2_oos_raw) & abs(R2_oos_raw) >= 0.1,
                          sprintf("%.1f", R2_oos_raw), ""),
        colour   = signif_flag,
        fontface = ifelse(signif_flag, "bold", "plain")
      ),
      size  = 3,
      alpha = alpha_txt
    ) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scico::scale_fill_scico(
      palette = "vikO",
      limits  = range(df$R2_cap, na.rm = TRUE),
      breaks  = scales::pretty_breaks(),
      name    = expression(R[OOS]^2~"(%)"),
      guide   = guide_colourbar(
        direction      = "horizontal",
        barwidth       = unit(7, "cm"),
        barheight      = unit(0.45, "cm"),
        title.position = "top",
        title.hjust    = 0.5
      )
    ) +
    labs(x = "Forecast horizon (months)", y = NULL) +
    facet_grid(. ~ predictor, labeller = labeller(predictor = pretty_pred)) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid        = element_blank(),
      axis.ticks        = element_blank(),
      legend.position   = "bottom",
      legend.key.height = unit(0.45, "cm"),
      strip.text        = element_text(face = "bold", size = 10),
      axis.text.x       = element_text(size = 8),
      axis.text.y       = element_text(size = 10)
    )
}


# ---- 10) Expanding Window Results ------------------------------

PREDICTORS_CROSS_EXP <- c("log_slope_quad", "skew_ratio_90_120")
WINDOWS_MONTHS_EXP   <- c(120)
HORIZONS_MONTHS_EXP  <- c(1, 3, 6, 12)

BENCHMARK_TYPE_EXP   <- "expanding"   
WINDOW_TYPE_OOS_EXP  <- "expanding"   
MIN_FORECASTS    <- 12            # min number of OOS points required
FDR_LEVEL        <- 0.1          # BH FDR threshold used for significance marking

results_exp <- purrr::map_dfr(PREDICTORS_CROSS_EXP, function(pred) {
  run_oos_across_indices(
    index_files     = INDEX_PATHS,
    predictor       = pred,
    windows_months  = WINDOWS_MONTHS_EXP,
    horizons_months = HORIZONS_MONTHS_EXP,
    benchmark       = BENCHMARK_TYPE_EXP,
    window_type     = WINDOW_TYPE_OOS_EXP,
    min_points      = MIN_FORECASTS
  )
})

# MULTIPLE TESTING  (same logic as your monthly suite
# Within each (benchmark, window_type, horizon) slice, adjust CW p-values
# across all indices × predictors × windows in that slice.
results_exp <- results_exp %>%
  group_by(benchmark, window_type, window_months, predictor,horizon) %>%
  mutate(cw_q_h = bh_adjust(cw_p_raw)) %>%
  ungroup()


fig17_cross_exp <- plot_heatmap_indices_both(
  results_exp,
  predictors_cross = PREDICTORS_CROSS_EXP,
  fdr_level = FDR_LEVEL
)

ggsave(
  filename = file.path(OUT_FIG_XM, "fig17_cross_heatmap_expanding.png"),
  plot     = fig17_cross_exp,
  width    = 8, height = 4.5, dpi = 300
)

#Appendix for significant results

sig_report_exp <- results_exp %>%
  mutate(signif_flag = is.finite(cw_q_h) & (cw_q_h <= FDR_LEVEL) & (R2_oos_raw > 0)) %>%
  filter(signif_flag) %>%
  transmute(
    index, predictor, window_months, horizon, benchmark, window_type,
    R2_oos   = round(R2_oos_raw, 2),
    beta_avg = round(beta_avg, 4),
    q_BH     = signif(cw_q_h, 3)
  ) %>%
  arrange(predictor, horizon, desc(R2_oos))

h_levels <- c("1m","3m","6m","12m")

sig_exp <- sig_report_exp %>%
  mutate(cell = sprintf("R2=%.1f, b=%.3f, q=%.3f", R2_oos, beta_avg, q_BH)) %>%
  select(index, predictor, window_months, horizon, cell) %>%
  tidyr::pivot_wider(names_from = horizon, values_from = cell) %>%
  dplyr::select(index, predictor, window_months, dplyr::any_of(h_levels))

TabC4 <- knitr::kable(
  sig_exp,
  format    = "latex",
  booktabs  = TRUE,
  longtable = TRUE,
  escape    = TRUE   
)

writeLines(TabC4, file.path(OUT_TAB_XM, "TabC2_cross_market_significant_exp.tex"))

# ---- 11) Rolling Window Results ------------------------------

PREDICTORS_CROSS_ROLL <- c("skew_ratio_80_120","skew_ratio_80_110","log_slope_quad")
WINDOWS_MONTHS_ROLL   <- c(60)
HORIZONS_MONTHS_ROLL  <- c(1, 3, 6, 12)

BENCHMARK_TYPE_ROLL   <- "rolling"   
WINDOW_TYPE_OOS_ROLL  <- "rolling"   
MIN_FORECASTS    <- 12            # min number of OOS points required
FDR_LEVEL        <- 0.1          # BH FDR threshold used for significance marking

results_roll <- purrr::map_dfr(PREDICTORS_CROSS_ROLL, function(pred) {
  run_oos_across_indices(
    index_files     = INDEX_PATHS,
    predictor       = pred,
    windows_months  = WINDOWS_MONTHS_ROLL,
    horizons_months = HORIZONS_MONTHS_ROLL,
    benchmark       = BENCHMARK_TYPE_ROLL,
    window_type     = WINDOW_TYPE_OOS_ROLL,
    min_points      = MIN_FORECASTS
  )
})

# MULTIPLE TESTING  (same logic as your monthly suite
# Within each (benchmark, window_type, horizon) slice, adjust CW p-values
# across all indices × predictors × windows in that slice.
results_roll <- results_roll %>%
  group_by(benchmark, window_type, window_months, predictor,horizon) %>%
  mutate(cw_q_h = bh_adjust(cw_p_raw)) %>%
  ungroup()


fig18_cross_roll <- plot_heatmap_indices_both(
  results_roll,
  predictors_cross = PREDICTORS_CROSS_ROLL,
  fdr_level = FDR_LEVEL
)

ggsave(
  filename = file.path(OUT_FIG_XM, "fig18_cross_heatmap_rolling.png"),
  plot     = fig18_cross_roll,
  width    = 8, height = 4.5, dpi = 300
)

#Appendix for significant results

sig_report_roll <- results_roll %>%
  mutate(signif_flag = is.finite(cw_q_h) & (cw_q_h <= FDR_LEVEL) & (R2_oos_raw > 0)) %>%
  filter(signif_flag) %>%
  transmute(
    index, predictor, window_months, horizon, benchmark, window_type,
    R2_oos   = round(R2_oos_raw, 2),
    beta_avg = round(beta_avg, 4),
    q_BH     = signif(cw_q_h, 3)
  ) %>%
  arrange(predictor, horizon, desc(R2_oos))

h_levels <- c("1m","3m","6m","12m")

sig_roll <- sig_report_roll %>%
  mutate(cell = sprintf("R2=%.1f, b=%.3f, q=%.3f", R2_oos, beta_avg, q_BH)) %>%
  select(index, predictor, window_months, horizon, cell) %>%
  tidyr::pivot_wider(names_from = horizon, values_from = cell) %>%
  dplyr::select(index, predictor, window_months, dplyr::any_of(h_levels))

TabC5 <- knitr::kable(
  sig_roll,
  format    = "latex",
  booktabs  = TRUE,
  longtable = TRUE,
  escape    = TRUE   
)

writeLines(TabC5, file.path(OUT_TAB_XM, "TabC3_cross_market_significant_roll.tex"))