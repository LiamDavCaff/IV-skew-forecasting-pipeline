############################################################
# MULTI-INDEX IV SKEW OOS FORECASTING WITH AUTO-ALIGNED DATES
# ------------------------------------------------------------
# - Uses full monthly builder (all IV composites, macro, wings)
# - Computes IV coverage (80–120%)
# - Finds common usable start month
# - Builds monthly panels for each index
# - Runs OOS R² + Clark–West for chosen IV skews
# - Produces heatmap: indices × horizons × window lengths,
#   faceted by predictor (auto-labelled)
# - Multiple-testing: BH-adjust CW p-values within each (benchmark, window_type, horizon)
#   -> significance = (q_BH <= FDR_LEVEL) & (R2 > 0)
############################################################


###########################
# 0) USER SETTINGS
###########################

PREDICTORS_CROSS <- c( "skew_ratio_80_120", "skew_ratio_80_110","log_slope_quad")
WINDOWS_MONTHS   <- c(60)
HORIZONS_MONTHS  <- c(1, 3, 6, 12)

BENCHMARK_TYPE   <- "rolling"   # "rolling" or "expanding"
WINDOW_TYPE_OOS  <- "rolling"   # "rolling" or "expanding"
MIN_FORECASTS    <- 12            # min number of OOS points required

FDR_LEVEL        <- 0.1          # BH FDR threshold used for significance marking

# ---- Paths to RDS files ----
INDEX_FILES <- c(
  SPX    = "C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/spx_iv_modelling_data.rds",
  FTSE   = "C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/ftse_iv_modelling_data.rds",
  DAX    = "C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/dax_iv_modelling_data.rds",
  ASX    = "C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/asx_iv_modelling_data.rds",
  ESTOX  = "C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/eurostox_iv_modelling_data.rds",
  KOSPI  = "C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/kospi_iv_modelling_data.rds",
  NASDQ  = "C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/nasdaq_iv_modelling_data.rds",
  NIFTY  = "C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/nifty_iv_modelling_data.rds",
  NIKKEI = "C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/nikkei_iv_modelling_data.rds",
  SMI    = "C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/smi_iv_modelling_data.rds"
)


###########################
# 1) PACKAGES
###########################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(lubridate)
  library(zoo)
  library(slider)
  library(scico)
  library(stringr)
  library(sandwich)
  library(lmtest)
  library(readr)
  library(grid)   # unit()
})


############################################################
# 2) SMALL HELPERS
############################################################

last_non_na <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) NA_real_ else tail(y, 1)
}

snap_back_then_ffill <- function(x, date, maxgap = 3) {
  x_nocb <- zoo::na.locf(x, fromLast = TRUE, maxgap = maxgap, na.rm = FALSE)
  zoo::na.locf(x_nocb, na.rm = FALSE)
}

roll_mean_lag1 <- function(y, k) {
  y2  <- replace(y, !is.finite(y), 0)
  cnt <- as.integer(is.finite(y))
  cs  <- c(0, cumsum(y2))
  cc  <- c(0, cumsum(cnt))
  m   <- (cs[-1] - dplyr::lag(cs, k, default = 0)[-1]) /
    (cc[-1] - dplyr::lag(cc, k, default = 0)[-1])
  dplyr::lag(m, 1)
}

exp_mean_lag1 <- function(y) {
  y2  <- ifelse(is.finite(y), y, 0)
  cnt <- as.integer(is.finite(y))
  cs  <- cumsum(y2)
  cc  <- cumsum(cnt)
  m   <- ifelse(cc > 0, cs / cc, NA_real_)
  dplyr::lag(m, 1)
}

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


############################################################
# 3) FULL DAILY → MONTHLY BUILDER (YOUR VERSION)
############################################################

prepare_monthly_for_forecasting <- function(df_daily,
                                            date_col  = "date",
                                            price_col = "total_price", # TRI
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
  
  eps   <- .Machine$double.eps
  m_vec <- c(0.8, 0.9, 1.0, 1.1, 1.2)
  
  Dd <- df_daily %>%
    dplyr::arrange(.data[[date_col]]) %>%
    dplyr::mutate(
      date        = as.Date(.data[[date_col]]),
      price       = .data[[price_col]],
      spot_price  = .data[[spot_col]],
      logp_d      = log(price),
      
      rf_ann      = zoo::na.locf(.data[[rf_col]] / 100, na.rm = FALSE),
      next_date   = dplyr::lead(date),
      dcf         = as.numeric(next_date - date) / rf_basis,
      rf_piece    = dplyr::if_else(is.na(dcf), NA_real_, log1p(rf_ann * dcf)),
      
      month_start = floor_date(date, "month"),
      month_end   = floor_date(next_date, "month")
    )
  
  rf_by_endmonth <- Dd %>%
    dplyr::group_by(month = month_end) %>%
    dplyr::summarise(
      rf_m = {
        x <- stats::na.omit(rf_piece)
        if (length(x) == 0) NA_real_ else sum(x)
      },
      .groups = "drop"
    )
  
  M <- Dd %>%
    dplyr::mutate(
      eps_ttm     = zoo::na.locf(.data[[eps_col]], na.rm = FALSE),
      book        = zoo::na.locf(.data[[bvps_col]], na.rm = FALSE),
      t3m_yield   = zoo::na.locf(.data[[t3m_col]],   na.rm = FALSE),
      lty_yield   = zoo::na.locf(.data[[lty_col]],   na.rm = FALSE),
      term_spread = zoo::na.locf(.data[[ts_col]],    na.rm = FALSE),
      infl_step   = snap_back_then_ffill(.data[[infl_col]], date, maxgap = 3),
      
      log_ep_d    = log(pmax(eps_ttm, eps) / pmax(spot_price, eps)),
      log_bm_d    = log(pmax(book,   eps) / pmax(spot_price, eps)),
      
      iv_80_d     = .data[[m80_col]]  / 100,
      iv_90_d     = .data[[m90_col]]  / 100,
      iv_100_d    = .data[[m100_col]] / 100,
      iv_110_d    = .data[[m110_col]] / 100,
      iv_120_d    = .data[[m120_col]] / 100,
      vix_ann_d   = .data[[vix_col]]  / 100,
      
      var_1m_d    = .data[[var1m_col]],
      rv_1m_d     = sqrt(var_1m_d)
    ) %>%
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
      
      price         = last_non_na(price),
      rfr_1m        = last_non_na(.data[[rf_col]]),
      logp          = log(price),
      log_ep        = last_non_na(log_ep_d),
      log_bm        = last_non_na(log_bm_d),
      t3m_yield     = last_non_na(t3m_yield),
      lty_yield     = last_non_na(lty_yield),
      term_spread   = last_non_na(term_spread),
      inflation_yoy = last_non_na(infl_step),
      .groups       = "drop"
    ) %>%
    dplyr::arrange(month) %>%
    dplyr::left_join(rf_by_endmonth, by = "month")
  
  M <- M %>%
    dplyr::mutate(
      iv_skew_80_100    = iv_80 - iv_100,
      iv_skew_80_110    = iv_80 - iv_110,
      iv_skew_80_120    = iv_80 - iv_120,
      iv_skew_90_100    = iv_90 - iv_100,
      iv_skew_90_110    = iv_90 - iv_110,
      iv_skew_90_120    = iv_90 - iv_120,
      
      skew_ratio_80_100 = log(pmax(iv_80,  eps) / pmax(iv_100, eps)),
      skew_ratio_90_100 = log(pmax(iv_90,  eps) / pmax(iv_100, eps)),
      skew_ratio_80_110 = log(pmax(iv_80,  eps) / pmax(iv_110, eps)),
      skew_ratio_90_110 = log(pmax(iv_90,  eps) / pmax(iv_110, eps)),
      skew_ratio_80_120 = log(pmax(iv_80,  eps) / pmax(iv_120, eps)),
      skew_ratio_90_120 = log(pmax(iv_90,  eps) / pmax(iv_120, eps)),
      
      wing_slope  = (iv_120 - iv_80) / 0.4,
      wing_curve  = (iv_120 - 2*iv_100 + iv_80) / (0.2^2),
      
      iv_imp_1m   = vix_ann * sqrt(21/252),
      var_imp_1m  = vix_ann^2 * (21/252),
      vol_prem_1m = iv_imp_1m  - rv_1m_mth,
      var_prem_1m = var_imp_1m - var_1m_mth
    )
  
  mk <- log(m_vec)
  coefs_log <- purrr::pmap(list(M$iv_80, M$iv_90, M$iv_100, M$iv_110, M$iv_120), ~{
    v <- log(pmax(c(..1, ..2, ..3, ..4, ..5), eps))
    if (all(is.finite(v))) lm(v ~ mk + I(mk^2))$coefficients
    else c("(Intercept)" = NA_real_, "mk" = NA_real_, "I(mk^2)" = NA_real_)
  })
  b_log <- purrr::map_dbl(coefs_log, ~ .x[["mk"]])
  c_log <- purrr::map_dbl(coefs_log, ~ .x[["I(mk^2)"]])
  
  M$log_slope_quad <- b_log
  M$log_curve_quad <- 2 * c_log
  
  M
}


############################################################
# 4) COVERAGE CHECK & COMMON START MONTH
############################################################

INDEX_LIST <- lapply(INDEX_FILES, readRDS)

iv_cols <- c("80%", "90%", "100%", "110%", "120%")

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

print(iv_coverage_daily)

latest_first_full  <- max(iv_coverage_daily$first_full_80_120, na.rm = TRUE)
COMMON_START_MONTH <- floor_date(latest_first_full, "month")
message(">>> COMMON_START_MONTH = ", COMMON_START_MONTH)


############################################################
# 5) ADD EXCESS RETURNS + MONTHLY BUILDER WITH TRUNCATION
############################################################

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

build_monthly_for_index <- function(df_daily,
                                    common_start_month = COMMON_START_MONTH,
                                    horizons_months    = HORIZONS_MONTHS) {
  df_daily %>%
    prepare_monthly_for_forecasting() %>%
    add_excess_returns_monthly(horizons_months) %>%
    dplyr::filter(month >= common_start_month)
}


############################################################
# 6) OOS FORECASTER
############################################################

get_ols_forecasts_univariate <- function(y, x, h_months, window_months, window_type) {
  n      <- length(y)
  f_mod  <- rep(NA_real_, n)
  beta_t <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    idx <- get_window_indices(i, window_months, window_type)
    if (is.null(idx)) next
    yi <- y[idx]; xi <- x[idx]
    ok <- is.finite(yi) & is.finite(xi)
    if (sum(ok) < window_months) next
    
    fit <- lm(yi[ok] ~ xi[ok])
    cf  <- coef(fit)
    if (length(cf) < 2L) next
    
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
  
  vapply(
    horizons_months,
    function(h) {
      
      y <- M[[paste0("excess_", h, "m")]]
      x <- M[[pred]]
      
      ols   <- get_ols_forecasts_univariate(y, x, h, window_months, window_type)
      f_mod <- ols$f_mod
      b_t   <- ols$beta_t
      
      f_bm <- if (benchmark == "rolling") {
        roll_mean_lag1(y, window_months)
      } else {
        exp_mean_lag1(y)
      }
      
      ok <- which(is.finite(f_mod) & is.finite(f_bm) & is.finite(y))
      if (length(ok) < min_points) {
        return(c(R2_raw=NA_real_, cw_stat=NA_real_, cw_p=NA_real_, beta_avg=NA_real_))
      }
      
      y_ok <- y[ok]
      e_m  <- y_ok - f_mod[ok]
      e_b  <- y_ok - f_bm[ok]
      
      R2_raw <- 100 * (1 - sum(e_m^2) / sum(e_b^2))
      
      adj_loss <- e_b^2 - (e_m^2 - (f_bm[ok] - f_mod[ok])^2)
      cw_fit   <- lm(adj_loss ~ 1)
      V        <- sandwich::NeweyWest(cw_fit, lag = h, prewhite = FALSE)
      cw_stat  <- coef(cw_fit)[1] / sqrt(V[1, 1])
      cw_p     <- pnorm(cw_stat, lower.tail = FALSE)
      
      beta_avg <- mean(b_t[ok], na.rm = TRUE)
      
      c(R2_raw=R2_raw, cw_stat=cw_stat, cw_p=cw_p, beta_avg=beta_avg)
    },
    FUN.VALUE = c(R2_raw=NA_real_, cw_stat=NA_real_, cw_p=NA_real_, beta_avg=NA_real_)
  )
}


############################################################
# 7) RUN FORECASTS ACROSS ALL INDICES
############################################################

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


############################################################
# 8) HEATMAP FUNCTION (AUTO PRETTY FACET LABELS)
############################################################

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


############################################################
# 9) EXECUTE FULL PIPELINE + BH ADJUSTMENT + PLOT
############################################################

results_both <- purrr::map_dfr(PREDICTORS_CROSS, function(pred) {
  run_oos_across_indices(
    index_files     = INDEX_FILES,
    predictor       = pred,
    windows_months  = WINDOWS_MONTHS,
    horizons_months = HORIZONS_MONTHS,
    benchmark       = BENCHMARK_TYPE,
    window_type     = WINDOW_TYPE_OOS,
    min_points      = MIN_FORECASTS
  )
})

# ---------- MULTIPLE TESTING FIX (same logic as your monthly suite) ----------
# Within each (benchmark, window_type, horizon) slice, adjust CW p-values
# across all indices × predictors × windows in that slice.
results_both <- results_both %>%
  group_by(benchmark, window_type, window_months, predictor,horizon) %>%
  mutate(cw_q_h = bh_adjust(cw_p_raw)) %>%
  ungroup()

print(results_both, n = nrow(results_both))

fig_cross_both <- plot_heatmap_indices_both(
  results_both,
  predictors_cross = PREDICTORS_CROSS,
  fdr_level = FDR_LEVEL
)

print(fig_cross_both)

# Optional save
ggsave("cross_market_roll_heatmap.png", fig_cross_both, width = 8, height = 4, dpi = 300)


#Table for significance

sig_report <- results_both %>%
  mutate(signif_flag = is.finite(cw_q_h) & (cw_q_h <= FDR_LEVEL) & (R2_oos_raw > 0)) %>%
  filter(signif_flag) %>%
  transmute(
    index, predictor, window_months, horizon, benchmark, window_type,
    R2_oos   = round(R2_oos_raw, 2),
    beta_avg = round(beta_avg, 4),
    q_BH     = signif(cw_q_h, 3)
  ) %>%
  arrange(predictor, horizon, desc(R2_oos))

print(sig_report, n = nrow(sig_report))

h_levels <- c("1m","3m","6m","12m")

sig_pretty <- sig_report %>%
  mutate(cell = sprintf("R2=%.1f, b=%.3f, q=%.3f", R2_oos, beta_avg, q_BH)) %>%
  select(index, predictor, window_months, horizon, cell) %>%
  tidyr::pivot_wider(names_from = horizon, values_from = cell) %>%
  dplyr::select(index, predictor, window_months, dplyr::any_of(h_levels))

sig_pretty







