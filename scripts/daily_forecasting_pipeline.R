# ===========================================================
#  Daily Forecasting Pipeline 
#  - Daily panel cached once
#  - Rolling / Expanding OLS via helper
#  - Benchmark via rolling / expanding mean
#  - Clark–West test
#  - Multiple-testing correction using BH
# ===========================================================

# ---- 0) Packages ---------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(zoo)
  library(slider)
  library(scales)
  library(stringr)
  library(cowplot)
  library(scico)
  library(ggtext)
  library(lmtest)
  library(sandwich)
  library(broom)
  library(np)
  suppressWarnings({ if (!requireNamespace("ggridges", quietly = TRUE)) NULL else library(ggridges) })
  library(knitr); library(kableExtra); library(tibble); library(purrr); library(readr)
})
if (!requireNamespace("roll", quietly = TRUE)) install.packages("roll")
library(roll)


# ---- 0.1) Output folders ------------------------------------------------------

FIG_DAILY <- file.path("outputs", "daily", "figures")
dir.create(FIG_DAILY, recursive = TRUE, showWarnings = FALSE)

# ---- 1) Load data --------------------------------------------------------

# Update path if needed:
DATA_FILE <- file.path("data", "bbg_spx_data.rds")

if (!file.exists(DATA_FILE)) {
  stop(
    paste0("Missing ", DATA_FILE, ".\n",
           "Place the file in the repo's data/ folder (it stays local / ignored)."),
    call. = FALSE
  )
}

data <- readRDS(DATA_FILE)
stopifnot("date" %in% names(data))

# ---- 2) Helpers ----------------------------------------------------------

last_non_na <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) NA_real_ else tail(y, 1)
}

# Snap EOM macro prints back ≤ maxgap days (handles weekend EOM), then LOCF (no peek)
snap_back_then_ffill <- function(x, date, maxgap = 3) {
  x_nocb <- zoo::na.locf(x, fromLast = TRUE, maxgap = maxgap, na.rm = FALSE) # NOCB
  zoo::na.locf(x_nocb, na.rm = FALSE)                                        # LOCF
}

# Rolling mean using cumsums, then lag by 1 to avoid peek
roll_mean_lag1 <- function(y, k) {
  y2  <- replace(y, !is.finite(y), 0)
  cnt <- as.integer(is.finite(y))
  cs  <- c(0, cumsum(y2))
  cc  <- c(0, cumsum(cnt))
  m <- (cs[-1] - dplyr::lag(cs, k, default = 0)[-1]) /
    (cc[-1] - dplyr::lag(cc, k, default = 0)[-1])
  dplyr::lag(m, 1)
}

# Expanding mean using only info through t-1 (no peek)
exp_mean_lag1 <- function(y) {
  y2  <- ifelse(is.finite(y), y, 0)
  cnt <- as.integer(is.finite(y))
  cs  <- cumsum(y2)
  cc  <- cumsum(cnt)
  m   <- ifelse(cc > 0, cs/cc, NA_real_)
  dplyr::lag(m, 1)
}

# Helper: indices for rolling or expanding estimation window
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

#  OLS forecasts for a given window type (rolling / expanding)
get_ols_forecasts <- function(
    y,
    X,
    h_days,
    window_days,
    window_type = c("rolling", "expanding"),
    min_obs     = window_days
) {
  window_type <- match.arg(window_type)
  n <- length(y)
  k <- ncol(X)
  
  f_mod <- rep(NA_real_, n)
  betas <- matrix(NA_real_, nrow = n, ncol = k)
  
  if (window_type == "rolling") {
    coef_mat <- roll::roll_lm(
      X, y,
      width        = window_days,
      intercept    = TRUE,
      min_obs      = min_obs,
      complete_obs = TRUE
    )$coef
    
    b0 <- dplyr::lag(coef_mat[, 1], 1)
    Bx <- dplyr::lag(coef_mat[, -1, drop = FALSE], 1)
    
  } else { # expanding
    coef_mat <- roll::roll_lm(
      X, y,
      width        = n,
      intercept    = TRUE,
      min_obs      = min_obs,
      complete_obs = TRUE
    )$coef
    
    b0 <- dplyr::lag(coef_mat[, 1], 1)
    Bx <- dplyr::lag(coef_mat[, -1, drop = FALSE], 1)
  }
  
  # forecasts (where we have coefficients and X is finite)
  ok <- which(is.finite(b0) & apply(Bx, 1, function(r) all(is.finite(r))) &
                apply(X, 1, function(r) all(is.finite(r))))
  
  f_mod[ok] <- b0[ok] + rowSums(Bx[ok, , drop = FALSE] * X[ok, , drop = FALSE])
  betas[ok, ] <- Bx[ok, , drop = FALSE]
  
  list(f_mod = f_mod, betas = betas)
}

# Multiple testing bias: Used to correct clark-west p values in OOS
bh_adjust <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

# ---- 3) Daily panel (predictors + targets) ----------------------------

# Build month-end panel from daily data, including RF accrual timing

prepare_daily_for_forecasting <- function(df_daily,
                                          date_col  = "date",
                                          price_col = "total_price", #Total Return Index (SPXT) - for returns
                                          spot_col  = "price",   #Price Index (SPX) - to calculate fundamentals
                                          rf_col    = "rfr_1m",     # annualized rate in %
                                          vix_col   = "vix",
                                          m80_col   = "80%", m90_col = "90%", m100_col = "100%",
                                          m110_col  = "110%", m120_col = "120%",
                                          eps_col   = "eps_ttm",
                                          bvps_col  = "book",
                                          t3m_col   = "tbl_3m",
                                          lty_col   = "lty",
                                          infl_col  = "inflation_yoy",
                                          ts_col    = "term_spread",
                                          var1m_col = "stock_var_1m",  # realized ~21d variance column
                                          rf_basis  = 360) {           # 360 (USD/EUR/CHF), 365 (GBP/HKD/JPY/INR)
  eps   <- .Machine$double.eps
  m_vec <- c(0.8, 0.9, 1.0, 1.1, 1.2)
  
  D <- df_daily %>%
    dplyr::arrange(.data[[date_col]]) %>%
    dplyr::mutate(
      date      = as.Date(.data[[date_col]]),
      price     = .data[[price_col]],
      spot_price =.data[[spot_col]],
      logp      = log(price),
      
      # risk-free accrual to next trading day (assumes simple annualized yield)
      rf_ann    = zoo::na.locf(.data[[rf_col]] / 100, na.rm = FALSE),  # % → decimal
      next_date = dplyr::lead(date),
      dcf       = as.numeric(next_date - date) / rf_basis,
      rf_piece  = dplyr::if_else(is.na(dcf), NA_real_, log1p(rf_ann * dcf)),
      
      # carry low-frequency fundamentals and macro
      eps_ttm   = zoo::na.locf(.data[[eps_col]], na.rm = FALSE),
      book      = zoo::na.locf(.data[[bvps_col]], na.rm = FALSE),
      t3m_yield = zoo::na.locf(.data[[t3m_col]],   na.rm = FALSE),
      lty_yield = zoo::na.locf(.data[[lty_col]],   na.rm = FALSE),
      term_spread = zoo::na.locf(.data[[ts_col]],  na.rm = FALSE),
      infl_step   = snap_back_then_ffill(.data[[infl_col]], date, maxgap = 3),
      inflation_yoy = infl_step,
      
      # valuation ratios
      log_ep    = log(pmax(eps_ttm, eps) / pmax(spot_price, eps)),
      log_bm    = log(pmax(book,   eps) / pmax(spot_price, eps)),
      
      # IV smile, levels (decimals)
      iv_80   = .data[[m80_col]]  / 100,
      iv_90   = .data[[m90_col]]  / 100,
      iv_100  = .data[[m100_col]] / 100,
      iv_110  = .data[[m110_col]] / 100,
      iv_120  = .data[[m120_col]] / 100,
      vix_ann = .data[[vix_col]]  / 100,
      
      # realized ~21d variance
      var_1m  = .data[[var1m_col]],
      rv_1m   = sqrt(var_1m),
      
      # IV composites (diffs/ratios + simple wing stats)
      iv_skew_80_100 = iv_80 - iv_100,
      iv_skew_80_110 = iv_80 - iv_110,
      iv_skew_80_120 = iv_80 - iv_120,
      iv_skew_90_100 = iv_90 - iv_100,
      iv_skew_90_110 = iv_90 - iv_110,
      iv_skew_90_120 = iv_90 - iv_120,
      skew_ratio_80_100 = log(pmax(iv_80,  eps) / pmax(iv_100, eps)),
      skew_ratio_90_100 = log(pmax(iv_90,  eps) / pmax(iv_100, eps)),
      skew_ratio_80_110 = log(pmax(iv_80,  eps) / pmax(iv_110, eps)),
      skew_ratio_90_110 = log(pmax(iv_90,  eps) / pmax(iv_110, eps)),
      skew_ratio_80_120 = log(pmax(iv_80,  eps) / pmax(iv_120, eps)),
      skew_ratio_90_120 = log(pmax(iv_90,  eps) / pmax(iv_120, eps)),
      wing_slope        = (iv_120 - iv_80) / (0.4),
      wing_curve        = (iv_120 - 2*iv_100 + iv_80) / (0.2^2),
      
      # implied ~21d vol/var + premia (aligned basis with realized)
      iv_imp_1m  = vix_ann * sqrt(21/252),
      var_imp_1m = vix_ann^2 * (21/252),
      vol_prem_1m = iv_imp_1m  - rv_1m,
      var_prem_1m = var_imp_1m - var_1m
    )
  
  # ---------- quadratic smile fits ----------
  # A) In levels 
  D$coefs <- purrr::pmap(list(D$iv_80, D$iv_90, D$iv_100, D$iv_110, D$iv_120), ~{
    v <- c(..1, ..2, ..3, ..4, ..5)
    if (all(is.finite(v))) lm(v ~ m_vec + I(m_vec^2))$coefficients
    else c("(Intercept)"=NA_real_, "m_vec"=NA_real_, "I(m_vec^2)"=NA_real_)
  })
  D$b_quad <- purrr::map_dbl(D$coefs, ~ .x[["m_vec"]])
  D$c_quad <- purrr::map_dbl(D$coefs, ~ .x[["I(m_vec^2)"]])
  D$skew_slope_quad <- D$b_quad + 2*D$c_quad
  D$skew_curve_quad <- 2*D$c_quad
  
  # B) NEW: In logs — log(IV) ~ a + b * log(m) + c * [log(m)]^2
  mk <- log(m_vec)
  D$coefs_log <- purrr::pmap(list(D$iv_80, D$iv_90, D$iv_100, D$iv_110, D$iv_120), ~{
    v <- log(pmax(c(..1, ..2, ..3, ..4, ..5), eps))  # protect log
    if (all(is.finite(v))) lm(v ~ mk + I(mk^2))$coefficients
    else c("(Intercept)"=NA_real_, "mk"=NA_real_, "I(mk^2)"=NA_real_)
  })
  D$b_log <- purrr::map_dbl(D$coefs_log, ~ .x[["mk"]])
  D$c_log <- purrr::map_dbl(D$coefs_log, ~ .x[["I(mk^2)"]])
  D$log_slope_quad <- D$b_log              # derivative at log(m)=0 (ATM)
  D$log_curve_quad <- 2 * D$c_log          # curvature at ATM
  
  D %>% dplyr::select(-coefs, -b_quad, -c_quad, -coefs_log, -b_log, -c_log)
}


add_excess_returns_daily <- function(D, horizons_days = c(1,5,21,42,63,126,252)) {
  out <- D
  for (h in horizons_days) {
    # cumulative log rf from t..t+h-1, robust to gaps via slide_index
    cum_rf <- slider::slide_index_dbl(out$rf_piece, out$date, sum,
                                      .before = 0, .after = h - 1, .complete = TRUE)
    out[[paste0("excess_", h, "d")]] <- dplyr::lead(out$logp, h) - out$logp - cum_rf
  }
  out
}

# PC1 (level factor) from monthly smile (EDA-only)
add_pc1_to_daily <- function(D) {
  iv_cols <- c("iv_80","iv_90","iv_100","iv_110","iv_120")
  ok <- stats::complete.cases(D[, iv_cols, drop = FALSE])
  pc1 <- rep(NA_real_, nrow(D))
  if (sum(ok) >= length(iv_cols) + 2) {
    X  <- D[ok, iv_cols, drop = FALSE] |> scale() |> as.matrix()
    pc <- prcomp(X, center = FALSE, scale. = FALSE)
    pc1_raw <- as.numeric(pc$x[,1])
    sgn <- ifelse(stats::cor(pc1_raw, D$vix_ann[ok], use = "pairwise") < 0, -1, 1)
    pc1[ok] <- sgn * pc1_raw
  }
  dplyr::mutate(D, pc1 = pc1)
}

# ---- 4) Measure families (daily) --------------------------------------
iv_levels_d   <- c("iv_80","iv_90","iv_100","iv_110","iv_120","vix_ann")
realised_d    <- c("var_1m_mth","rv_1m_mth")
skew_diff_d   <- c("iv_skew_80_100","iv_skew_80_110","iv_skew_80_120",
                   "iv_skew_90_100","iv_skew_90_110","iv_skew_90_120")
skew_ratio_d  <- c("skew_ratio_80_100","skew_ratio_80_110","skew_ratio_80_120",
                   "skew_ratio_90_100","skew_ratio_90_110","skew_ratio_90_120")
wing_d        <- c("wing_slope","wing_curve","skew_slope_quad","skew_curve_quad","log_slope_quad","log_curve_quad")
implied_d     <- c("iv_imp_1m","var_imp_1m")
premia_d      <- c("vol_prem_1m","var_prem_1m")
valuation_macro_d <- c("log_ep","log_bm","t3m_yield","lty_yield","term_spread","inflation_yoy")

# ---- 5) Build daily panel (cache once) --------------------------------

D_full <- data %>%
  prepare_daily_for_forecasting() %>%
  add_excess_returns_daily(c(1,2,3,5,10,15,21,63,126,252)) %>%
  add_pc1_to_daily()

d_pct <- D_full %>% mutate(across(starts_with("iv_"), ~ . * 100),
                           vix_ann = vix_ann * 100)


# ---- 6) OOS + Clark West (daily) --------------------------------------

oos_stats_fast_daily <- function(
    D,
    horizons_days,
    window_days = 252,
    pred       = "iv_skew_80_100",
    benchmark  = c("rolling","expanding"),
    window_type= c("rolling","expanding")
) {
  benchmark   <- match.arg(benchmark)
  window_type <- match.arg(window_type)
  
  vapply(
    horizons_days,
    function(h) {
      y <- D[[paste0("excess_", h, "d")]]
      X <- as.matrix(D[, pred, drop = FALSE])
      
      ols <- get_ols_forecasts(
        y, X, h_days = h,   # argument name irrelevant; just pass h
        window_days = window_days,
        window_type   = window_type
      )
      f_mod <- ols$f_mod
      betas <- ols$betas
      
      f_bm <- if (benchmark == "rolling") {
        roll_mean_lag1(y, window_days)
      } else {
        exp_mean_lag1(y)
      }
      
      ok <- which(is.finite(f_mod) & is.finite(f_bm) & is.finite(y))
      if (length(ok) < 40)
        return(c(R2_raw=NA, cw_stat=NA, cw_p=NA, beta_avg=NA))
      
      y_act <- y[ok]
      err_m <- y_act - f_mod[ok]
      err_b <- y_act - f_bm[ok]
      
      R2_raw <- 100 * (1 - sum(err_m^2) / sum(err_b^2))
      
      adj_loss <- err_b^2 - (err_m^2 - (f_bm[ok] - f_mod[ok])^2)
      cw_var   <- sandwich::NeweyWest(lm(adj_loss ~ 1), lag = h, prewhite = FALSE)[1,1]
      cw_stat  <- mean(adj_loss) / sqrt(cw_var)
      cw_p     <- pnorm(cw_stat, lower.tail = FALSE)
      
      if (is.null(dim(betas))) {
        beta_avg <- mean(betas[ok], na.rm = TRUE)
      } else {
        beta_avg <- mean(rowMeans(betas[ok, , drop = FALSE], na.rm = TRUE), na.rm = TRUE)
      }
      
      c(R2_raw=R2_raw, cw_stat=cw_stat, cw_p=cw_p, beta_avg=beta_avg)
    },
    FUN.VALUE = c(R2_raw=NA_real_, cw_stat=NA_real_, 
                  cw_p=NA_real_, beta_avg=NA_real_)
  )
}

# ---- 7) OOS Master Compute of Results (daily) ------------------------

compute_oos_paths_raw_daily <- function(
    D,
    pred,
    h_days,
    window_days,
    benchmark   = c("rolling","expanding"),
    window_type = c("rolling","expanding"),
    min_points  = 12,
    eps         = 1e-12
) {
  benchmark   <- match.arg(benchmark)
  window_type <- match.arg(window_type)
  
  y <- D[[paste0("excess_", h_days, "d")]]
  X <- as.matrix(D[, pred, drop = FALSE])
  
  ols <- get_ols_forecasts(
    y          = y,
    X          = X,
    h_days     = h_days,
    window_days= window_days,
    window_type= window_type
  )
  f_mod <- ols$f_mod
  
  f_bm <- if (benchmark == "rolling") {
    roll_mean_lag1(y, window_days)
  } else {
    exp_mean_lag1(y)
  }
  
  ok_idx <- which(is.finite(y) & is.finite(f_mod) & is.finite(f_bm))
  if (length(ok_idx) < min_points) return(tibble())
  
  y_ok <- y[ok_idx]
  em   <- y_ok - f_mod[ok_idx]
  eb   <- y_ok - f_bm[ok_idx]
  
  cm <- cumsum(em^2); cb <- cumsum(eb^2)
  start <- min_points
  idx   <- start:length(ok_idx)
  R2cum <- 100 * (1 - cm[idx] / pmax(cb[idx], eps))
  
  tibble(
    date         = D$date[ok_idx][idx],
    cum_R2_raw   = R2cum,
    horizon      = paste0(h_days, "d"),
    window_days  = window_days,
    predictor    = pred,
    benchmark    = benchmark,
    window_type  = window_type
  )
}


compute_oos_suite_daily <- function(
    D,
    predictors,
    windows_days,
    horizons_days,
    benchmarks   = c("rolling","expanding"),
    window_type  = c("rolling","expanding"),
    min_points   = 12
) {
  benchmarks  <- match.arg(benchmarks, several.ok = TRUE)
  window_type <- match.arg(window_type)
  
  # Summary table
  results_table <- tidyr::crossing(
    predictor    = predictors,
    window_days  = windows_days,
    horizon      = horizons_days,
    benchmark    = benchmarks
  ) %>%
    dplyr::arrange(predictor, window_days, horizon, benchmark) %>%
    purrr::pmap_dfr(function(predictor, window_days, horizon, benchmark) {
      
      o <- oos_stats_fast_daily(
        D,
        horizons_days = horizon,
        window_days   = window_days,
        pred          = predictor,
        benchmark     = benchmark,
        window_type   = window_type
      )
      
      tibble::tibble(
        predictor    = predictor,
        window_days  = window_days,
        horizon      = paste0(horizon, "d"),
        benchmark    = benchmark,
        window_type  = window_type,
        R2_oos_raw   = round(o["R2_raw", ], 3),
        cw_stat      = as.numeric(o["cw_stat", ]),
        cw_p_raw     = as.numeric(o["cw_p", ]),   # <-- keep RAW p
        beta_avg     = round(o["beta_avg",], 6)
      )
    }) %>%
    # ---------- MULTIPLE TESTING FIX ----------
  # Within each horizon (and benchmark + window_type),
  # adjust across predictor × window_days
  dplyr::group_by(benchmark, window_type, horizon) %>%
    dplyr::mutate(cw_q_h = bh_adjust(cw_p_raw)) %>%
    dplyr::ungroup() %>%
    # Optional: global within benchmark+window_type across ALL horizons too
    dplyr::group_by(benchmark, window_type) %>%
    dplyr::mutate(cw_q_global_bench = bh_adjust(cw_p_raw)) %>%
    dplyr::ungroup() 
  
  # Cum R² paths (unchanged)
  paths_df <- tidyr::crossing(
    predictor    = predictors,
    window_days  = windows_days,
    horizon      = horizons_days,
    benchmark    = benchmarks
  ) %>%
    purrr::pmap_dfr(function(predictor, window_days, horizon, benchmark) {
      compute_oos_paths_raw_daily(
        D            = D,
        pred         = predictor,
        h_days       = horizon,
        window_days  = window_days,
        benchmark    = benchmark,
        window_type  = window_type,
        min_points   = min_points
      )
    }) %>%
    dplyr::mutate(
      window_lab = factor(
        paste0(window_days, "d Window"),
        levels = paste0(sort(unique(window_days)), "d Window")
      )
    )
  
  list(results_table = results_table, paths_df = paths_df)
}


# ---- 8) OOS Plots (daily) --------------------------------------------

plot_heatmap_raw_daily <- function(results, bm = NULL, alpha_txt = 0.85,
                                   sig_level = 0.10,
                                   p_col = c("cw_q_h", "cw_p_raw"),
                                   predictor_order = NULL) {
  
  p_col <- match.arg(p_col)
  
  df <- results %>%
    { if (!is.null(bm)) dplyr::filter(., benchmark == bm) else . } %>%
    dplyr::mutate(
      horizon_num  = readr::parse_number(horizon),
      window_days  = factor(window_days, levels = sort(unique(window_days))),
      predictor    = if (!is.null(predictor_order)) {
        factor(predictor, levels = predictor_order)
      } else {
        factor(predictor, levels = unique(predictor))
      },
      p_use        = .data[[p_col]],
      signif_flag  = is.finite(p_use) & (p_use < sig_level) & (R2_oos_raw > 0),
      R2_cap       = pmax(R2_oos_raw, -50)
    )
  
  win_lab <- labeller(window_days = function(x) paste0(x, "d Window"))
  
  gg <- ggplot(df, aes(x = factor(horizon_num, levels = sort(unique(horizon_num))),
                       y = predictor, fill = R2_cap)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label   = ifelse(is.finite(R2_oos_raw) & abs(R2_oos_raw) >= 0.1,
                                   sprintf("%.1f", R2_oos_raw), ""),
                  colour  = signif_flag,
                  fontface= ifelse(signif_flag, "bold", "plain")),
              size = 3.2, alpha = alpha_txt) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scico::scale_fill_scico(palette = "vikO",
                            limits  = range(df$R2_cap, na.rm = TRUE),
                            breaks  = scales::pretty_breaks(),
                            name    = expression(R[OOS]^2~"(%)")) +
    labs(x = "Forecast horizon (days)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.ticks  = element_blank(),
          legend.position = "bottom",
          legend.key.width = unit(2.5, "cm"),
          strip.text = element_text(face = "bold"))
  
  if (is.null(bm)) gg + facet_grid(benchmark ~ window_days, labeller = win_lab)
  else             gg + facet_grid(. ~ window_days, labeller = win_lab)
}


plot_oos_trend_through_time_raw_daily <- function(paths_df,
                                                  facet_by = c("window", "predictor")) {
  facet_by <- match.arg(facet_by)
  
  df <- paths_df %>%
    dplyr::mutate(
      horizon_num = readr::parse_number(horizon),
      horizon     = factor(
        horizon,
        levels = paste0(sort(unique(horizon_num)), "d")
      )
    )
  
  g <- ggplot(df, aes(date, cum_R2_raw, colour = horizon, group = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, alpha = 0.6) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    labs(
      y      = expression("Cumulative " * R[OOS]^2 * " (%)"),
      x      = NULL,
      colour = "Horizon"
    ) +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
  
  if (facet_by == "window") {
    g + facet_wrap(~ window_lab, ncol = 2, scales = "free_x")
  } else {
    g + facet_wrap(~ predictor, ncol = 2, scales = "free_x")
  }
}

# ---- 9) Build OOS  plots (daily) ----------------------

# 9.1 Expanding OOS result

pred_vars_exp_d  <- c("skew_ratio_90_120", "log_slope_quad")
windows_exp_d    <- c(1008, 1260, 1512, 2520)
horizons_exp_d   <- c(1,5,21,63,126,252)
benchmarks_exp_d <- "expanding"
window_type_exp  <- "expanding"

oos_all_exp_d <- compute_oos_suite_daily(
  D_full,
  predictors    = pred_vars_exp_d,
  windows_days  = windows_exp_d,
  horizons_days = horizons_exp_d,
  benchmarks    = benchmarks_exp_d,
  window_type   = window_type_exp
)

results_table_exp_d <- oos_all_exp_d$results_table
paths_df_exp_d     <- oos_all_exp_d$paths_df

# Expanding Heatmap
fig15_d <- plot_heatmap_raw_daily(
  results_table_exp_d,
  bm = benchmarks_exp_d,
  predictor_order = pred_vars_exp_d,
  sig_level = 0.10,
  p_col = "cw_q_h"     # <-- FDR-corrected significance
)

ggsave(
  filename = file.path(FIG_DAILY, "fig15_heatmap_oos_monthly_Expanding_Daily.png"),
  plot     = fig15_d,
  width    = 8, height = 4.5, dpi = 300
)

# 9.2 Rolling OOS result

pred_vars_roll_d  <- c("iv_skew_90_100", "log_slope_quad", "skew_ratio_80_110", "skew_ratio_80_120")
windows_roll_d    <- c(1008, 1260, 1512, 2520)
horizons_roll_d   <- c(1,5,21,63,126,252)
benchmarks_roll_d <- "rolling"
window_type_roll  <- "rolling"

oos_all_roll_d <- compute_oos_suite_daily(
  D_full,
  predictors    = pred_vars_roll_d,
  windows_days  = windows_roll_d,
  horizons_days = horizons_roll_d,
  benchmarks    = benchmarks_roll_d,
  window_type   = window_type_roll
)

results_table_roll_d <- oos_all_roll_d$results_table
paths_df_roll_d      <- oos_all_roll_d$paths_df

# Expanding Heatmap
fig16_d <- plot_heatmap_raw_daily(
  results_table_roll_d,
  bm = benchmarks_roll_d,
  predictor_order = pred_vars_roll_d,
  sig_level = 0.10,
  p_col = "cw_q_h"     # <-- FDR-corrected significance
)

ggsave(
  filename = file.path(FIG_DAILY, "fig16_heatmap_oos_monthly_Rolling_Daily.png"),
  plot     = fig16_d,
  width    = 8, height = 4.5, dpi = 300
)

S