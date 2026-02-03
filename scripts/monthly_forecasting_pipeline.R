# ===========================================================
#  Monthly Forecasting Pipeline (1m, 3m, 6m, 12m)
#  EDA → IS → IS|Level → Rolling Beta → OOS (Clark West)
#  - Monthly panel cached once
#  - Rolling / Expanding OLS via helper
#  - Benchmark via rolling / expanding mean
#  - Clark–West test
#  - NH adjustment for multiple testing
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


# ---- 0.1) Output folders ------------------------------------------------

OUT_FIG <- file.path("outputs", "monthly", "figures")
OUT_TAB <- file.path("outputs", "monthly", "tables")

FIG_EDA  <- file.path(OUT_FIG, "01_eda")
FIG_BETA <- file.path(OUT_FIG, "02_rolling_betas")
FIG_HM   <- file.path(OUT_FIG, "03_oos_heatmaps")
FIG_PATH <- file.path(OUT_FIG, "04_oos_paths")
FIG_APP  <- file.path(OUT_FIG, "05_appendix")

TAB_EDA  <- file.path(OUT_TAB, "01_eda")
TAB_IS   <- file.path(OUT_TAB, "02_insample")
TAB_APP  <- file.path(OUT_TAB, "03_appendix")

invisible(lapply(
  c(FIG_EDA, FIG_BETA, FIG_HM, FIG_PATH, FIG_APP, TAB_EDA, TAB_IS, TAB_APP),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

# ---- 1) Load data --------------------------------------------------------

# Update path if needed:
data <- readRDS("C:/Users/Liam/Documents/Dissertation/3. Code/2. Data/2. Data Sets for Modelling/spx_iv_modelling_data.rds")
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

# OLS forecasts for a given window type (rolling / expanding)
get_ols_forecasts <- function(
    y,
    X,
    h_months,
    window_months,
    window_type = c("rolling", "expanding"),
    min_obs     = window_months
) {
  window_type <- match.arg(window_type)
  n <- length(y)
  k <- ncol(X)
  
  f_mod <- rep(NA_real_, n)
  betas <- matrix(NA_real_, nrow = n, ncol = k)
  
  if (window_type == "rolling") {
    coef_mat <- roll::roll_lm(
      X, y,
      width        = window_months,
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

# ---- 3) Monthly panel (predictors + targets) ----------------------------

# Build month-end panel from daily data, including RF accrual timing
prepare_monthly_for_forecasting <- function(df_daily,
                                            date_col  = "date",
                                            price_col = "total_price", # Total Return Index (SPXT)
                                            spot_col  = "price",       # Price Index (SPX)
                                            rf_col    = "rfr_1m",      # annualized rate in %
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
  
  # Daily → with RF accrual to next calendar day
  Dd <- df_daily %>%
    dplyr::arrange(.data[[date_col]]) %>%
    dplyr::mutate(
      date        = as.Date(.data[[date_col]]),
      price       = .data[[price_col]],
      spot_price  = .data[[spot_col]],
      logp_d      = log(price),
      
      rf_ann      = zoo::na.locf(.data[[rf_col]] / 100, na.rm = FALSE),  # % → decimal
      next_date   = dplyr::lead(date),
      dcf         = as.numeric(next_date - date) / rf_basis,
      rf_piece    = dplyr::if_else(is.na(dcf), NA_real_, log1p(rf_ann * dcf)),
      
      # month tags
      month_start = lubridate::floor_date(date, "month"),
      month_end   = lubridate::floor_date(next_date, "month")  # attribute accrual to END month
    )
  
  # Risk-free accrued by END month (includes cross-month piece)
  rf_by_endmonth <- Dd %>%
    dplyr::group_by(month = month_end) %>%
    dplyr::summarise(
      rf_m = {
        x <- stats::na.omit(rf_piece)
        if (length(x) == 0) NA_real_ else sum(x)
      },
      .groups = "drop"
    )
  
  # Monthly features keyed by START month (last observation within month)
  M <- Dd %>%
    dplyr::mutate(
      eps_ttm     = zoo::na.locf(.data[[eps_col]], na.rm = FALSE),
      book        = zoo::na.locf(.data[[bvps_col]], na.rm = FALSE),
      t3m_yield   = zoo::na.locf(.data[[t3m_col]],   na.rm = FALSE),
      lty_yield   = zoo::na.locf(.data[[lty_col]],   na.rm = FALSE),
      term_spread = zoo::na.locf(.data[[ts_col]],    na.rm = FALSE),
      infl_step   = snap_back_then_ffill(.data[[infl_col]], date, maxgap = 3),
      
      # valuation ratios
      log_ep_d    = log(pmax(eps_ttm, eps) / pmax(spot_price, eps)),
      log_bm_d    = log(pmax(book,   eps) / pmax(spot_price, eps)),
      
      # IV smile, levels (decimals)
      iv_80_d     = .data[[m80_col]]  / 100,
      iv_90_d     = .data[[m90_col]]  / 100,
      iv_100_d    = .data[[m100_col]] / 100,
      iv_110_d    = .data[[m110_col]] / 100,
      iv_120_d    = .data[[m120_col]] / 100,
      vix_ann_d   = .data[[vix_col]]  / 100,
      
      # realized ~21d variance
      var_1m_d    = .data[[var1m_col]],
      rv_1m_d     = sqrt(var_1m_d)
    ) %>%
    dplyr::group_by(month = month_start) %>%
    dplyr::summarise(
      # last observation within month (EOM snapshot)
      iv_80        = last_non_na(iv_80_d),
      iv_90        = last_non_na(iv_90_d),
      iv_100       = last_non_na(iv_100_d),
      iv_110       = last_non_na(iv_110_d),
      iv_120       = last_non_na(iv_120_d),
      vix_ann      = last_non_na(vix_ann_d),
      var_1m_mth   = last_non_na(var_1m_d),
      rv_1m_mth    = last_non_na(rv_1m_d),
      
      price        = last_non_na(price),
      rfr_1m       = last_non_na(.data[[rf_col]]),   # original % reference
      logp         = log(price),
      log_ep       = last_non_na(log_ep_d),
      log_bm       = last_non_na(log_bm_d),
      t3m_yield    = last_non_na(t3m_yield),
      lty_yield    = last_non_na(lty_yield),
      term_spread  = last_non_na(term_spread),
      inflation_yoy= last_non_na(infl_step),
      .groups      = "drop"
    ) %>%
    dplyr::arrange(month) %>%
    dplyr::left_join(rf_by_endmonth, by = "month") %>%
    dplyr::mutate(
      # IV composites (diffs/ratios + wing)
      iv_skew_80_100   = iv_80 - iv_100,
      iv_skew_80_110   = iv_80 - iv_110,
      iv_skew_80_120   = iv_80 - iv_120,
      iv_skew_90_100   = iv_90 - iv_100,
      iv_skew_90_110   = iv_90 - iv_110,
      iv_skew_90_120   = iv_90 - iv_120,
      skew_ratio_80_100= log(pmax(iv_80,  eps) / pmax(iv_100, eps)),
      skew_ratio_90_100= log(pmax(iv_90,  eps) / pmax(iv_100, eps)),
      skew_ratio_80_110= log(pmax(iv_80,  eps) / pmax(iv_110, eps)),
      skew_ratio_90_110= log(pmax(iv_90,  eps) / pmax(iv_110, eps)),
      skew_ratio_80_120= log(pmax(iv_80,  eps) / pmax(iv_120, eps)),
      skew_ratio_90_120= log(pmax(iv_90,  eps) / pmax(iv_120, eps)),
      wing_slope       = (iv_120 - iv_80) / 0.4,
      wing_curve       = (iv_120 - 2*iv_100 + iv_80) / (0.2^2),
      
      # implied ~21d vol/var + premia (aligned basis with realized)
      iv_imp_1m   = vix_ann * sqrt(21/252),
      var_imp_1m  = vix_ann^2 * (21/252),
      vol_prem_1m = iv_imp_1m  - rv_1m_mth,
      var_prem_1m = var_imp_1m - var_1m_mth
    )
  
  # Quadratic smile fits (levels & logs)
  m_vec <- c(0.8, 0.9, 1.0, 1.1, 1.2); mk <- log(m_vec)
  
  M$coefs <- purrr::pmap(list(M$iv_80, M$iv_90, M$iv_100, M$iv_110, M$iv_120), ~{
    v <- c(..1, ..2, ..3, ..4, ..5)
    if (all(is.finite(v))) lm(v ~ m_vec + I(m_vec^2))$coefficients
    else c("(Intercept)"=NA_real_, "m_vec"=NA_real_, "I(m_vec^2)"=NA_real_)
  })
  M$b_quad <- purrr::map_dbl(M$coefs, ~ .x[["m_vec"]])
  M$c_quad <- purrr::map_dbl(M$coefs, ~ .x[["I(m_vec^2)"]])
  M$skew_slope_quad <- M$b_quad + 2*M$c_quad
  M$skew_curve_quad <- 2*M$c_quad
  
  M$coefs_log <- purrr::pmap(list(M$iv_80, M$iv_90, M$iv_100, M$iv_110, M$iv_120), ~{
    v <- log(pmax(c(..1, ..2, ..3, ..4, ..5), eps))
    if (all(is.finite(v))) lm(v ~ mk + I(mk^2))$coefficients
    else c("(Intercept)"=NA_real_, "mk"=NA_real_, "I(mk^2)"=NA_real_)
  })
  M$b_log <- purrr::map_dbl(M$coefs_log, ~ .x[["mk"]])
  M$c_log <- purrr::map_dbl(M$coefs_log, ~ .x[["I(mk^2)"]])
  M$log_slope_quad <- M$b_log
  M$log_curve_quad <- 2 * M$c_log
  
  M %>% dplyr::select(-coefs, -b_quad, -c_quad, -coefs_log, -b_log, -c_log)
}

# Excess returns with correct RF timing (sum rf_{m+1..m+h})
add_excess_returns_monthly <- function(M, horizons_months = c(1,3,6,12,24,36,60)) {
  out <- M
  rf_forward <- dplyr::lead(out$rf_m, 1)
  for (h in horizons_months) {
    cum_rf <- slider::slide_index_dbl(rf_forward, out$month, sum,
                                      .before = 0, .after = h - 1, .complete = TRUE)
    out[[paste0("excess_", h, "m")]] <- dplyr::lead(out$logp, h) - out$logp - cum_rf
  }
  out
}

# PC1 (level factor) from monthly smile (EDA-only)
add_pc1_to_monthly <- function(M) {
  iv_cols <- c("iv_80","iv_90","iv_100","iv_110","iv_120")
  ok <- stats::complete.cases(M[, iv_cols, drop = FALSE])
  pc1 <- rep(NA_real_, nrow(M))
  if (sum(ok) >= length(iv_cols) + 2) {
    X  <- M[ok, iv_cols, drop = FALSE] |> scale() |> as.matrix()
    pc <- prcomp(X, center = FALSE, scale. = FALSE)
    pc1_raw <- as.numeric(pc$x[,1])
    sgn <- ifelse(stats::cor(pc1_raw, M$vix_ann[ok], use = "pairwise") < 0, -1, 1)
    pc1[ok] <- sgn * pc1_raw
  }
  dplyr::mutate(M, pc1 = pc1)
}

# ---- 4) Measure families (monthly) --------------------------------------
iv_levels_m   <- c("iv_80","iv_90","iv_100","iv_110","iv_120","vix_ann")
realised_m    <- c("var_1m_mth","rv_1m_mth")
skew_diff_m   <- c("iv_skew_80_100","iv_skew_80_110","iv_skew_80_120",
                   "iv_skew_90_100","iv_skew_90_110","iv_skew_90_120")
skew_ratio_m  <- c("skew_ratio_80_100","skew_ratio_80_110","skew_ratio_80_120",
                   "skew_ratio_90_100","skew_ratio_90_110","skew_ratio_90_120")
wing_m        <- c("wing_slope","wing_curve","skew_slope_quad","skew_curve_quad","log_slope_quad","log_curve_quad")
implied_m     <- c("iv_imp_1m","var_imp_1m")
premia_m      <- c("vol_prem_1m","var_prem_1m")
valuation_macro_m <- c("log_ep","log_bm","t3m_yield","lty_yield","term_spread","inflation_yoy")

# ---- 4b) Global Config for Forecasting -----------------------------------

WINDOW_TYPE_BETA <- "rolling"   # "rolling" or "expanding" for rolling betas
WINDOW_TYPE_OOS  <- "rolling"   # "rolling" or "expanding" for OOS regression
BENCHMARK_OOS    <- "rolling"   # "rolling" or "expanding" for benchmark mean

# ---- 5) Build monthly panel (cache once) --------------------------------

M_full <- data %>%
  prepare_monthly_for_forecasting() %>%
  add_excess_returns_monthly(c(1,3,6,12)) %>%
  add_pc1_to_monthly() %>%
  arrange(month) %>%
  mutate(
    vix_ma3 = slider::slide_dbl(
      vix_ann,
      ~ mean(.x, na.rm = TRUE),
      .before = 2,
      .complete = FALSE
    )
  ) %>%
  { 
    p80 <- quantile(.$vix_ma3, probs = 0.80, na.rm = TRUE, names = FALSE)
    message("80th percentile of vix_ma3 = ", round(p80, 4))
    
    mutate(
      .,
      vix_p80 = p80,
      regime  = if_else(vix_ma3 >= p80, "crisis", "calm")
    )
  }


M_pct <- M_full %>% mutate(across(starts_with("iv_"), ~ . * 100),
                           vix_ann = vix_ann * 100)

# ---- 6) EDA (+ Appendix figs) ---------------------------------------------

# Table 1 (MAIN TEXT): Full-sample summary stats
vars_stats <- c(
  paste0("excess_", c(1,3,6,12), "m"),
  iv_levels_m, realised_m, skew_diff_m, skew_ratio_m, wing_m, implied_m, premia_m, valuation_macro_m
)

summ_one <- function(df, vars = vars_stats){
  vars <- intersect(vars, names(df))
  df %>%
    reframe(
      across(all_of(vars), \(x){
        x <- x[is.finite(x)]
        tibble(
          n    = length(x),
          mean = mean(x),
          sd   = stats::sd(x),
          p10  = stats::quantile(x, 0.10),
          med  = stats::median(x),
          p90  = stats::quantile(x, 0.90),
          min  = min(x),
          max  = max(x),
          ar1  = if (length(x) > 2) stats::acf(x, plot = FALSE, na.action = na.pass)$acf[2] else NA_real_
        )
      })
    ) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "stats") %>%
    unnest_wider(stats)
}

S_all <- summ_one(M_full) %>%
  mutate(variable = factor(variable, levels = vars_stats)) %>%
  arrange(variable)

# Print as HTML
S_all %>%
  mutate(across(c(mean, sd, p10, med, p90, min, max, ar1), ~ round(.x, 4))) %>%
  kable(format = "html", caption = "Full-sample summary statistics (monthly, end-of-month)") %>%
  kable_styling(full_width = FALSE)

# Pretty labels for LaTeX
var_labels_math <- c(
  "excess_1m"="Excess return (1m)", "excess_3m"="Excess return (3m)",
  "excess_6m"="Excess return (6m)", "excess_12m"="Excess return (12m)",
  "iv_80"="$\\mathrm{IV}_{80}$", "iv_90"="$\\mathrm{IV}_{90}$", "iv_100"="$\\mathrm{IV}_{100}$",
  "iv_110"="$\\mathrm{IV}_{110}$", "iv_120"="$\\mathrm{IV}_{120}$",
  "vix_ann"="$\\mathrm{VIX}_{\\text{ann}}$",
  "iv_skew_80_100"="$\\mathrm{iv\\_skew}_{80,100}$",
  "skew_ratio_80_100"="$\\mathrm{skew\\_ratio}_{80,100}$",
  "wing_slope"="wing\\_slope", "wing_curve"="wing\\_curve",
  "var_prem_1m"="var\\_prem\\_1m", "vol_prem_1m"="vol\\_prem\\_1m"
)

tab1 <- S_all %>%
  dplyr::mutate(variable = as.character(variable)) %>%
  dplyr::mutate(variable = dplyr::recode(variable, !!!var_labels_math)) %>%
  dplyr::mutate(dplyr::across(c(mean, sd, p10, med, p90, min, max, ar1), ~ round(.x, 3)))

col_heads_math <- c(
  "$\\text{Variable}$", "$n$", "$\\mu$", "$\\sigma$",
  "$Q_{10\\%}$", "$\\tilde{x}$", "$Q_{90\\%}$", "$\\min$", "$\\max$", "$\\phi_{1}$"
)

knitr::kable(
  tab1,
  format    = "latex",
  booktabs  = TRUE,
  longtable = FALSE,
  caption   = "Full-sample summary statistics (monthly, end-of-month).",
  label     = "tab:summary_stats",
  col.names = col_heads_math,
  align     = c("l", rep("r", ncol(tab1)-1)),
  escape    = FALSE
) |>
  kableExtra::kable_styling(latex_options = c("hold_position")) |>
  kableExtra::row_spec(0, align = "c") |>
  kableExtra::add_header_above(c(" " = 1, "Summary statistics" = ncol(tab1)-1)) |>
  kableExtra::footnote(
    general = "\\(\\mu\\)=mean, \\(\\sigma\\)=standard deviation, \\(Q_{p}\\)=quantile at probability \\(p\\), \\(\\tilde{x}\\)=median, \\(\\phi_{1}\\)=first-order autocorrelation.",
    threeparttable = TRUE, escape = FALSE
  ) |>
  kableExtra::save_kable("table1_summary_stats.tex")

# Figure 1: Mean IV by moneyness
iv_nodes <- c("iv_80","iv_90","iv_100","iv_110","iv_120")

iv_means <- S_all %>%
  filter(variable %in% iv_nodes) %>%
  mutate(
    variable  = factor(variable, levels = iv_nodes),
    moneyness = as.numeric(sub("iv_", "", as.character(variable)))
  )

fig1_main <- ggplot(iv_means, aes(x = moneyness, y = mean * 100)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = c(80,90,100,110,120),
                     labels = c("80%","90%","100%","110%","120%")) +
  labs(x = "Moneyness", y = "Mean IV (%)") +
  theme_minimal()

# Figure 2: Monthly IV by moneyness over time
iv_monthly_long <- M_pct %>%
  transmute(month, `80%`=iv_80,`90%`=iv_90,`100%`=iv_100,`110%`=iv_110,`120%`=iv_120) %>%
  pivot_longer(-month, names_to = "moneyness", values_to="iv") %>%
  mutate(moneyness = factor(moneyness, levels = c("80%","90%","100%","110%","120%")))

fig2_main <- ggplot(iv_monthly_long, aes(x = month, y = iv, colour = moneyness)) +
  geom_line(alpha = 0.65) +
  geom_smooth(se = FALSE, span = 0.2, linewidth = 0.6) +
  scale_y_continuous(name = "Implied Volatility (%)") +
  labs(x = NULL, colour = "Moneyness") +
  theme_minimal(base_size = 11)
print(fig2_main)

# Figure 3: skew vs 12m returns
df3 <- M_full %>%
  transmute(month,
            skew  = iv_skew_80_100,
            ret12 = excess_12m) %>%
  mutate(skew_z = as.numeric(scale(skew)),
         ret12_z = as.numeric(scale(ret12)),
         skew_z_lag1 = dplyr::lag(skew_z,1))

fig3_main <- ggplot(df3, aes(month)) +
  geom_line(aes(y = ret12_z, colour = "12m excess return"), linewidth = 0.6, na.rm = TRUE) +
  geom_line(aes(y = skew_z_lag1, colour = "iv_skew_80_100"), linewidth = 0.6, alpha = 0.9, na.rm = TRUE) +
  labs(x = NULL, y = "z-score") +
  scale_colour_manual(NULL, values = c("12m excess return" = "black",
                                       "iv_skew_80_100"     = "steelblue")) +
  theme_minimal(base_size = 11)

# Figure 4: Volatility feedback
Mm <- M_full %>%
  transmute(month, ret_12m = excess_12m,
            d_iv_pp = 100*(dplyr::lead(iv_100) - iv_100))

fig4_main <- ggplot(dplyr::filter(Mm, is.finite(ret_12m), is.finite(d_iv_pp)),
                    aes(ret_12m, d_iv_pp)) +
  geom_point(alpha = .25) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(x = "12-month excess return", y = "Next-month Δ ATM IV") +
  theme_minimal(base_size = 11)
print(fig4_main)

fit_m <- lm(d_iv_pp ~ ret_12m, data = Mm)
cat("\n[INFO] HAC slope (ΔIVpp ~ 12m return):\n")
print(lmtest::coeftest(fit_m, vcov = sandwich::NeweyWest(fit_m, lag = 12)))

# Figure 5a/b: PCA diagnostics
iv_cols <- c("iv_80","iv_90","iv_100","iv_110","iv_120")
ok_pca  <- stats::complete.cases(M_full[, iv_cols])
Xm      <- M_full[ok_pca, iv_cols, drop = FALSE] %>% scale() %>% as.matrix()
pcm     <- prcomp(Xm, center = FALSE, scale. = FALSE)
summary(pcm)

loadings_m <- as.data.frame(pcm$rotation) %>%
  tibble::rownames_to_column("moneyness") %>%
  mutate(moneyness = factor(moneyness, levels = c("iv_80","iv_90","iv_100","iv_110","iv_120")))
load_long_m <- loadings_m %>% select(moneyness, PC1, PC2, PC3) %>%
  pivot_longer(-moneyness, names_to="PC", values_to="loading")

fig5a_main <- ggplot(load_long_m, aes(moneyness, loading, group = PC, colour = PC)) +
  geom_line(aes(linetype = PC)) + geom_point(size = 2) +
  labs(x = "Moneyness node", y = "Loading") +
  theme_minimal(base_size = 11)
print(fig5a_main)

scores_m <- as_tibble(pcm$x) %>% mutate(month = M_full$month[ok_pca])
ov_m <- tibble(
  month = M_full$month[ok_pca],
  PC1_z = as.numeric(scale(scores_m$PC1)),
  VIX_z = as.numeric(scale(M_full$vix_ann[ok_pca]))
)

fig5b_main <- ggplot(ov_m, aes(month)) +
  geom_line(aes(y = PC1_z, colour = "PC1")) +
  geom_line(aes(y = VIX_z, colour = "VIX")) +
  labs(y = "z-score", x = NULL, colour = NULL) +
  theme_minimal()
print(fig5b_main)

cat("\n[INFO] PC1–VIX correlation (z): ",
    round(cor(scores_m$PC1, M_full$vix_ann[ok_pca], use="pairwise.complete.obs"), 3), "\n")

# Appendix Figure A1: smile fit R² over time
m_vec <- c(0.8,0.9,1.0,1.1,1.2)
m_fit_m <- M_full %>%
  rowwise() %>%
  mutate(
    fit_r2 = {
      y <- c(iv_80,iv_90,iv_100,iv_110,iv_120)
      if (any(!is.finite(y))) NA_real_ else summary(lm(y ~ poly(m_vec, 2)))$r.squared
    }
  ) %>% ungroup()

figA1_appx <- ggplot(m_fit_m, aes(month, fit_r2)) +
  geom_line() +
  labs(y="R² (quadratic skew fit)", x = NULL) +
  theme_minimal(base_size = 11)
print(figA1_appx)


# Appendix Figure A2: IV distribution by moneyness
figA2_appx <- ggplot(iv_monthly_long, aes(x = moneyness, y = iv, fill = moneyness)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "black") +
  scale_fill_scico_d(palette = "batlow", guide = "none") +
  labs(x = "Moneyness", y = "Implied Volatility (%)")  +
  theme_minimal(base_size = 11)
print(figA2_appx)
ggsave("Appendix_Dist_Moneyness.png", figA2_appx, width=8, height=4, dpi=300)

# Appendix Table A1: PCA variance explained
imp_m <- summary(pcm)$importance[2:3, 1:3]
imp_tbl <- as.data.frame(imp_m) %>%
  tibble::rownames_to_column("stat") %>%
  mutate(stat = dplyr::recode(
    stat,
    "Proportion of Variance" = "Proportion",
    "Cumulative Proportion"  = "Cumulative"
  )) %>%
  rename(PC1 = `PC1`, PC2 = `PC2`, PC3 = `PC3`)

kable(imp_tbl, format = "html",
      caption = "Table A1: PCA variance explained (monthly IV surface)") %>%
  kable_styling(full_width = FALSE)

# Optional exports
#ggsave("fig1_iv_skew_mean.png", fig1_main, width=8, height=4.5, dpi=300)
#ggsave("fig2_monthly_iv.png",  fig2_main, width=8, height=4.5, dpi=300)
#ggsave("fig3_skew_12m_returns.png", fig3_main, width=8, height=4.5, dpi=300)
#ggsave("fig4_vol_feedback.png",     fig4_main, width=8, height=4.5, dpi=300)
#ggsave("fig5a_pca_loadings.png",    fig5a_main, width=8, height=4.5, dpi=300)
#ggsave("fig5b_pc1_vix.png",         fig5b_main, width=8, height=4.5, dpi=300)
# ggsave("figA1_smile_fit.png", figA1_appx, width=7, height=3.8, dpi=300)
# ggsave("figA2_iv_violin.png",  figA2_appx, width=6.5, height=4.5, dpi=300)

# ---- 7) IN-SAMPLE (UNIVARIATE) ------------------------------------------

insample_stats_monthly <- function(M, horizons_months, pred = "iv_skew_80_100") {
  vapply(horizons_months, function(h) {
    y <- paste0("excess_", h, "m")
    ok <- is.finite(M[[pred]]) & is.finite(M[[y]])
    if (sum(ok) < 120/10) {
      return(c(R2_in=NA, adj_R2_in=NA, Beta=NA, Beta_1sd_y=NA, Beta_sd_sd=NA, t_HH=NA, p_HH=NA))
    }
    
    fit <- lm(reformulate(pred, y), data = M[ok, , drop = FALSE])
    R2        <- 100 * summary(fit)$r.squared
    adj_R2_in <- 100 * summary(fit)$adj.r.squared
    
    V   <- sandwich::NeweyWest(fit, lag = h, prewhite = FALSE, adjust = TRUE)
    se  <- sqrt(pmax(diag(V), 0))[pred]
    
    beta <- coef(fit)[pred]
    t_HH <- beta / se
    p_HH <- 2 * stats::pt(abs(t_HH), df = fit$df.residual, lower.tail = FALSE)
    
    sd_x <- sd(M[[pred]][ok]); sd_y <- sd(M[[y]][ok])
    beta_1sd_y <- beta * sd_x
    beta_sd_sd <- if (is.finite(sd_y) && sd_y > 0) beta_1sd_y / sd_y else NA_real_
    
    c(R2_in=R2, adj_R2_in=adj_R2_in,
      Beta=beta, Beta_1sd_y=beta_1sd_y, Beta_sd_sd=beta_sd_sd,
      t_HH=t_HH, p_HH=p_HH)
  },
  FUN.VALUE=c(R2_in=NA_real_, adj_R2_in=NA_real_, Beta=NA_real_, Beta_1sd_y=NA_real_,
              Beta_sd_sd=NA_real_, t_HH=NA_real_, p_HH=NA_real_))
}

in_sample_results_multi_from_M <- function(M,
                                           horizons_months  = c(1,3,6,12),
                                           pred_vars = c(skew_diff_m, skew_ratio_m, wing_m,
                                                         iv_levels_m, premia_m, valuation_macro_m),
                                           p_adjust_method = "BH") {
  
  out <- purrr::map_dfr(pred_vars, function(p) {
    mat <- insample_stats_monthly(M, horizons_months, pred = p)
    tibble::tibble(
      predictor   = p,
      horizon_m   = horizons_months,
      R2_in       = mat["R2_in", ],
      adj_R2_in   = mat["adj_R2_in",],
      Beta        = mat["Beta", ],
      Beta_1sd_y  = mat["Beta_1sd_y", ],
      Beta_sd_sd  = mat["Beta_sd_sd", ],
      t_HH        = mat["t_HH", ],
      p_HH        = mat["p_HH", ]
    )
  })
  
  # BH-FDR adjustment within each horizon across predictors
  out <- out %>%
    dplyr::group_by(horizon_m) %>%
    dplyr::mutate(
      p_FDR = {
        p <- p_HH
        outp <- rep(NA_real_, length(p))
        ok <- is.finite(p)
        if (any(ok)) outp[ok] <- p.adjust(p[ok], method = p_adjust_method)
        outp
      }
    ) %>%
    dplyr::ungroup()
  
  # Robustness: compute a global FDR across the whole grid
  out <- out %>%
    dplyr::mutate(
      p_FDR_global = {
        p <- p_HH
        outp <- rep(NA_real_, length(p))
        ok <- is.finite(p)
        if (any(ok)) outp[ok] <- p.adjust(p[ok], method = p_adjust_method)
        outp
      }
    )
  
  out
}

insample_tbl_m <- in_sample_results_multi_from_M(M_full, c(1,3,6,12)) %>%
  dplyr::mutate(
    adj_R2_in     = round(adj_R2_in,2),
    Beta          = round(Beta, 5),
    Beta_1sd_bps  = round(1e4 * Beta_1sd_y, 1),
    Beta_sd_sd    = round(Beta_sd_sd, 3),
    t_HH          = round(t_HH, 2),
    p_HH          = signif(p_HH, 3),
    p_FDR         = signif(p_FDR, 3),
    p_FDR_global  = signif(p_FDR_global, 3),
    horizon       = paste0(horizon_m, "m")
  ) 

print(insample_tbl_m, n = nrow(insample_tbl_m))

# LATEX rows for in-sample table
panel <- insample_tbl_m %>%
  select(predictor, horizon, adj_R2_in, p_FDR, Beta_1sd_bps) %>%
  mutate(
    Beta_1sd_bps = round(Beta_1sd_bps, 1),
    p_FDR        = round(p_FDR, 2),
    adj_R2_in    = round(adj_R2_in, 2)
  ) %>%
  pivot_wider(
    names_from  = horizon,
    values_from = c(Beta_1sd_bps, p_FDR, adj_R2_in),
    names_glue  = "{horizon}_{.value}"
  ) %>%
  arrange(predictor)

escape_pred <- function(x) {
  paste0("\\texttt{", gsub("_", "\\\\_", x), "}")
}

for (i in seq_len(nrow(panel))) {
  row <- panel[i, ]
  pred_ltx <- escape_pred(row$predictor)
  cat(
    pred_ltx, " & ",
    row$`1m_Beta_1sd_bps`, " & ", row$`1m_p_FDR`,  " & ", row$`1m_adj_R2_in`,  " & ",
    row$`3m_Beta_1sd_bps`, " & ", row$`3m_p_FDR`,  " & ", row$`3m_adj_R2_in`,  " & ",
    row$`6m_Beta_1sd_bps`, " & ", row$`6m_p_FDR`,  " & ", row$`6m_adj_R2_in`,  " & ",
    row$`12m_Beta_1sd_bps`," & ", row$`12m_p_FDR`, " & ", row$`12m_adj_R2_in`,
    " \\\\\n",
    sep = ""
  )
}

# ---- 8) IN-SAMPLE (CONDITIONAL on level: VIX) ---------------------------

insample_hh_multi_monthly <- function(M, h_months, main_pred, control = c("vix_ann")) {
  y_col  <- paste0("excess_", h_months, "m")
  x_vars <- c(main_pred, control)
  
  keep <- stats::complete.cases(M[, c(y_col, x_vars), drop = FALSE])
  if (sum(keep) < 12) {
    return(tibble::tibble(
      horizon = h_months, main_pred = main_pred, control = paste(control, collapse = "+"),
      dAIC = NA_real_, dBIC = NA_real_,
      t_main = NA_real_, p_main = NA_real_,
      beta_main_1sd_bps = NA_real_
    ))
  }
  
  df <- M[keep, , drop = FALSE]
  
  fit_full <- stats::lm(stats::reformulate(x_vars, response = y_col), data = df)
  fit_ctrl <- stats::lm(stats::reformulate(control, response = y_col), data = df)
  
  AIC_full <- stats::AIC(fit_full); AIC_ctrl <- stats::AIC(fit_ctrl)
  BIC_full <- stats::BIC(fit_full); BIC_ctrl <- stats::BIC(fit_ctrl)
  
  V_full <- sandwich::NeweyWest(fit_full, lag = h_months, prewhite = FALSE, adjust = TRUE)
  b      <- stats::coef(fit_full)
  se     <- sqrt(pmax(diag(V_full), 0))
  
  b_main  <- unname(b[main_pred])
  se_main <- unname(se[main_pred])
  
  t_main <- b_main / se_main
  p_main <- 2 * stats::pt(abs(t_main), df = fit_full$df.residual, lower.tail = FALSE)
  
  sd_x <- stats::sd(df[[main_pred]])
  beta_main_1sd_bps <- 1e4 * (b_main * sd_x)
  
  tibble::tibble(
    horizon = h_months,
    main_pred = main_pred,
    control = paste(control, collapse = "+"),
    dAIC = AIC_full - AIC_ctrl,
    dBIC = BIC_full - BIC_ctrl,
    t_main = t_main,
    p_main = p_main,
    beta_main_1sd_bps = beta_main_1sd_bps
  )
}


insample_hh_multi_batch_monthly <- function(M, main_preds, horizons_months = c(1,3,6,12), control = "vix_ann") {
  purrr::map_dfr(main_preds, \(p)
                 purrr::map_dfr(horizons_months, \(h) insample_hh_multi_monthly(M, h, p, control = control)))
}

main_preds_m <- c("iv_skew_80_100", "iv_skew_90_120", "skew_ratio_80_120", "skew_ratio_90_110", 
                  "skew_ratio_90_120", "log_slope_quad", "iv_110", "iv_120")

res_pc1_m <- insample_hh_multi_batch_monthly(M_full, main_preds_m, horizons_months = c(6,12), control = "vix_ann") |>
  dplyr::mutate(
    t_main = round(t_main, 2)
  )

print(res_pc1_m, n = nrow(res_pc1_m))

# ---- 9) Rolling/expanding betas (monthly) -----------------------------------------

roll_beta_monthly <- function(
    M,
    horizons_months = c(1, 3, 6, 12),
    window_months   = 120,
    pred            = "skew_ratio_90_100",
    window_type     = c("rolling", "expanding")
) {
  window_type <- match.arg(window_type)
  
  purrr::map_dfr(horizons_months, function(h) {
    
    y <- M[[paste0("excess_", h, "m")]]
    x <- M[[pred]]
    n <- length(y)
    
    beta_raw <- se_hac <- sd_x <- sd_y <- rep(NA_real_, n)
    
    # window indices helper
    get_idx <- function(i) {
      if (window_type == "rolling") {
        j <- i - window_months + 1
        if (j < 1) return(NULL)
        j:i
      } else {
        # uses your existing helper
        get_window_indices(i, window_months, window_type = "expanding")
      }
    }
    
    for (i in seq_len(n)) {
      idx <- get_idx(i)
      if (is.null(idx)) next
      
      yi <- y[idx]
      xi <- x[idx]
      ok <- is.finite(yi) & is.finite(xi)
      
      # keep your existing rule: require full usable window
      if (sum(ok) < window_months) next
      
      yi <- yi[ok]
      xi <- xi[ok]
      
      sd_x[i] <- stats::sd(xi)
      sd_y[i] <- stats::sd(yi)
      if (!is.finite(sd_x[i]) || sd_x[i] == 0) next
      if (!is.finite(sd_y[i]) || sd_y[i] == 0) next
      
      fit_i <- stats::lm(yi ~ xi)
      
      beta_raw[i] <- stats::coef(fit_i)[2]
      
      # Newey-West with lag = h (as you requested)
      V_i <- sandwich::NeweyWest(
        fit_i,
        lag      = h,
        prewhite = FALSE,
        adjust   = TRUE
      )
      se_hac[i] <- sqrt(pmax(diag(V_i)[2], 0))
    }
    
    # Standardised betas: per 1 SD move in predictor
    beta_sd <- beta_raw * sd_x
    se_sd   <- se_hac   * sd_x
    
    tibble::tibble(
      month       = M$month,
      regime      = M$regime,   # <<< bring regime through
      horizon     = paste0(h, "m"),
      pred        = pred,
      window      = window_months,
      window_type = window_type,
      beta_raw    = beta_raw,
      se_hac      = se_hac,
      sd_x        = sd_x,
      sd_y        = sd_y,
      beta_sd     = beta_sd,
      se_sd       = se_sd,
      t_hac       = beta_raw / se_hac
    )
  })
}

# Build regime spans for shading (contiguous crisis blocks

make_regime_spans <- function(M) {
  M %>%
    distinct(month, regime) %>%
    arrange(month) %>%
    mutate(
      regime = as.character(regime),
      chg = regime != lag(regime, default = first(regime)),
      grp = cumsum(chg)
    ) %>%
    group_by(grp, regime) %>%
    summarise(
      start = min(month),
      end   = max(month),
      .groups = "drop"
    ) %>%
    mutate(end = end %m+% months(1))  # cover the full last month
}

# rolling betas for selected predictors/windows/horizons

pred_for_beta_m <- c("iv_skew_90_100","log_slope_quad", "skew_ratio_80_110", "skew_ratio_80_120")
windows_beta_m  <- c(60)
horizons_beta_m <- c(3, 6, 12)

# WINDOW_TYPE_BETA is set to "rolling" or "expanding"
betas_all_m <- tidyr::crossing(pred = pred_for_beta_m, window = windows_beta_m) |>
  purrr::pmap_dfr(function(pred, window) {
    roll_beta_monthly(
      M_full,
      horizons_months = horizons_beta_m,
      window_months   = window,
      pred            = pred,
      window_type     = WINDOW_TYPE_BETA
    )
  }) |>
  dplyr::mutate(
    pred_label = stringr::str_replace_all(pred, "_", " "),
    horizon    = factor(horizon, levels = paste0(horizons_beta_m, "m")),
    window     = factor(window, levels = windows_beta_m,
                        labels = paste0(windows_beta_m, "-month window"))
  )

# regime spans (from the full panel)
regime_spans <- make_regime_spans(M_full)


# PLOT: beta_sd with crisis shading

fig_beta_sd_m <- betas_all_m |>
  ggplot2::ggplot(ggplot2::aes(month, beta_sd, colour = horizon, group = horizon)) +
  
  # shade crisis periods behind everything
  ggplot2::geom_rect(
    data = dplyr::filter(regime_spans, regime == "crisis"),
    ggplot2::aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.10
  ) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, alpha = 0.6) +
  ggplot2::geom_line(na.rm = TRUE, linewidth = 0.8) +
  ggplot2::facet_wrap(~ pred_label, ncol = 2, scales = "fixed") +
  ggplot2::labs(
    x = NULL,
    y = expression(beta~"(per 1 SD move in predictor)"),
    colour = "Horizon",
    fill   = "Horizon"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

fig_beta_sd_m
ggsave("fig_roll_beta_60m.png", fig_beta_sd_m, width=8, height=4, dpi=300)

# ---- 10) OOS + Clark West (monthly) --------------------------------------

oos_stats_fast_monthly <- function(
    M,
    horizons_months,
    window_months = 120,
    pred         = "iv_skew_80_100",
    benchmark    = c("rolling","expanding"),
    window_type  = c("rolling","expanding")
) {
  benchmark   <- match.arg(benchmark)
  window_type <- match.arg(window_type)
  
  vapply(
    horizons_months,
    function(h) {
      y <- M[[paste0("excess_", h, "m")]]
      X <- as.matrix(M[, pred, drop = FALSE])
      
      ols <- get_ols_forecasts(
        y, X, h_months = h,
        window_months = window_months,
        window_type   = window_type
      )
      f_mod <- ols$f_mod
      betas <- ols$betas
      
      f_bm <- if (benchmark == "rolling") {
        roll_mean_lag1(y, window_months)
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

# ---- 11) OOS Master Compute of Results (monthly) ------------------------

compute_oos_paths_raw_monthly <- function(
    M,
    pred,
    h_months,
    window_months,
    benchmark   = c("rolling","expanding"),
    window_type = c("rolling","expanding"),
    min_points  = 12,         # <- required number of forecasts
    eps         = 1e-12
) {
  benchmark   <- match.arg(benchmark)
  window_type <- match.arg(window_type)
  
  y <- M[[paste0("excess_", h_months, "m")]]
  X <- as.matrix(M[, pred, drop = FALSE])
  
  # Your helper for rolling/expanding OLS forecasts
  ols <- get_ols_forecasts(
    y,
    X,
    h_months      = h_months,
    window_months = window_months,
    window_type   = window_type
  )
  f_mod <- ols$f_mod
  
  # Benchmark forecasts
  f_bm <- if (benchmark == "rolling") {
    roll_mean_lag1(y, window_months)
  } else {
    exp_mean_lag1(y)
  }
  
  # All dates where y, model forecast and benchmark forecast exist
  ok_idx <- which(is.finite(y) & is.finite(f_mod) & is.finite(f_bm))
  if (length(ok_idx) < min_points) {
    # Not enough forecasts to start a cum-R2 path
    return(tibble())
  }
  
  y_ok <- y[ok_idx]
  em   <- y_ok - f_mod[ok_idx]
  eb   <- y_ok - f_bm[ok_idx]
  
  cm <- cumsum(em^2)
  cb <- cumsum(eb^2)
  
  # Cumulative R2 is computed from the start,
  # but we only *keep* points starting at `min_points`
  start <- min_points
  idx   <- start:length(ok_idx)
  R2cum <- 100 * (1 - cm[idx] / pmax(cb[idx], eps))
  
  tibble(
    month         = M$month[ok_idx][idx],
    cum_R2_raw    = R2cum,
    horizon       = paste0(h_months, "m"),
    window_months = window_months,
    predictor     = pred,
    benchmark     = benchmark,
    window_type   = window_type
  )
}

compute_oos_suite_monthly <- function(
    M,
    predictors,
    windows_months,
    horizons_months,
    benchmarks   = c("rolling","expanding"),
    window_type  = c("rolling","expanding"),
    min_points   = 12
) {
  benchmarks  <- match.arg(benchmarks, several.ok = TRUE)
  window_type <- match.arg(window_type)
  
  # ----- Summary table -----
  results_table <- tidyr::crossing(
    predictor     = predictors,
    window_months = windows_months,
    horizon       = horizons_months,
    benchmark     = benchmarks
  ) %>%
    dplyr::arrange(predictor, window_months, horizon, benchmark) %>%
    purrr::pmap_dfr(function(predictor, window_months, horizon, benchmark) {
      
      o <- oos_stats_fast_monthly(
        M,
        horizons_months = horizon,
        window_months   = window_months,
        pred            = predictor,
        benchmark       = benchmark,
        window_type     = window_type
      )
      
      tibble::tibble(
        predictor     = predictor,
        window_months = window_months,
        horizon       = paste0(horizon, "m"),
        benchmark     = benchmark,
        window_type   = window_type,
        R2_oos_raw    = round(o["R2_raw", ], 3),
        cw_stat       = as.numeric(o["cw_stat", ]),
        cw_p_raw      = as.numeric(o["cw_p", ]),      # <-- KEEP RAW
        beta_avg      = round(o["beta_avg",], 6)
      )
    }) %>%
    # ---------- MULTIPLE TESTING FIX ----------
  # Heatmap is read "within each horizon" while scanning predictors × windows.
  # So, for each (benchmark, window_type, horizon), adjust across ALL cells
  # in that panel: predictor × window_months.
  dplyr::group_by(benchmark, window_type, horizon) %>%
    dplyr::mutate(cw_q_h = bh_adjust(cw_p_raw)) %>%
    dplyr::ungroup() %>%
    # Robustness: global within benchmark+window_type across all horizons as well
    dplyr::group_by(benchmark, window_type) %>%
    dplyr::mutate(cw_q_global_bench = bh_adjust(cw_p_raw)) %>%
    dplyr::ungroup() 
  
  # ----- Cumulative R² paths (unchanged) -----
  paths_df <- tidyr::crossing(
    predictor     = predictors,
    window_months = windows_months,
    horizon       = horizons_months,
    benchmark     = benchmarks
  ) %>%
    purrr::pmap_dfr(function(predictor, window_months, horizon, benchmark) {
      compute_oos_paths_raw_monthly(
        M            = M,
        pred         = predictor,
        h_months     = horizon,
        window_months= window_months,
        benchmark    = benchmark,
        window_type  = window_type,
        min_points   = min_points
      )
    }) %>%
    dplyr::mutate(
      window_lab = factor(
        paste0(window_months, "m Window"),
        levels = paste0(sort(unique(window_months)), "m Window")
      )
    )
  
  list(
    results_table = results_table,
    paths_df      = paths_df
  )
}

# ---- 12) OOS Plots (monthly) --------------------------------------------

plot_heatmap_raw_monthly <- function(results, bm = NULL, alpha_txt = 0.85,
                                     sig_level = 0.10,
                                     p_col = c("cw_q_h", "cw_p_raw"),
                                     predictor_order = NULL) {
  
  p_col <- match.arg(p_col)
  
  df <- results %>%
    { if (!is.null(bm)) dplyr::filter(., benchmark == bm) else . } %>%
    dplyr::mutate(
      horizon_num   = readr::parse_number(horizon),
      window_months = factor(window_months, levels = sort(unique(window_months))),
      predictor     = if (!is.null(predictor_order)) {
        factor(predictor, levels = predictor_order)
      } else {
        factor(predictor, levels = unique(predictor))
      },
      p_use        = .data[[p_col]],
      signif_flag  = is.finite(p_use) & (p_use < sig_level) & (R2_oos_raw > 0),
      R2_cap       = pmax(R2_oos_raw, -50)
    )
  
  win_lab <- labeller(window_months = function(x) paste0(x, "m Window"))
  
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
    labs(x = "Forecast horizon (months)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.ticks  = element_blank(),
          legend.position = "bottom",
          legend.key.width = unit(2.5, "cm"),
          strip.text = element_text(face = "bold", size = 11))
  
  if (is.null(bm)) gg + facet_grid(benchmark ~ window_months, labeller = win_lab)
  else             gg + facet_grid(. ~ window_months, labeller = win_lab)
}



plot_oos_trend_by_horizon_raw_monthly <- function(results_table,
                                                  sig_level = 0.10,
                                                  p_col = c("cw_q_h", "cw_p_raw")) {
  p_col <- match.arg(p_col)
  
  df <- results_table %>%
    dplyr::mutate(
      horizon_num   = readr::parse_number(horizon),
      window_lab    = factor(
        window_months,
        levels = sort(unique(window_months)),
        labels = paste0(sort(unique(window_months)), "m")
      ),
      predictor_lab = stringr::str_replace_all(predictor, "_", " "),
      p_use         = .data[[p_col]],
      sig           = is.finite(p_use) & (p_use < sig_level) & (R2_oos_raw > 0)
    )
  
  ggplot(
    df,
    aes(
      x      = factor(horizon_num, levels = sort(unique(horizon_num))),
      y      = R2_oos_raw,
      group  = window_lab,
      colour = window_lab
    )
  ) +
    geom_line() +
    geom_point(aes(shape = sig), size = 2, na.rm = TRUE) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       guide = "none", na.translate = TRUE) +
    labs(
      x      = "Horizon (months)",
      y      = expression(R[OOS]^2~"(%, raw)"),
      colour = "Window",
      title  = "Raw OOS R² vs. horizon (facet = predictor) — Monthly"
    ) +
    facet_wrap(~ predictor_lab, scales = "free_y") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
}


plot_oos_trend_through_time_raw_monthly <- function(paths_df,
                                                    facet_by = c("window", "predictor")) {
  facet_by <- match.arg(facet_by)
  
  df <- paths_df %>%
    dplyr::mutate(
      horizon_num = readr::parse_number(horizon),
      horizon     = factor(horizon, levels = paste0(sort(unique(horizon_num)), "m"))
    )
  
  g <- ggplot(df, aes(x = .data$month, y = .data$cum_R2_raw,
                      colour = .data$horizon, group = .data$horizon)) +
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

# ---- 13) Build OOS suite + example plots (monthly) ----------------------

pred_vars_all_m <- c("iv_skew_90_100","log_slope_quad", "skew_ratio_80_110",  "skew_ratio_80_120")
windows_all_m   <- c(60)
horizons_all_m  <- c(3,6,12)
benchmarks_all_m <- BENCHMARK_OOS   # "rolling" or "expanding"

oos_all_m <- compute_oos_suite_monthly(
  M_full,
  predictors      = pred_vars_all_m,
  windows_months  = windows_all_m,
  horizons_months = horizons_all_m,
  benchmarks      = benchmarks_all_m,
  window_type     = WINDOW_TYPE_OOS
)

results_table_m <- oos_all_m$results_table
paths_df_m      <- oos_all_m$paths_df
print(results_table_m, n = nrow(results_table_m))
#Robust Variables
robust_tbl_m <- results_table_m %>%
  mutate(
    window_months = as.numeric(window_months),
    pass = (R2_oos_raw > 0) & is.finite(cw_q_h) & (cw_q_h <= 0.10)
  ) %>%
  group_by(window_type, benchmark, predictor, horizon) %>%
  arrange(window_months, .by_group = TRUE) %>%
  summarise(
    n_windows_pass = sum(pass, na.rm = TRUE),
    # Adjacent in the ordered window grid: 48->60, 60->72, 72->120
    robust_adj = any(pass & dplyr::lead(pass), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(robust_adj)

# Heatmap
pred_vars_all_heat <- c(iv_levels_m, skew_ratio_m, skew_diff_m, wing_m)

fig13_m <- plot_heatmap_raw_monthly(
  results_table_m,
  bm = BENCHMARK_OOS,
  predictor_order = pred_vars_all_heat,
  sig_level = 0.1,
  p_col = "cw_q_h"        # <-- IMPORTANT
)

#ggsave("fig1_heatmap_exp.png", fig13_m, width=8, height=6, dpi=300)

# Raw R² vs horizon
fig14_m <- plot_oos_trend_by_horizon_raw_monthly(
  results_table_m %>% dplyr::filter(benchmark == BENCHMARK_OOS)
)

# Cumulative raw R² through time for one key predictor
fig15_m <- plot_oos_trend_through_time_raw_monthly(
  paths_df_m %>% dplyr::filter(benchmark == BENCHMARK_OOS,
                               predictor %in% c("iv_skew_90_100","log_slope_quad", "skew_ratio_80_110",  "skew_ratio_80_120")),
  facet_by = "predictor"
)

ggsave("fig_roll_r2_60m.png", fig15_m, width=8, height=4, dpi=300)

# ---- 14) Forecast vs benchmark Plot (monthly) ---------------------------

plot_oos_vs_benchmark_monthly <- function(
    M,
    h_months     = 12,
    window_months = 120,
    pred_vars    = "skew_ratio_90_100",
    benchmark    = c("rolling","expanding"),
    window_type  = c("rolling","expanding")
) {
  benchmark   <- match.arg(benchmark)
  window_type <- match.arg(window_type)
  
  y_col <- paste0("excess_", h_months, "m")
  
  fc_df <- M %>%
    dplyr::select(month, actual = dplyr::all_of(y_col)) %>%
    dplyr::mutate(
      hist_mean = if (benchmark == "rolling")
        roll_mean_lag1(actual, window_months)
      else
        exp_mean_lag1(actual),
      forecast  = NA_real_
    )
  
  usable <- which(is.finite(fc_df$actual))
  for (i in usable) {
    idx <- get_window_indices(i, window_months, window_type = window_type)
    if (is.null(idx)) next
    if (length(idx) < window_months) next
    
    fit <- lm(stats::reformulate(pred_vars, y_col), data = M[idx, ])
    fc_df$forecast[i] <- predict(fit, newdata = M[i, , drop = FALSE])
  }
  
  fc_long <- tidyr::pivot_longer(
    fc_df, c(actual, hist_mean, forecast),
    names_to = "series", values_to = "value"
  ) %>%
    dplyr::mutate(series = factor(series, levels = c("actual", "hist_mean", "forecast")))
  
  bm_label <- if (benchmark == "rolling")
    paste0("Rolling mean (", window_months, "m)")
  else
    paste0("Expanding mean (", window_months, "m)")
  
  ggplot(fc_long, aes(x = month, y = value, colour = series)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    scale_colour_manual(
      name   = NULL,
      breaks = c("actual", "hist_mean", "forecast"),
      values = c(actual = "black", hist_mean = "steelblue", forecast = "firebrick"),
      labels = c("Actual", bm_label, paste(pred_vars, collapse = " + "))
    ) +
    labs(
      title = sprintf(
        "%s-Month Excess-Return Forecast — %s window %s vs %s benchmark",
        h_months, window_type, window_months, benchmark
      ),
      y = "Log excess return", x = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
}

# Example: keep your original illustrative plots
fig16m_roll <- plot_oos_vs_benchmark_monthly(
  M_full, h_months = 12, window_months = 12,
  pred_vars = "iv_skew_80_120",
  benchmark = "rolling",
  window_type = "rolling"
)

fig16m_exp  <- plot_oos_vs_benchmark_monthly(
  M_full, h_months = 12, window_months = 12,
  pred_vars = "skew_ratio_80_120",
  benchmark = "expanding",
  window_type = "expanding"
)

# ---- Final figure prints (optional) --------------------------------------
print(fig1_main); print(fig2_main); print(figA2_appx); print(fig3_main); print(fig4_main)
print(fig5a_main); print(fig5b_main)
purrr::walk(fig6_m, print); purrr::walk(fig7_m, print)
print(fig13_m); print(fig14_m); print(fig15_m)
print(fig16m_roll); print(fig16m_exp)
