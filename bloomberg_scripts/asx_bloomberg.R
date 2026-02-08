###############################################################################
# ASX-200 : 30-day implied-vol surface + macro / market panel
###############################################################################
library(dplyr)
library(tidyr)
library(purrr)
library(zoo)
library(lubridate)
library(Rblpapi)

blpConnect()

DATA_DIR <- "data"
dir.create(DATA_DIR, showWarnings = FALSE)

OUT_FILE <- file.path(DATA_DIR, "bbg_asx_data.rds")

# ── 1.  Define moneyness levels ──────────────────────────────────────────────
moneyness_levels <- c("MONEY_LVL_80_0","MONEY_LVL_90_0",
                      "MONEY_LVL_100_0","MONEY_LVL_110_0","MONEY_LVL_120_0")
moneyness_labels <- c("80%","90%","100%","110%","120%")
names(moneyness_labels) <- moneyness_levels

# ── 2.  Pull IV surface (underlier = ASX 200) ───────────────────────────────
get_iv_for_moneyness <- function(level_code){
  bdh(securities = "AS51 Index",              
      fields     = "IVOL_MONEYNESS",
      start.date = as.Date("2005-01-01"),
      end.date   = as.Date("2025-07-31"),
      overrides  = c(
        "IVOL_MATURITY"        = "MATURITY_30D",
        "IVOL_MONEYNESS_LEVEL" = level_code
      )) %>% 
    rename(date        = date,
           implied_vol = IVOL_MONEYNESS) %>% 
    mutate(moneyness = moneyness_labels[level_code])
}

iv_surface_df <- map_dfr(moneyness_levels, get_iv_for_moneyness)

iv_wide <- iv_surface_df %>% 
  select(date, moneyness, implied_vol) %>% 
  pivot_wider(names_from = moneyness, values_from = implied_vol)

# ── 3.  Market & macro data (AUS) ──────────────────────────────────────────
asx_tot  <- bdh("AS51T Index", "PX_LAST", as.Date("2005-01-01"), as.Date("2025-07-31")) %>% 
  rename(date = date, total_price = PX_LAST)

asx_spot <- bdh("AS51 Index",  "PX_LAST", as.Date("2005-01-01"), as.Date("2025-07-31")) %>% 
  rename(date = date, price = PX_LAST)

vix_df   <- bdh("AS51VIX Index", "PX_LAST", as.Date("2005-01-01"), as.Date("2025-07-31")) %>% 
  rename(date = date, vix = PX_LAST)

rfr_df   <- bdh("BBSW1M Index", "PX_LAST", as.Date("2005-01-01"), as.Date("2025-07-31")) %>% 
  rename(date = date, rfr_1m = PX_LAST)

tbl_df   <- bdh("BBSW3M Index", "PX_LAST", as.Date("2005-01-01"), as.Date("2025-07-31")) %>% 
  rename(date = date, tbl_3m = PX_LAST)

lty_df   <- bdh("GACGB10 Index","PX_LAST", as.Date("2005-01-01"), as.Date("2025-07-31")) %>% 
  rename(date = date, lty = PX_LAST)

book_df  <- bdh("AS51 Index", "BOOK_VAL_PER_SH",
                as.Date("2005-01-01"), as.Date("2025-07-31")) %>% 
  rename(date = date, book = BOOK_VAL_PER_SH)

infl_df  <- bdh("AUCPIYOY Index", "PX_LAST",
                as.Date("2005-01-01"), as.Date("2025-07-31")) %>% 
  rename(date = date, inflation_yoy = PX_LAST)

# ── 5.  SPX fundamentals (PE, EPS, PB) ─────────────────────────────────────
fund_df <- bdh("AS51 Index",
               c("PE_RATIO","TRAIL_12M_EPS","PX_TO_BOOK_RATIO"),
               as.Date("2005-01-01"), as.Date("2025-07-31")) %>% 
  rename(date     = date,
         pe_ratio = PE_RATIO,
         eps_ttm  = TRAIL_12M_EPS,
         pb_ratio = PX_TO_BOOK_RATIO)

# ── 6.  Merge everything ─────────────────────────────────────────────────────
master_df <- iv_wide %>% 
  left_join(asx_tot,  by = "date") %>% 
  left_join(asx_spot, by = "date") %>% 
  left_join(vix_df,   by = "date") %>% 
  left_join(rfr_df,   by = "date") %>% 
  left_join(tbl_df,   by = "date") %>% 
  left_join(lty_df,   by = "date") %>% 
  left_join(book_df,  by = "date") %>% 
  left_join(infl_df,  by = "date") %>% 
  left_join(fund_df,  by = "date") %>% 
  arrange(date) %>% 
  fill(pe_ratio, eps_ttm, pb_ratio, .direction = "down") %>% 
  mutate(
    earnings_yield = eps_ttm / price,
    term_spread    = lty - tbl_3m
  )

# ── 7.  Feature engineering ──────────────────────────────────────────────────
master_df <- master_df %>% 
  mutate(
    ret_daily        = log(total_price / lag(total_price)),
    return_1m_fwd    = lead(log(total_price), 21) - log(total_price),
    ret_1m_past      = rollapply(ret_daily, 21, sum, fill = NA, align = "right"),
    ret_3m_past      = rollapply(ret_daily, 63, sum, fill = NA, align = "right"),
    realized_vol_1m  = rollapply(ret_daily, 21, sd,  fill = NA, align = "right"),
    stock_var_1m     = rollapply(ret_daily^2, 21, sum, fill = NA, align = "right"),
    stock_var_12m    = rollapply(ret_daily^2,252, sum, fill = NA, align = "right")
  )

head(master_df)
saveRDS(master_df, OUT_FILE)
message("Saved: ", OUT_FILE)
