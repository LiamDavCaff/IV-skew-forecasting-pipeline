# Load packages
library(Rblpapi)
library(dplyr)
library(tidyr)
library(purrr)
library(zoo)
library(lubridate)

# Connect to Bloomberg
blpConnect()

DATA_DIR <- "data"
dir.create(DATA_DIR, showWarnings = FALSE)
OUT_FILE <- file.path(DATA_DIR, "spx_iv_surface_long.rds")

# --- 1. Define moneyness levels ---
moneyness_levels <- c("MONEY_LVL_80_0", "MONEY_LVL_90_0", "MONEY_LVL_100_0", 
                      "MONEY_LVL_110_0", "MONEY_LVL_120_0")
moneyness_labels <- c("80%", "90%", "100%", "110%", "120%")
names(moneyness_labels) <- moneyness_levels

# --- 2. Pull IV surface ---
get_iv_for_moneyness <- function(level_code) {
  iv_data <- bdh(
    securities = "SPX Index",
    fields = "IVOL_MONEYNESS",
    start.date = as.Date("2020-02-20"),
    end.date = as.Date("2020-03-10"),
    overrides = c(
      "IVOL_MATURITY" = "MATURITY_720D",
      "IVOL_MONEYNESS_LEVEL" = level_code
    )
  )
  iv_data$moneyness <- moneyness_labels[level_code]
  return(iv_data)
}

iv_surface_df <- map_dfr(moneyness_levels, get_iv_for_moneyness) %>%
  rename(date = date, implied_vol = IVOL_MONEYNESS)

iv_wide_720D <- iv_surface_df %>%
  pivot_wider(names_from = moneyness, values_from = implied_vol)

iv_wide_360D <- iv_surface_df %>%
  pivot_wider(names_from = moneyness, values_from = implied_vol)

iv_wide_180D <- iv_surface_df %>%
  pivot_wider(names_from = moneyness, values_from = implied_vol)

iv_wide_90D <- iv_surface_df %>%
  pivot_wider(names_from = moneyness, values_from = implied_vol)

iv_wide_60D <- iv_surface_df %>%
  pivot_wider(names_from = moneyness, values_from = implied_vol)

iv_wide_30D <- iv_surface_df %>%
  pivot_wider(names_from = moneyness, values_from = implied_vol)


iv_long <- bind_rows(
  iv_wide_30D  %>% mutate(tenor = "30D"),
  iv_wide_60D  %>% mutate(tenor = "60D"),
  iv_wide_90D  %>% mutate(tenor = "90D"),
  iv_wide_180D %>% mutate(tenor = "180D"),
  iv_wide_360D %>% mutate(tenor = "360D"),
  iv_wide_720D %>% mutate(tenor = "720D")
) %>% 
  relocate(tenor, .after = date)            

saveRDS(iv_long, file = OUT_FILE)
message("Saved: ", OUT_FILE)  

