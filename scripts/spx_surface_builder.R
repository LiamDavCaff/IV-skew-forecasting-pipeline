# ===========================================================
# Build an IV surface with 3 axes
# IV, Moneyness and Tenor (time to maturity)
# Using select day during COVID
# ===========================================================

# Load Packages
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(plot3D)
  library(viridisLite)
})

# ===========================
# user inputs
# ===========================
DATA_FILE   <- file.path("data", "spx_iv_surface_long.rds")
target_date <- as.Date("2020-03-03")

# output folder (monthly EDA)
OUT_DIR <- file.path("outputs", "monthly", "figures", "01_eda")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(
  OUT_DIR,
  sprintf("spx_iv_surface_3axes_%s.png", format(target_date, "%Y-%m-%d"))
)

# sanity check
if (!file.exists(DATA_FILE)) {
  stop(
    paste0(
      "Missing ", DATA_FILE, "\n",
      "Run the Bloomberg pull script to generate it first."
    ),
    call. = FALSE
  )
}

# (labels for months
label_tenor_mo <- function(t_mo) {
  sapply(t_mo, function(tm) {
    if (is.na(tm)) return(NA_character_)
    if (abs(tm - 1)  < 0.05) return("1M")
    if (abs(tm - 2)  < 0.05) return("2M")
    if (abs(tm - 3)  < 0.05) return("3M")
    if (abs(tm - 6)  < 0.10) return("6M")
    if (abs(tm - 12) < 0.20) return("1Y")
    if (abs(tm - 24) < 0.30) return("2Y")
    paste0(round(tm), "M")
  })
}

# ===========================
# load & tidy
# ===========================
iv_surface_spx <- readRDS(DATA_FILE)

iv_long <- iv_surface_spx %>%
  filter(date == target_date) %>%
  pivot_longer(
    cols = -c(date, tenor),
    names_to = "moneyness_lbl",
    values_to = "iv"
  ) %>%
  mutate(
    moneyness = as.numeric(str_remove(moneyness_lbl, "%"))/100,
    tenor_months = case_when(
      str_detect(tenor, "D$") ~ as.numeric(str_remove(tenor, "D$"))/30,
      str_detect(tenor, "M$") ~ as.numeric(str_remove(tenor, "M$")),
      str_detect(tenor, "Y$") ~ 12 * as.numeric(str_remove(tenor, "Y$")),
      TRUE ~ suppressWarnings(as.numeric(tenor))  # if already numeric
    )
  ) %>%
  select(tenor_months, moneyness, iv) %>%
  arrange(tenor_months, moneyness)

stopifnot(n_distinct(iv_long$tenor_months) >= 2, n_distinct(iv_long$moneyness) >= 2)

# ===========================
# grid -> matrix (rows = tenor, cols = moneyness)
# ===========================
tenors    <- sort(unique(iv_long$tenor_months))
moneys    <- sort(unique(iv_long$moneyness))

grid <- tidyr::complete(iv_long, tenor_months = tenors, moneyness = moneys) %>%
  arrange(tenor_months, moneyness)

Z <- matrix(grid$iv, nrow = length(tenors), ncol = length(moneys), byrow = TRUE)

# ===========================
# swap axes for plotting:
# want x = moneyness, y = tenor
# -> transpose Z so rows match length(x)
# ===========================
Zt <- t(Z)
stopifnot(length(moneys) == nrow(Zt), length(tenors) == ncol(Zt))

# ===========================
# plot & save
# ===========================
png(out_png, width = 1600, height = 1100, res = 150)
par(mar = c(4.8, 5.4, 3.2, 5.2))

persp3D(
  x = moneys,               # X = Moneyness (rows of Zt)
  y = tenors,               # Y = Tenor (cols of Zt)
  z = Zt,
  colvar = Zt,
  col = viridisLite::viridis(200),
  NAcol = "grey85",
  xlab = "Moneyness (K/S)",
  ylab = "Tenor (months)",
  zlab = "Implied Volatility (%)",
  ticktype = "detailed",
  phi = 20, theta = 45, expand = 0.6,
  bty = "b2", shade = 0.2, border = "grey40",
  colkey = list(clab = "IV", length = 0.5, width = 0.7, dist = 0.08)
)

mtext(paste0("SPX IV Surface — ", format(target_date, "%Y-%m-%d")),
      side = 3, line = 0.6, cex = 1)

dev.off()

cat("Saved:", normalizePath(out_png), "\n",
    "Moneyness (X):", paste(sprintf("%.2f", moneys), collapse = ", "), "\n",
    "Tenor months (Y):", paste(tenors, collapse = ", "), "\n")
