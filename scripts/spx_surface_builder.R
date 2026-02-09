# ===========================================================
# Build an IV surface with 3 axes (IV, moneyness, tenor)
# ===========================================================

# ---- 0) Load Packages -----------------------------------------------------------

pkgs <- c("dplyr","tidyr","stringr","plot3D","viridisLite","readr","here")

to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1) Project root + output folders -----------------------------------------------------------

ROOT        <- here::here()  # repo root
OUT_DIR <- file.path(ROOT, "outputs", "monthly", "figures", "01_eda")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ---- 2) Load Data ------------------------------------------------------------

DATA_FILE   <- file.path(ROOT, "data", "spx_iv_surface_long.rds")
target_date <- as.Date("2020-03-03") #COVID

# sanity check
if (!file.exists(DATA_FILE)) {
  stop(
    paste0(
      "Missing ", DATA_FILE, "\n",
      "Generate it first (Bloomberg pull / surface builder)."
    ),
    call. = FALSE
  )
}

iv_surface_spx <- readRDS(DATA_FILE)

out_png <- file.path(
  OUT_DIR,
  sprintf("fig1_spx_iv_surface_%s.png", format(target_date, "%Y-%m-%d"))
)

# ---- 3) Helpers -----------------------------------------------------------------

# To convert tenor days to months
parse_tenor_months <- function(tenor) {
  t <- toupper(trimws(as.character(tenor)))
  num <- readr::parse_number(t)
  
  dplyr::case_when(
    is.na(t) ~ NA_real_,
    str_detect(t, "D$") ~ num / 30,
    str_detect(t, "W$") ~ (num * 7) / 30,
    str_detect(t, "M$") ~ num,
    str_detect(t, "Y$") ~ 12 * num,
    TRUE ~ num  # if already numeric-ish
  )
}

# ---- 4) Tidy data set  ----------------------------------------------------------

iv_long <- iv_surface_spx %>%
  filter(date == target_date) %>%
  pivot_longer(
    cols = -c(date, tenor),
    names_to = "moneyness_lbl",
    values_to = "iv"
  ) %>%
  mutate(
    moneyness    = as.numeric(str_remove(moneyness_lbl, "%")) / 100,
    tenor_months = parse_tenor_months(tenor)
  ) %>%
  select(tenor_months, moneyness, iv) %>%
  arrange(tenor_months, moneyness)

stopifnot(
  n_distinct(iv_long$tenor_months) >= 2,
  n_distinct(iv_long$moneyness)   >= 2
)


# ---- 5) Convert grid to matrix for plotting -------------------------------------

tenors <- sort(unique(iv_long$tenor_months))
moneys <- sort(unique(iv_long$moneyness))

grid <- tidyr::complete(iv_long, tenor_months = tenors, moneyness = moneys) %>%
  arrange(tenor_months, moneyness)

Z  <- matrix(grid$iv, nrow = length(tenors), ncol = length(moneys), byrow = TRUE)
Zt <- t(Z)  # x = moneyness, y = tenor
stopifnot(length(moneys) == nrow(Zt), length(tenors) == ncol(Zt))

# ---- 6) Plot and save surface ---------------------------------------------------

png(out_png, width = 1600, height = 1100, res = 150)
par(mar = c(4.8, 5.4, 3.2, 5.2))

persp3D(
  x = moneys,
  y = tenors,
  z = Zt,
  colvar = Zt,
  col = viridisLite::viridis(200),
  NAcol = "grey85",
  xlab = "Moneyness (K/S)",
  ylab = "Tenor (months)",
  zlab = "Implied Volatility",
  ticktype = "detailed",
  phi = 20, theta = 45, expand = 0.6,
  bty = "b2", shade = 0.2, border = "grey40",
  colkey = list(length = 0.5, width = 0.7, dist = 0.08) # removed clab to avoid warnings
)

dev.off()

cat("Saved:", out_png, "\n")
