pacman::p_load(pacman,
       ecb,
       eurostat,
       fredr,
       readr,
       dplyr,
       lubridate,
       tidyr,
       httr,
       jsonlite,
       readxl,
       here,
       httr2,
       purrr,
       zoo,
       eurostat,
       fredr
       )



# Download EA-MD-QD (February 2026 version)
url      <- "https://zenodo.org/records/18804061/files/EA-MD-QD-2026-02.zip?download=1"
destfile <- here("Data", "EA-MD-QD", "EA-MD-QD-2026-02.zip")

download.file(url, destfile, mode = "wb")

# Unzip into a folder
unzip(destfile, exdir = here("Data", "EA-MD-QD"))

list.files(here("Data", "EA-MD-QD"), recursive = TRUE)

# Load the raw EA aggregate data
ea_data <- read_excel(here("Data", "EA-MD-QD", "EA-MD-QD-2026-02", "EAdata.xlsx"))

# Check structure
head(ea_data)
dim(ea_data)
names(ea_data)


# Ran the MATLAB routine_data.m for data transformation:
# Settings: country = EA, frequency = M, transformation = light, imputation = none
# Output saved to: Data/EA-MD-QD/EA-MD-QD-2026-02/data_TR2/eadataM_NA_TR2.xlsx

# Load the preprocessed monthly panel
panel_a <- read_excel(here("Data", "EA-MD-QD", "EA-MD-QD-2026-02", 
                         "data_TR2", "eadataM_NA_TR2.xlsx"))

# Check start date and whether target variable is present
range(panel_a$Time)
"HICPNEF_EA" %in% names(panel_a)
"HICPOV_EA" %in% names(panel_a)
names(panel_a)
saveRDS(panel_a, here("Data", "panel_a.rds"))
# ──────────────────────────────────────────────────────────────────────────────
# Adding ECB Data sources separately

ECB_BASE <- "https://data-api.ecb.europa.eu/service/data"

# ── Helper ────────────────────────────────────────────────────────────────────
ecb_get <- function(dataflow, key, name, start = "1997-01") {
  url <- paste(ECB_BASE, dataflow, key, sep = "/")
  
  resp <- tryCatch(
    request(url) %>%
      req_url_query(format = "csvdata", startPeriod = start) %>%
      req_error(is_error = \(r) FALSE) %>%
      req_perform(),
    error = function(e) { message("FAIL [", name, "] ", e$message); NULL }
  )
  
  if (is.null(resp)) return(NULL)
  
  if (resp_status(resp) != 200) {
    message("FAIL [", name, "] HTTP ", resp_status(resp))
    return(NULL)
  }
  
  resp |>
    resp_body_string() %>%
    read_csv(show_col_types = FALSE) %>%
    select(date = TIME_PERIOD, value = OBS_VALUE) %>%
    mutate(value = as.numeric(value), series = name) %>%
    filter(!is.na(value))
}

# ── B1. Credit aggregates ─────────────────────────────────────────────────────
loans_nfc <- ecb_get("BSI", "M.U2.N.A.A20.A.1.U2.2240.Z01.E",
                     "Loans_NFC",   start = "1997-09")

loans_hh  <- ecb_get("BSI", "M.U2.N.A.A20.A.1.U2.2250.Z01.E",
                     "Loans_HH",    start = "1997-09")

loans_pvt <- ecb_get("BSI", "M.U2.N.A.A20.A.1.U2.2200.Z01.E",
                     "Loans_Total", start = "1997-09")

# ── B2. Macroeconomics ────────────────────────────────────────────────────────
m3        <- ecb_get("BSI", "M.U2.Y.V.M30.X.1.U2.2300.Z01.E",
                     "M3",          start = "1980-01")

neg_wages <- ecb_get("STS", "Q.U2.N.INWR.000000.3.ANR",
                     "Neg_Wages",   start = "1991-01")

# ── B3. Inflation expectations (SPF) ─────────────────────────────────────────
spf_1y <- ecb_get("SPF", "M.U2.HICP.POINT.P12M.Q.AVG",
                  "SPF_1Y", start = "1999-12")

spf_2y <- ecb_get("SPF", "M.U2.HICP.POINT.P24M.Q.AVG",
                  "SPF_2Y", start = "2000-12")

spf_lt    <- ecb_get("SPF", "Q.U2.HICP.POINT.LT.Q.AVG",
                     "SPF_LT",      start = "1999-Q1")

# ── B4. Global sentiment ──────────────────────────────────────────────────────
ciss      <- ecb_get("CISS", "M.U2.Z0Z.4F.EC.SOV_EW.IDX",
                     "CISS",        start = "1999-01")

# ── B5. Price passthrough ─────────────────────────────────────────────────────
nec       <- ecb_get("NEC", "M.I9.N.ECPE.CTOTNE.4F0.N.IX",
                     "NEC",         start = "1996-01")


# ──────────────────────────────────────────────────────────────────────────────

# ── 1. Wu-Xia European Shadow Rate (Google Drive) ────────────────────────────
# Note: data runs to April 2021 — extend with ECB policy rate post-2021 if needed
shadow_rate <- read_excel(here("Data", "shadowrate_ECB.xls"), 
                          col_names = FALSE, col_types = "text") %>%
  setNames(c("date", "shadow_rate")) %>%
  mutate(
    date        = ymd(paste0(as.integer(date), "01")),
    shadow_rate = as.numeric(shadow_rate)
  ) |>
  filter(!is.na(shadow_rate), !is.na(date))


# ── 2. BCS Capacity Utilisation — Quarterly ─────────
cap_util_raw <- get_eurostat("ei_bsin_q_r2", time_format = "date")

cap_util <- cap_util_raw %>%
  filter(
    geo   == "EA20",
    indic == "BS-ICU-PC",
    s_adj == "SA"
  ) %>%
  arrange(TIME_PERIOD) %>%
  transmute(
    date     = TIME_PERIOD,
    cap_util = values
  )


# ── 3. BCS Short-Term Consumer Inflation Expectations — Monthly ───────────────

infl_exp_st_raw <- get_eurostat("ei_bsco_m", time_format = "date")

infl_exp_st <- infl_exp_st_raw %>%
  filter(
    geo   == "EA20",
    indic == "BS-PT-NY",
    s_adj == "SA"
  ) %>%
  arrange(TIME_PERIOD) %>%
  transmute(
    date        = TIME_PERIOD,
    infl_exp_st = values
  )


# ── 4. NY Fed GSCPI — Monthly ─────────────────────────────────────────────────
# Copy with correct extension and read

file.copy(here("Data", "gscpi_data.xls"), here("Data", "gscpi_data.xlsx"))

gscpi <- read_excel(here("Data", "gscpi_data.xlsx"), 
                    sheet = "GSCPI Monthly Data",
                    skip = 5, col_names = FALSE) %>%
  select(1, 2) %>%
  setNames(c("date", "gscpi")) %>%
  mutate(
    date  = floor_date(dmy(date), "month"),   # "31-Jan-1998" → 1998-01-01
    gscpi = as.numeric(gscpi)
  ) %>%
  filter(!is.na(gscpi), !is.na(date))


# ── 5. Eurostat Import Prices in Industry — Monthly ───────────────────────────
import_prices_raw <- get_eurostat("sts_inpi_m", time_format = "date")
names(import_prices_raw)

import_prices <- import_prices_raw %>%
  filter(
    geo      == "EA20",
    indic_bt == "PRC_IMP_NEU",
    cpa2_1   == "CPA_B_C_X_MIG_NRG", # mining and manufacturing excluding energy
    unit     == "I15", # index 2015=100
    s_adj    == "NSA"
  ) %>%
  arrange(TIME_PERIOD) %>%
  transmute(
    date          = TIME_PERIOD,
    import_prices = values
  ) %>%
  filter(!is.na(import_prices))


# ── 6. Brent Crude Oil (FRED: DCOILBRENTEU) — Daily → Monthly mean ────────────
fredr_set_key("c0ba7ce01b21d7020a10c57b9578e589")

brent <- fredr(
  series_id        = "DCOILBRENTEU",
  observation_start = as.Date("1987-01-01"),
  frequency        = "m",
  aggregation_method = "avg"
) %>%
  transmute(
    date  = date,
    brent = value
  ) %>%
  filter(!is.na(brent))

# ── 7. Global Price of Natural Gas (FRED: PNGASEUUSDM) — Daily → Monthly mean ────────────
ttf_gas <- fredr(
  series_id         = "PNGASEUUSDM",      # EU natural gas price, monthly
  observation_start = as.Date("1999-01-01")
) %>%
  transmute(
    date    = date,
    ttf_gas = value
  ) %>%
  filter(!is.na(ttf_gas))





# ── Monthly ECB series → wide ─────────────────────────────────────────────────
panel_ecb_monthly <- list(loans_nfc, loans_hh, loans_pvt, m3, ciss, nec) %>%
  discard(is.null) %>%
  bind_rows() %>%
  mutate(date = ym(date)) %>%
  summarise(value = mean(value, na.rm = TRUE), .by = c(date, series)) %>%
  pivot_wider(names_from = series, values_from = value)

saveRDS(panel_ecb_monthly, here("Data", "panel_ecb_monthly.rds"))

# ── Smart ECB date parser (handles both "YYYY-MM" and "YYYY-Qq") ─────────────
parse_ecb_date <- function(x) {
  x      <- as.character(x)
  is_qtr <- grepl("[Qq]\\d", x)
  result <- as.Date(rep(NA, length(x)))
  
  # "YYYY-Qq"  →  first month of that quarter
  if (any(is_qtr)) {
    year <- as.integer(regmatches(x[is_qtr], regexpr("\\d{4}", x[is_qtr])))
    q    <- as.integer(regmatches(x[is_qtr], regexpr("(?<=[Qq])\\d", x[is_qtr], perl = TRUE)))
    result[is_qtr] <- as.Date(paste0(year, "-", (q - 1L) * 3L + 1L, "-01"))
  }
  
  # "YYYY-MM"  →  first day of that month
  if (any(!is_qtr)) {
    result[!is_qtr] <- as.Date(paste0(x[!is_qtr], "-01"))
  }
  result
}

# ── Quarterly ECB series → wide (without forward-filling) ─────────────────────────
panel_ecb_quarterly <- list(spf_1y, spf_2y, spf_lt, neg_wages) %>%
  discard(is.null) %>%
  bind_rows() %>%
  mutate(date = parse_ecb_date(date)) %>%
  summarise(value = mean(value, na.rm = TRUE), .by = c(date, series)) %>%
  pivot_wider(names_from = series, values_from = value)

saveRDS(panel_ecb_quarterly, here("Data", "panel_ecb_quarterly.rds"))

# ── Join everything ───────────────────────────────────────────────────────────
panel_b_full <- panel_ecb_monthly %>%
  left_join(panel_ecb_quarterly, by = "date") %>%
  left_join(cap_util,            by = "date") %>%
  left_join(shadow_rate,         by = "date") %>%
  left_join(infl_exp_st,         by = "date") %>%
  left_join(gscpi,               by = "date") %>%
  left_join(import_prices,       by = "date") %>%
  left_join(brent,               by = "date") %>%
  left_join(ttf_gas,             by = "date") %>%
  filter(date >= as.Date("1999-01-01")) %>%
  arrange(date) 

glimpse(panel_b_full)

panel_b_filled <- panel_b_full %>%
  fill(SPF_1Y, SPF_2Y, SPF_LT, Neg_Wages, cap_util, .direction = "down")


saveRDS(panel_b_full, here("Data", "panel_b_full.rds"))
saveRDS(panel_b_filled, here("Data", "panel_ecb_mfilled.rds"))

