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
       here
       )

here()

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
panel <- read_excel(here("Data", "EA-MD-QD", "EA-MD-QD-2026-02", 
                         "data_TR2", "eadataM_NA_TR2.xlsx"))

# Check start date and whether target variable is present
range(panel$Time)
"HICPNEF_EA" %in% names(panel)
"HICPOV_EA" %in% names(panel)
names(panel)