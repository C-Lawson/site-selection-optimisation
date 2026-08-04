# Run this file from the project folder:
# source("Scripts/00_MASTER.R")

# ---------------------------------------------------------------------------
# Input data files
# ---------------------------------------------------------------------------
data_files <- list(
  coral_cover = "Data/coral_cover.csv",
  cyclone = "Data/GBR_PAST_CYCLONES_2008-2026.csv",
  dhw = "Data/GBR_PAST_DHW_1985-2026.csv",
  lookup_cyclone = "Data/lookup_cyclone.RDS",
  lookup_dhw = "Data/lookup_dhw.RDS",
  manta_sem_model = "Data/2.01_AIMS_Manta_SEM_Model.RData",
  mmp_sem_model = "Data/2.02_MMP_SEM_Model.RData",
  transect_sem_model = "Data/2.03_AIMS transect_SEM_Model.RData",
  grc_sem_model = "Data/2.04_GRC_SEM_Model.RData",
  reefmod_sem_model = "Data/2.04.1_Reefmod_SEM_Model.RData",
  ksr_data = "Data/keysourcereefs2023.RData"
)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
# Number of uncertainty simulations for each reef.
lookup_repeats <- 10
# Final year for which coral-cover uncertainty is estimated.
lookup_last_year <- 2025
# First calendar year in the cyclone disturbance file.
cyclone_first_year <- 2008
# First calendar year in the DHW disturbance file.
dhw_first_year <- 1985

# Coral-cover percentage below which a reef is considered in need.
ksr_critic <- 20
# KSR weighting coefficient retained from the original calculation (advised to set at 1. See Mumby et al 2021 LOM).
ksr_alpha <- 1
# Number of processor cores used for KSR calculations.
ksr_parallel_cores <- max(1, parallel::detectCores() - 1)

# Percentage of top-ranked reefs considered within the KSR zone.
KSR_zone_threshold <- 5
# Number of reefs reported as potential KSR opportunities.
opp_number <- 100
# Number of reefs reported as at risk of leaving the KSR zone.
risk_number <- 100
# Plot colour for reefs at risk of leaving the KSR zone.
col_risk <- "#CC6677"
# Plot colour for reefs with an opportunity to enter the KSR zone.
col_opp <- "#44AA99"
# Plot colour for all other reefs.
col_default <- "#DDDDDD"

# Minimum percentage of a source's uncertainty that surveys must resolve.
optimisation_threshold <- 80
# Number of reefs selected for survey by the optimisation (this should match your survey capacity).
optimisation_reef_count <- 100
# Minimum coral-cover range required for a reef to be eligible for selection.
optimisation_min_coral_range <- 5
# Whether to display detailed optimisation-solver messages.
optimisation_verbose <- FALSE

# ---------------------------------------------------------------------------
# Load packages
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(foreach)
  library(doSNOW)
  library(ggimage)
  library(zoo)
  library(visreg)
  library(ompr)
  library(ompr.roi)
  library(ROI.plugin.highs)
})

# ---------------------------------------------------------------------------
# Read input data
# ---------------------------------------------------------------------------
missing_files <- unlist(data_files)[!file.exists(unlist(data_files))]
if (length(missing_files) > 0) {
  stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
}

coral_cover_data <- read.csv(data_files$coral_cover, header = TRUE, check.names = FALSE)
cyclone_data <- read.csv(data_files$cyclone, header = FALSE)
dhw_data <- read.csv(data_files$dhw, header = FALSE)
lookup_cyclone_data <- readRDS(data_files$lookup_cyclone)
lookup_dhw_data <- readRDS(data_files$lookup_dhw)

required_coral_columns <- c("Reef.Number", "Year", "Coral.Cover", "Method")
missing_coral_columns <- setdiff(required_coral_columns, names(coral_cover_data))
if (length(missing_coral_columns) > 0) {
  stop(
    "The coral-cover file is missing required column(s): ",
    paste(missing_coral_columns, collapse = ", ")
  )
}

load(data_files$manta_sem_model)
load(data_files$mmp_sem_model)
load(data_files$transect_sem_model)
load(data_files$grc_sem_model)
load(data_files$reefmod_sem_model)
load(data_files$ksr_data)

# ---------------------------------------------------------------------------
# Prepare output folders
# ---------------------------------------------------------------------------

dir.create("Outputs/Use lookup tables", recursive = TRUE, showWarnings = FALSE)
dir.create("Outputs/Key source reef uncertainty", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Estimate coral-cover uncertainty
# ---------------------------------------------------------------------------
source("Scripts/Use lookup tables/3.2_Results_Use_lookup_tables.R")

# ---------------------------------------------------------------------------
# 2. Calculate KSR score ranges
# ---------------------------------------------------------------------------
source("Scripts/Key source reef uncertainty/01_recoveryreef2_test_for_min_and_max.R")

# ---------------------------------------------------------------------------
# 3. Identify KSRs affected by uncertainty
# ---------------------------------------------------------------------------
source("Scripts/Key source reef uncertainty/02_results_get_uncertain_source_reefs.R")

# ---------------------------------------------------------------------------
# 4. Attribute uncertainty to individual reefs
# ---------------------------------------------------------------------------
source("Scripts/Key source reef uncertainty/03_results_flip_sink_values_ksmin.R")
source("Scripts/Key source reef uncertainty/03_1_results_flip_sink_values_ksmax.R")

# ---------------------------------------------------------------------------
# 5. Summarise the attribution results
# ---------------------------------------------------------------------------
source("Scripts/Key source reef uncertainty/04_results_identify_reefs_of_interest.R")

# ---------------------------------------------------------------------------
# 6. Select survey sites and save the results
# ---------------------------------------------------------------------------
source("Scripts/Key source reef uncertainty/05.2_optimiser_function.R")

optimisation_result <- run_ksr_optimization(
  T_thr = optimisation_threshold,
  K = optimisation_reef_count,
  MIN_CORAL_RANGE = optimisation_min_coral_range,
  verbose = optimisation_verbose
)

write.csv(
  optimisation_result$selected_reefs,
  file.path(
    "Outputs",
    "Key source reef uncertainty",
    paste0("o5.2_ILP_optimisation_selected_reefs_K", optimisation_reef_count, ".csv")
  ),
  row.names = FALSE
)

write.csv(
  optimisation_result$resolved_sources,
  file.path(
    "Outputs",
    "Key source reef uncertainty",
    paste0("o5.2_ILP_optimisation_resolved_reefs_K", optimisation_reef_count, ".csv")
  ),
  row.names = FALSE
)

cat("Workflow complete.\n")
