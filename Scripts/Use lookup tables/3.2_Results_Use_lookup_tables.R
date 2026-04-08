# Load libraries and functions
library(tidyverse)
library(ggimage)
library(zoo)
library(visreg)

source("Scripts/Use lookup tables/Functions/function_extract sequences.R") # This function is used to extract the sequences of recovery periods based on the disturbance regime.
source("Scripts/Use lookup tables/Functions/07.1_FUNCTION_Results_get_uncertainty_using_trajectory_lookup_tables.R") # This function runs the lookup tables to get the coral cover estimates based on the disturbance regime and starting coral cover.

# Load data.
dat <- read.csv("Data/coral_cover.csv", header = F) # coral cover data
dat.cyclone <- read.csv("Data/GBR_PAST_CYCLONES_2008-2024.csv", header = F) # cyclone data
dat.dhw <- read.csv("Data/GBR_PAST_DHW_1985-2024.csv", header = F) # dhw data

lookup_cyclone <- readRDS("Data/lookup_cyclone.RDS") #bigger than dhw because it contains "no disturbance" trajectories. 
lookup_dhw <- readRDS("Data/lookup_dhw.RDS")

# Parameters
repeats <- 100 # Number of times the simulation will be repeated to account for uncertainty in the results
last.year <- 2024 # The most recent year of survey and disturbance data. I.e. if you want uncertainty in coral cover for 2025, this should be 2024; if a reef has been surveyed in last.year then just use that observed value, but if a reef was last surveyed prior to last.year then run it through the lookup tables to get current uncertainty in coral cover. 

# Store output summaries for all rows (including multiple rows per Reef.Number)
all_summaries <- vector("list", length = nrow(dat))

# Start for loop over all rows in dat
for (i in seq_len(nrow(dat))) {
  if (dat$Year[i] == last.year) {
    all_summaries[[i]] <- NA  # Store NA if last survey was last year (i.e. nothing to simulate and we just want to use the most recent observed value) 
    next
  }
  
  print(i)
  this_dat <- dat[i, ]
  
  starting_cover <- this_dat$Coral.Cover  # Starting coral cover (i.e. coral cover at last survey). To simulate uncertainty in sampling protocol, this value can be used as a mean that is sampled from (+/- precision of sampling method) in each repeat run.
 
  Method <- this_dat$Method  # Survey method used to get last coral cover value.
  
  # Build disturbance regime to follow for this reef
  # Use Reef.Number to index cyclone/DHW matrices (not row index i, which may repeat for multiple sources per reef)
  reef_num <- this_dat$Reef.Number
  regime <- data.frame(
    Year = (this_dat$Year + 1):last.year,
    Cyclone = as.numeric(dat.cyclone[reef_num, (ncol(dat.cyclone) - (last.year - this_dat$Year) + 1):ncol(dat.cyclone)]),
    DHW = as.numeric(dat.dhw[reef_num, (ncol(dat.dhw) - (last.year - this_dat$Year) + 1):ncol(dat.dhw)])
  )
  
  # Run lookup tables
  output <- run_lookup_tables(
    regime = regime,
    starting_cover = starting_cover,
    Method = Method,
    repeats = repeats,
    lookup_cyclone = lookup_cyclone,
    lookup_dhw = lookup_dhw,
    extract_sequences = extract_sequences
  )
  
  all_summaries[[i]] <- output$summary
}


######
# Append lower and upper estimates to original file for each row.

# Create empty vectors to store the new columns
lower_estimate_vec <- numeric(nrow(dat))
upper_estimate_vec <- numeric(nrow(dat))
mean_vec <- numeric(nrow(dat))

# Loop over rows of dat and fill in values
for (i in seq_len(nrow(dat))) {
  if (dat$Year[i] == last.year) {
    # If the last survey was the final year, then just use the observed Coral.Cover value
    lower_estimate_vec[i] <- dat$Coral.Cover[i]
    upper_estimate_vec[i] <- dat$Coral.Cover[i]
    mean_vec[i] <- dat$Coral.Cover[i]
  } else {
    # Otherwise, use lookup summary uncertainty range of coral cover
    summary_i <- all_summaries[[i]]
    lower_estimate_vec[i] <- summary_i$lower_estimate
    upper_estimate_vec[i] <- summary_i$upper_estimate
    mean_vec[i] <- summary_i$mean_cover
  }
}

# Append new columns to dat
dat$Mean <- mean_vec
dat$CoverLower_estimate <- lower_estimate_vec
dat$CoverUpper_estimate <- upper_estimate_vec
dat$Range <- dat$CoverUpper_estimate - dat$CoverLower_estimate

# Select the row with the lowest Range for each Reef.Number to produce final 3806-row output
dat_final <- dat %>%
  group_by(Reef.Number) %>%
  slice_min(order_by = Range, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Reef.Number)


# Save to file (final output with 3806 rows - one per reef, selected by lowest Range)
write.csv(dat_final, "Outputs/Use lookup tables/3.2_lower_and_upper_coral_cover_estimates.csv", row.names = F)