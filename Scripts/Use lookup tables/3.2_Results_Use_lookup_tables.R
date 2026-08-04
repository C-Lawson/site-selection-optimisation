# Load functions
source("Scripts/Use lookup tables/Functions/function_extract sequences.R") # This function is used to extract the sequences of recovery periods based on the disturbance regime.
source("Scripts/Use lookup tables/Functions/07.1_FUNCTION_Results_get_uncertainty_using_trajectory_lookup_tables.R") # This function runs the lookup tables to get the coral cover estimates based on the disturbance regime and starting coral cover.

# Use input data loaded by the master script.
if (!exists("coral_cover_data") ||
    !exists("cyclone_data") ||
    !exists("dhw_data") ||
    !exists("lookup_cyclone_data") ||
    !exists("lookup_dhw_data") ||
    !exists("lookup_repeats") ||
    !exists("lookup_last_year") ||
    !exists("cyclone_first_year") ||
    !exists("dhw_first_year")) {
  stop("Run Scripts/00_MASTER.R instead of this script.")
}
dat <- coral_cover_data
dat.cyclone <- cyclone_data
dat.dhw <- dhw_data

cyclone_last_year <- cyclone_first_year + ncol(dat.cyclone) - 1
dhw_last_year <- dhw_first_year + ncol(dat.dhw) - 1

lookup_cyclone <- lookup_cyclone_data
lookup_dhw <- lookup_dhw_data

# Parameters
repeats <- lookup_repeats # Number of times the simulation will be repeated to account for uncertainty in the results
last.year <- lookup_last_year # The final year for which coral-cover uncertainty is estimated.

if (last.year > cyclone_last_year || last.year > dhw_last_year) {
  stop(
    "lookup_last_year is not covered by both disturbance files. ",
    "Cyclone data cover ", cyclone_first_year, "-", cyclone_last_year,
    "; DHW data cover ", dhw_first_year, "-", dhw_last_year, "."
  )
}

# Store output summaries for all rows (including multiple rows per Reef.Number)
all_summaries <- vector("list", length = nrow(dat))

# Start for loop over all rows in dat
for (i in seq_len(nrow(dat))) {
  if (dat$Year[i] > last.year) {
    stop("Survey year is later than lookup_last_year for Reef.Number ", dat$Reef.Number[i], ".")
  }

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
  regime_years <- seq.int(this_dat$Year + 1, last.year)

  if (min(regime_years) < cyclone_first_year || min(regime_years) < dhw_first_year) {
    stop(
      "Disturbance data do not cover all years required for Reef.Number ", reef_num,
      ". Cyclone data begin in ", cyclone_first_year,
      " and DHW data begin in ", dhw_first_year, "."
    )
  }

  regime <- data.frame(
    Year = regime_years,
    Cyclone = as.numeric(dat.cyclone[reef_num, regime_years - cyclone_first_year + 1]),
    DHW = as.numeric(dat.dhw[reef_num, regime_years - dhw_first_year + 1])
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