# Aim: Get the range of uncertainty for a given a) starting coral cover (i.e. last known coral cover) and b) disturbance regime since last known coral cover including recovery periods following each disturbance. 

# Note that coral cover values are given for the END of the associated year and at the end of each sampled trajectory.

# Use SEM precision models loaded by the master script.
if (!exists("model.manta") ||
    !exists("model.mmp") ||
    !exists("model.ltmp") ||
    !exists("model.grc") ||
    !exists("model.reefmod")) {
  stop("Run Scripts/00_MASTER.R instead of this script.")
}


##################################################################################################################################################################################################################################
run_lookup_tables <- function(regime, starting_cover, Method, repeats,
                                      lookup_cyclone, lookup_dhw,
                                      extract_sequences) {

# Correct any years that have both cyclones and DHW by setting DHW to 0.
regime$DHW[regime$Cyclone > 0] <- 0

# First convert exact DHW to DHW categories to match the lookup tables.
regime$DHW_cat <- cut(regime$DHW,
                      breaks = c(-Inf, 3, 7, 11, 15, 19, Inf),
                      labels = c(0, 3, 7, 11, 15, 19),
                      right = FALSE)

# Get sequence of disturbance periods (i.e. sequence of trajectories to sample).
sequences_raw <- extract_sequences(regime)


###############################################################################################################################################
# Run the lookup tables, modifying the starting coral cover based on parameter above and then subsequent coral values for next trajectory periods based on whatever the cover ends up as. 



#######################################################################################

# Count number of disturbances
trajectory_length <- length(sequences_raw)

# Create a list to store results
results_list <- vector("list", length = repeats)

# Get SEM to vary starting coral cover according to survey method of last known data point. 
if (Method == 'Manta' |
    Method == 'RJFMP' |
    Method == 'COTS') {
  a <- model.manta$coefficients[3]
  b <- model.manta$coefficients[2]
  c <- model.manta$coefficients[1]
} else if (Method == 'AIMS Transect') {
  a <- model.ltmp$coefficients[3]
  b <- model.ltmp$coefficients[2]
  c <- model.ltmp$coefficients[1]
} else if (Method == 'MMP') {
  a <- model.mmp$coefficients[3]
  b <- model.mmp$coefficients[2]
  c <- model.mmp$coefficients[1]
} else if (Method == 'RHIS' | 
           Method == 'Census') {
  a <- model.grc$coefficients[3]
  b <- model.grc$coefficients[2]
  c <- model.grc$coefficients[1]
} else if (Method == 'Reefmod') {
  a <- model.reefmod$coefficients[3]
  b <- model.reefmod$coefficients[2]
  c <- model.reefmod$coefficients[1]
} 

SEM <- max(0.1, a * starting_cover^2 + b * starting_cover + c)

# Loop through the number of repeats
for (i in 1:repeats) {
  # Modify starting_cover based on precision of survey method, bounded by 0.01 and 99.9
  starting_cover_i <- min(99.9,
                          max(0.01,
                              runif(1, min = starting_cover - SEM,
                                    max = starting_cover + SEM)))
  
  
  sequences_temp <- sequences_raw
  
  new_cover_values <- numeric(trajectory_length + 1)  # +1 to include step 0
  
  ## Modify the first sequence to include starting coral cover
  # Group starting_cover into bins of size 5
  bin_cover_temp <- floor(starting_cover_i / 5) * 5
  
  # Modify the first sequence to include the binned cover value
  sequences_temp[1] <- paste0(sequences_temp[1], "_Cover_", bin_cover_temp)
  
  new_cover_temp <- starting_cover_i
  
  # Step 0
  new_cover_values[1] <- starting_cover_i
  
  # Loop through each trajectory
  for (step in 1:trajectory_length) {
    
    key <- sequences_temp[step] # Check which lookup table the sequence appears in
    
    if (key %in% names(lookup_cyclone)) {
      matched_df <- lookup_cyclone[[key]]
    } else if (key %in% names(lookup_dhw)) {
      matched_df <- lookup_dhw[[key]]
    } else {
      stop(paste("Key not found in either lookup:", key))
    }
    
    if (nrow(matched_df) == 0 || length(matched_df$Cover_diff) == 0) { # error check
      stop(paste("Empty matched_df for key:", key))
    }
    
    sampled_diff <- sample(matched_df$Cover_diff, 1)
    new_cover_temp <- min(99.9, max(0.01, new_cover_temp + sampled_diff)) # limit coral cover to be between 0.01 and 99.9.
    
    new_cover_values[step + 1] <- new_cover_temp  # step + 1 because index 1 is step 0
    
    if (step < trajectory_length) {
      bin_cover_temp <- floor(new_cover_temp / 5) * 5
      sequences_temp[step + 1] <- paste0(sequences_temp[step + 1], "_Cover_", bin_cover_temp)
    }
  }
  
  results_list[[i]] <- tibble(
    run = i,
    step = 0:trajectory_length,
    cover = new_cover_values
  )
}

# Combine into one long data frame for plotting
results_df <- bind_rows(results_list)

# Add year to the results
## Get years with disturbance (Cyclone > 0 or DHW_cat > 0), and include Year 0 for step 0
event_years <- regime$Year[regime$Cyclone > 0 | as.numeric(as.character(regime$DHW_cat)) > 0]

## Create step-to-year map so plots correctly map year end with coral cover. 
if (min(event_years) == min(regime$Year) && length(event_years) > 1) { # for when disturbance is on first year followed by more
  step_to_year <- tibble(
    step = c(0, 1:length(event_years)),
    year = c(min(regime$Year) - 1, (event_years - 1)[-1], max(regime$Year))
  )
} else  if (min(event_years) == min(regime$Year)) { # for when disturbance is on first year only
  step_to_year <- tibble(
    step = c(0, 1:length(event_years)),
    year = c(min(regime$Year) - 1, max(regime$Year))
  )
} else { # for all other scenarios
  step_to_year <- tibble(
    step = c(0, 1:(length(event_years) + 1)),
    year = c(min(regime$Year) - 1, event_years - 1, max(regime$Year))
  )
}


# Join to results_df
results_df <- results_df %>%
  left_join(step_to_year, by = "step")

############
# Summary of uncertainty
# Find the final year (should be the max in the step-to-year mapping)
final_year <- max(results_df$year, na.rm = TRUE)

# Calculate mean and SD
final_summary <- results_df %>%
  filter(year == final_year) %>%
  summarise(
    mean_cover = mean(cover, na.rm = TRUE),
    sd_cover = sd(cover, na.rm = TRUE),
    lowest_cover = min(cover, na.rm = TRUE),
    highest_cover = max(cover, na.rm = TRUE),
    lower_estimate = as.numeric(quantile(cover, 0.16, na.rm = TRUE)),
    upper_estimate = as.numeric(quantile(cover, 0.84, na.rm = TRUE)),
    n = n()
  )

return(list(results_df = results_df, summary = final_summary))

}
