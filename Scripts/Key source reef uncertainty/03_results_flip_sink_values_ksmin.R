# Aim: Step 4 in our workflow. For each reef identified in previous script
# as either top risk or opportunity score, we iteratively flip its sink reefs' coral cover
# from min to max (and vice verse) to assess the portion of difference in
# importance score of the source reef is attributable to uncertainty in
# each of its sinks. This needs to be done for both KSmin and KSmax starting scenarios - this script only does KSmin.  


## Data
uncertain_reefs <- read.csv(file.path("Outputs","Key source reef uncertainty", "o2_KSR_risky_and_opportunity_reefs.csv"), stringsAsFactors = FALSE)
if (!exists("MEAN_CONNECT") ||
    !exists("reefsizes") ||
    !exists("region") ||
    !exists("ksr_alpha") ||
    !exists("ksr_critic") ||
    !exists("ksr_parallel_cores")) {
  stop("Run Scripts/00_MASTER.R instead of this script.")
}
coral_covers <- read.csv(file.path("Outputs", "Use lookup tables" , "3.2_lower_and_upper_coral_cover_estimates.csv"), stringsAsFactors = FALSE)
cover_min <- coral_covers$CoverLower_estimate # if cover is zero then algorithm says source reefs aren't connected to any sinks so make minimum possible value slightly higher than zero.
cover_max <- coral_covers$CoverUpper_estimate # if cover is zero then algorithm says source reefs aren't connected to any sinks so make minimum possible value slightly higher than zero.
cover_mean <- coral_covers$Mean


con <- MEAN_CONNECT
con <- con - Matrix::Diagonal(x = diag(con))  # Remove self-seeding.
reefsizes <- reefsizes
alpha <- ksr_alpha
critic <- ksr_critic
region <- region


## 
cat('Starting computation...\n')
start_time <- Sys.time()
    
row_col <- dim(con)
row <- row_col[1]
col <- row_col[2]
out <- matrix(0, nrow = row, ncol = 11)

# calculate weights of each reef based on current cover
need <- rep(0, row)

# Each cell holds all rows for one uncertain reef
out_local <- vector("list", nrow(uncertain_reefs))

# Source helper function
source('Scripts/Key source reef uncertainty/Functions/run_KS_algorithm_flip.R')

# Set up parallel processing
cl <- makeCluster(ksr_parallel_cores)
registerDoSNOW(cl)

# Export necessary variables and function to workers
clusterExport(cl, c('con', 'reefsizes', 'critic', 'col', 'cover_mean', 'cover_min', 'cover_max', 
                    'uncertain_reefs', 'run_KS_algorithm_flip'), 
              envir = environment())

# Progress bar setup
pb <- txtProgressBar(max = nrow(uncertain_reefs), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

out_local <- foreach(i = 1:nrow(uncertain_reefs), 
                     .options.snow = opts,
                     .packages = c('Matrix')) %dopar% {
  # Calculate score for each reef.
  cat(sprintf('Uncertain source reef: %d\n', i)) # visualise current iteration
  cover <- cover_mean    # <- reset every time we start reef i
  

  # Get reef number in matrix for source reef i (because
  # uncertain_reefs is filtered, i is the row number of this filtered
  # list so it does not correspond to reef number out of 3806.  
  source_ID <- uncertain_reefs[i, "Reef"]

  # Find sinks (non-zero columns in row i of 'con')
  sink_columns <- which(con[source_ID, ] != 0)  # these are the reef numbers that are sinks for source reef i (source_ID)
  sink_columns <- sink_columns[sink_columns != source_ID]  # remove occurrences equal to i (i.e. don't list the source as a sink to itself even if there is self-seeding, we manually flip the value of the source itself later).
  
  # Preallocate space: 1 row for source flip + N rows for each sink
  num_sinks <- uncertain_reefs[i, "Number_of_sinks"]
  tmp_out <- matrix(0, nrow = num_sinks + 1, ncol = 12)


  ## Set starting parameters before flipping
  for (j in 1:col) { # Modify coral cover for each sink of reef i so the score can be calculated.
    if (source_ID == j) {
      cover[j] <- cover_min[j] # if 'sink' reef is the source (reef i) then it gets minimum cover value (to minimise importance score).
    } else if (con[source_ID, j] > 0) { # if reef j is a sink to reef i
      cover[j] <- cover_max[j] # then the sink gets maximum cover value to minimise score of reef i.
    } else {
      cover[j] <- cover_mean[j]  # if reef j is not a sink of reef i then it just gets its normal value.
    }
  }


  ## Flip source coral cover on its own
  cover_flip <- cover
  cover_flip[source_ID] <- cover_max[source_ID]
  tmp_out[1, ] <- run_KS_algorithm_flip(cover_flip, source_ID, NA_real_, critic, con, reefsizes, col, uncertain_reefs)


  ## Iteratively flip each sink and run algorithm
  for (m in 1:uncertain_reefs[i,"Number_of_sinks"]) {
    cover_flip <- cover # Assign starting parameters for reef i
    
    # Find sink number m
    sink_ID <- sink_columns[m]

    # Now flip sink number m coral cover value
    cover_flip[sink_ID] <- cover_min[sink_ID]

    ## Coral covers are all good now for this 'flip iteration', now begin KS algorithm proper.
    tmp_out[m + 1, ] <- run_KS_algorithm_flip(cover_flip, source_ID, sink_ID, critic, con, reefsizes, col, uncertain_reefs)

  }

  # Store 
  out_local_i <- tmp_out
  
  # Return tmp_out for this source reef
  out_local_i

## Store output as minimum or maximum importance scores for each reef.
}

close(pb)
stopCluster(cl)

end_time <- Sys.time()
cat('Elapsed time:', as.numeric(end_time - start_time, units = "secs"), 'seconds\n')

## concatenate all outputs
out <- do.call(rbind, out_local)

## Label columns more clearly
out_table_ksmin <- as.data.frame(out)
colnames(out_table_ksmin) <- c('SourceReefID', 'Flipped_SinkReefID',
                               'EnrichUnweighted', 'Flipped_Score',
                               'SourceNeed', 'SinkNeedAfterFlip', 'NumSinks', 'NumNeedySinks',
                               'Source_KSmin', 'Diff_FlippedScore_KSmin', 'Source_KSdiff',
                               'Pc_KSdiff_accounted_by_flip')


# ---------- Add Source Reef Opportunity / Risk Scores ----------
# I need to add two new columns to out_table_ksmin: Source_opportunity and
# Source_risk. The values for these columns should be from
# uncertain_reefs.Opportunity and uncertain_reefs.Risk, where rows are
# matched according to uncertain_reefs.Reef and
# out_table_ksmin.SourceReefID. For each appearance in uncertain_reefs,
# there will be many appearances in out_table_ksmin

# Preallocate new columns
out_table_ksmin$Source_opportunity <- NA_real_
out_table_ksmin$Source_risk        <- NA_real_

# Match SourceReefID to uncertain_reefs.Reef
match_idx <- match(out_table_ksmin$SourceReefID, uncertain_reefs$Reef)

# Fill the new columns where matches exist
valid_matches <- !is.na(match_idx)
out_table_ksmin$Source_opportunity[valid_matches] <- uncertain_reefs$Opportunity[match_idx[valid_matches]]
out_table_ksmin$Source_risk[valid_matches]        <- uncertain_reefs$Risk[match_idx[valid_matches]]

## Write to file
write.csv(out_table_ksmin, file.path("Outputs", "Key source reef uncertainty", "o3_results_flip_KSmin.csv"), row.names = FALSE)
