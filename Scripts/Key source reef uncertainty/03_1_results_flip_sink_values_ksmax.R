
## Load data
uncertain_reefs <- read.csv(file.path("Outputs", "Key source reef uncertainty", "o2_KSR_risky_and_opportunity_reefs.csv"), stringsAsFactors = FALSE)
load('Data/keysourcereefs2023.RData')
library(Matrix)
coral_covers <- read.csv(file.path("Outputs", "Use lookup tables" , "3.2_lower_and_upper_coral_cover_estimates.csv"), stringsAsFactors = FALSE)
cover_min <- coral_covers$CoverLower_estimate
cover_max <- coral_covers$CoverUpper_estimate
cover_mean <- coral_covers$Mean


## KSR parameters (copied from o1)
con <- MEAN_CONNECT
# con <- con - diag(diag(con))  # optional self-seeding removal
reefsizes <- reefsizes
alpha <- 1
critic <- 20
region <- region

cat('Starting computation (KSmax)...\n')
start_time <- Sys.time()

row_col <- dim(con)
row <- row_col[1]
col <- row_col[2]
out <- matrix(0, nrow = row, ncol = 11)

# Placeholder for neediness initialization (computed inside function per flip)
need <- rep(0, row)

# Prepare list to collect rows per uncertain reef
out_local <- vector("list", nrow(uncertain_reefs))

# Source helper (KSmax version)
source('Scripts/Key source reef uncertainty/Functions/run_KS_algorithm_flip_max.R')

# Parallel setup
library(foreach)
library(doSNOW)
cl <- makeCluster(parallel::detectCores() - 1)
registerDoSNOW(cl)

clusterExport(cl, c('con','reefsizes','critic','col','cover_mean','cover_min','cover_max',
                    'uncertain_reefs','run_KS_algorithm_flip_max'), envir = environment())

# Progress bar
pb <- txtProgressBar(max = nrow(uncertain_reefs), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

out_local <- foreach(i = 1:nrow(uncertain_reefs), .options.snow = opts, .packages = c('Matrix')) %dopar% {
  # Reset cover for this source reef
  cover <- cover_mean

  # Identify source reef ID (GBR reef number)
  source_ID <- uncertain_reefs[i, "Reef"]

  # Identify sinks (non-zero connectivity) excluding self
  sink_columns <- which(con[source_ID, ] != 0)
  sink_columns <- sink_columns[sink_columns != source_ID]

  # Preallocate: 1 for source flip + N for sinks
  num_sinks <- uncertain_reefs[i, "Number_of_sinks"]
  tmp_out <- matrix(0, nrow = num_sinks + 1, ncol = 12)

  ## Set starting parameters for KSmax baseline
  for (j in 1:col) {
    if (source_ID == j) {
      cover[j] <- cover_max[j]   # source at max
    } else if (con[source_ID, j] > 0) {
      cover[j] <- cover_min[j]   # sinks at min
    } else {
      cover[j] <- cover_mean[j]  # others at mean
    }
  }

  ## Flip source coral cover (to min) and evaluate
  cover_flip <- cover
  cover_flip[source_ID] <- cover_min[source_ID]
  tmp_out[1, ] <- run_KS_algorithm_flip_max(cover_flip, source_ID, NA_real_, critic, con, reefsizes, col, uncertain_reefs)

  ## Iteratively flip each sink (to max) and evaluate
  for (m in 1:uncertain_reefs[i, "Number_of_sinks"]) {
    cover_flip <- cover
    sink_ID <- sink_columns[m]
    cover_flip[sink_ID] <- cover_max[sink_ID]
    tmp_out[m + 1, ] <- run_KS_algorithm_flip_max(cover_flip, source_ID, sink_ID, critic, con, reefsizes, col, uncertain_reefs)
  }

  tmp_out
}

close(pb)
stopCluster(cl)

end_time <- Sys.time()
cat('Elapsed time (KSmax):', as.numeric(end_time - start_time, units = "secs"), 'seconds\n')

# Concatenate results
out <- do.call(rbind, out_local)

# Label columns
out_table_ksmax <- as.data.frame(out)
colnames(out_table_ksmax) <- c('SourceReefID', 'Flipped_SinkReefID',
                               'EnrichUnweighted', 'Flipped_Score',
                               'SourceNeed', 'SinkNeedAfterFlip', 'NumSinks', 'NumNeedySinks',
                               'Source_KSmax', 'Diff_FlippedScore_KSmax', 'Source_KSdiff',
                               'Pc_KSdiff_accounted_by_flip')

# Add Opportunity / Risk
out_table_ksmax$Source_opportunity <- NA_real_
out_table_ksmax$Source_risk        <- NA_real_
match_idx <- match(out_table_ksmax$SourceReefID, uncertain_reefs$Reef)
valid_matches <- !is.na(match_idx)
out_table_ksmax$Source_opportunity[valid_matches] <- uncertain_reefs$Opportunity[match_idx[valid_matches]]
out_table_ksmax$Source_risk[valid_matches]        <- uncertain_reefs$Risk[match_idx[valid_matches]]

# Write output
write.csv(out_table_ksmax, file.path("Outputs", "Key source reef uncertainty", "o3_1_results_flip_KSmax.csv"), row.names = FALSE)
