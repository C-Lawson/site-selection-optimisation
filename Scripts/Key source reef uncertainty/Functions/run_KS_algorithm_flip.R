# Aim: run the key source reef algorithm specifically for the purpose of
# iteratively flipping coral covers between their lower and upper estimates
# to determine the impact that any one source or sink reef has the
# importance score of the source reef.

run_KS_algorithm_flip <- function(cover_flip, source_ID, sink_ID, critic, con, reefsizes, col, uncertain_reefs) {
  
  ## Calculate neediness
  need <- rep(0, col) # reset neediness scores
  # Identify needy reefs: apply logical index of where cover is less than critic
  below_critic <- cover_flip < critic
  # Now apply formula to calculate need scores for needy reefs, only at needy
  # reefs (non-needy reefs remain on 0 neediness).
  need[below_critic] <- (critic - cover_flip[below_critic]) / critic
  
  ## Larval input and flux
  # create weighted contributions to connectivity
  # Ensure con is sparse
  con <- as(con, "sparseMatrix")  # if it isn't already
  
  # Compute - scale each ROW of con by (cover_flip * reefsizes)
  larvae_export <- as.vector(cover_flip * reefsizes)
  larvalinputs <- Matrix::Diagonal(x = larvae_export) %*% con
  larvalinput_reef <- Matrix::colSums(larvalinputs)
  larvalflux_reef <- larvalinput_reef / reefsizes # convert to a flux based on size
  
  ## Enrichment calculations
  enrichunweight <- 0 # calc enrichment done by each source
  enrichweight <- 0
  
  # Get non-zero connections for this source reef (sparse-aware)
  sink_indices <- which(con[source_ID, ] > 0)
  
  for (j in sink_indices) {
    # add an if to ensure connectivity?
    inputflux <- larvalinputs[source_ID, j] / reefsizes[j] # convert to a flux
    propinputflux <- inputflux / larvalflux_reef[j] # express input as proprotion of total flux to that reef
    enrichunweight <- enrichunweight + propinputflux #create total
    enrichweight <- enrichweight + (propinputflux * need[j]) # weight by need at the sink
  }
  
  
  ## Count connected reefs
  # identify number of reefs each source contacts and then the weighted
  # number which ignores those that don't need it
  numreefs <- 0
  numreefsneeded <- 0
  for (j in sink_indices) {
    if (j != source_ID) {
      numreefs <- numreefs + 1
      if (cover_flip[j] <= critic) {
        numreefsneeded <- numreefsneeded + 1
      }
    }
  }
  
  
  # Metrics from uncertain_reefs
  Source_KSmin <- uncertain_reefs$Min[uncertain_reefs$Reef == source_ID]
  Diff_FlippedScore_KSmin <- enrichweight - Source_KSmin
  Source_KSdiff <- uncertain_reefs$Difference[uncertain_reefs$Reef == source_ID]
  Percent_KSdiff_accounted <- (Diff_FlippedScore_KSmin / Source_KSdiff) * 100
  
  # define sink_need for output, even if sink_need is NaN
  if (is.na(sink_ID)) {
    sink_need <- NaN  
  } else {
    sink_need <- need[sink_ID]
  }
  
  # Output row
  out_row <- c(source_ID, sink_ID, enrichunweight, enrichweight,
               need[source_ID], sink_need, numreefs, numreefsneeded,
               Source_KSmin, Diff_FlippedScore_KSmin, Source_KSdiff, Percent_KSdiff_accounted)
  
  return(out_row)
}
