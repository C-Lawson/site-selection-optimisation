# Aim: run the key source reef algorithm specifically for the purpose of
# iteratively flipping coral covers between their lower and upper estimates
# to determine the impact that any one source or sink reef has the
# importance score of the source reef (starting from KSmax scenario).

run_KS_algorithm_flip_max <- function(cover_flip, source_ID, sink_ID, critic, con, reefsizes, col, uncertain_reefs) {
  # Calculate neediness
  need <- rep(0, col) # reset neediness scores
  below_critic <- cover_flip < critic
  need[below_critic] <- (critic - cover_flip[below_critic]) / critic

  # Larval input and flux (sparse-aware)
  con <- as(con, "sparseMatrix")
  larvae_export <- as.vector(cover_flip * reefsizes)
  # Scale rows of con by larvae_export (cover*size per source)
  larvalinputs <- Matrix::Diagonal(x = larvae_export) %*% con
  larvalinput_reef <- Matrix::colSums(larvalinputs)
  larvalflux_reef <- larvalinput_reef / reefsizes

  # Enrichment calculations
  enrichunweight <- 0
  enrichweight <- 0
  sink_indices <- which(con[source_ID, ] > 0)
  for (j in sink_indices) {
    inputflux <- larvalinputs[source_ID, j] / reefsizes[j]
    propinputflux <- inputflux / larvalflux_reef[j]
    enrichunweight <- enrichunweight + propinputflux
    enrichweight <- enrichweight + (propinputflux * need[j])
  }

  # Count connected reefs
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
  Source_KSmax <- uncertain_reefs$Max[uncertain_reefs$Reef == source_ID]
  Diff_FlippedScore_KSmax <- Source_KSmax - enrichweight
  Source_KSdiff <- uncertain_reefs$Difference[uncertain_reefs$Reef == source_ID]
  Percent_KSdiff_accounted <- (Diff_FlippedScore_KSmax / Source_KSdiff) * 100

  # define sink_need for output 
  if (is.na(sink_ID)) {
    sink_need <- NaN
  } else {
    sink_need <- need[sink_ID]
  }

  # Output row
  out_row <- c(source_ID, sink_ID, enrichunweight, enrichweight,
               need[source_ID], sink_need, numreefs, numreefsneeded,
               Source_KSmax, Diff_FlippedScore_KSmax, Source_KSdiff, Percent_KSdiff_accounted)
  return(out_row)
}
