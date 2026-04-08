# Aim: After Identifying the KSR networks that are most at risk and at opportunity of moving in or out of the KSR zone, and iteratively flipping each reef j in that KSR network i to determine the percent of KSR i uncertainty given by the possible coral cover range in reef j, we now want to select the 100 reefs (or whatever number) which if their coral cover is known will clarify the state of the most number of KSR networks. This is different to simply selecting the most influential reef j in the most at risk or most opportunity KSR i (this was done in script 05.1).
library(tidyverse)
library(Matrix)
library(ompr)
library(ompr.roi)
library(ROI.plugin.highs)

# Function to run the optimization
run_ksr_optimization <- function(T_thr, K, MIN_CORAL_RANGE = 5, verbose = FALSE) {
  # Purpose: Run KSR network optimization to select K reefs that maximize total value of resolved sources
  # Inputs:
  #   T_thr: Resolution threshold (0-100)
  #   K: Number of reefs to survey
  #   MIN_CORAL_RANGE: Minimum coral cover range to consider a reef valid for selection
  #   verbose: Whether to print solver output
  # Returns: List with elements:
  #   - selected_reefs: data frame of selected reef IDs
  #   - resolved_sources: data frame of resolved source reef IDs
  #   - n_resolved: number of resolved sources
  #   - total_value: total value of resolved sources
  

  # Load data (these must exist in the working directory outputs)
  dat <- read.csv("Outputs/Key source reef uncertainty/o4_reefs_of_interest_summaries.csv")
  dat.coral <- read.csv('Outputs/Use lookup tables/3.2_lower_and_upper_coral_cover_estimates.csv')

  ## Add a column to dat called 'coral range' that takes value from dat.coral$Range, matching dat$Flipped_Reef to dat.coral$Reef.Number
  dat <- dat %>%
    left_join(dat.coral %>% select(Reef.Number, Range),
              by = c("Flipped_Reef" = "Reef.Number")) %>%
    rename(flipped_coral_range = Range)

  ## Remove reefs with coral range below the minimum threshold
  dat <- dat %>%
    filter(flipped_coral_range >= MIN_CORAL_RANGE)

  # Prepare value column
  ## Reverse the sign of dat$Source_risk
  dat$Source_risk <- -dat$Source_risk

  ## get M, the minimum value from Source_risk and Source_opportunity.
  M <- min(c(min(dat$Source_risk, na.rm = TRUE),
             min(dat$Source_opportunity, na.rm = TRUE)))

  dat$risk_value <- dat$Source_risk - M + 1
  dat$opp_value <- dat$Source_opportunity - M + 1
  dat <- dat %>%
    mutate(value = coalesce(risk_value, opp_value))

  # 1) Build (sparse) contribution matrix C: rows = sources i, cols = reefs j
  triples <- dat %>%
    transmute(SourceReefID,
              Flipped_Reef,
              c = ifelse(pmax(Pc_diff_mean, 0) < 1e-6, 0, pmax(Pc_diff_mean, 0))) %>%
    filter(!is.na(SourceReefID), !is.na(Flipped_Reef), !is.na(c), c > 0) %>%
    group_by(SourceReefID, Flipped_Reef) %>%
    summarise(c = sum(c), .groups = "drop")

  source_ids <- sort(unique(triples$SourceReefID))
  reef_ids   <- sort(unique(triples$Flipped_Reef))

  triples <- triples %>%
    mutate(i = match(SourceReefID, source_ids),
           j = match(Flipped_Reef, reef_ids))

  C <- sparseMatrix(i = triples$i,
                    j = triples$j,
                    x = triples$c,
                    dims = c(length(source_ids), length(reef_ids)))

  n_sources <- nrow(C)
  n_reefs   <- ncol(C)

  # Get value for each source reef
  source_values <- dat %>%
    filter(SourceReefID %in% source_ids) %>%
    group_by(SourceReefID) %>%
    summarise(value = first(value), .groups = "drop") %>%
    arrange(match(SourceReefID, source_ids)) %>%
    pull(value)

  # 2) ILP model: maximise value of resolved sources with K reefs selected
  model <- MIPModel() %>%
    add_variable(x[j], j = 1:n_reefs, type = "binary") %>%
    add_variable(y[i], i = 1:n_sources, type = "binary") %>%

    # Constraints
    add_constraint(sum_expr(x[j], j = 1:n_reefs) == K) %>%
    add_constraint(sum_expr(C[i, j] * x[j], j = 1:n_reefs) >= T_thr * y[i],
                   i = 1:n_sources) %>%

    # Objective: maximize total value
    set_objective(sum_expr(source_values[i] * y[i], i = 1:n_sources), "max")

  result <- solve_model(model, with_ROI(solver = "highs", verbose = verbose))

  # 3) Read solution
  selected_reefs <- get_solution(result, x[j]) %>%
    filter(value > 0.5) %>%
    transmute(reef_col = j,
              Flipped_Reef = reef_ids[j])

  resolved_sources <- get_solution(result, y[i]) %>%
    filter(value > 0.5) %>%
    transmute(src_row = i,
              SourceReefID = source_ids[i])

  # Extract total value
    total_val <- ompr::objective_value(result)

  # Return results as list
  return(list(
    selected_reefs = selected_reefs,
    resolved_sources = resolved_sources,
    n_resolved = nrow(resolved_sources),
    total_value = total_val
  ))
}


# ============================================================================
# STANDALONE USAGE (for interactive testing)
# ============================================================================
# Uncomment below to run the optimization with fixed parameters
# 
# # Parameters
# T_thr <- 80   # resolution threshold
# K <- 100      # number of reefs to survey
# MIN_CORAL_RANGE <- 5  # Minimum coral cover range to consider a reef valid for selection
# 
# result <- run_ksr_optimization(T_thr = T_thr, K = K, MIN_CORAL_RANGE = MIN_CORAL_RANGE, verbose = TRUE)
# 
# cat("Selected reefs:", nrow(result$selected_reefs), "\n",
#     "Resolved sources:", result$n_resolved, "\n",
#     "Total value:", result$total_value, "\n")
# 
# # Write outputs
# filename <- paste0("Outputs/Key source reef uncertainty/o5.2_ILP_optimisation_selected_reefs",
#                    "_K", K, ".csv")
# write.csv(result$selected_reefs, filename, row.names = FALSE)
# 
# filename_resolved <- paste0("Outputs/Key source reef uncertainty/o5.2_ILP_optimisation_resolved_reefs",
#                             "_K", K, ".csv")
# write.csv(result$resolved_sources, filename_resolved, row.names = FALSE)
