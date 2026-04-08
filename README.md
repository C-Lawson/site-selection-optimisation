# Site Selection Algorithms

Code for assessing uncertainty in Key Source Reef (KSR) rankings on the Great Barrier Reef. Accompanies the associated publication.

## Dependencies

R packages: `tidyverse`, `Matrix`, `foreach`, `doSNOW`, `ggplot2`, `ompr`, `ompr.roi`, `ROI.plugin.highs`

## Run Order

Scripts should be run sequentially. Outputs from each step are used as inputs to the next.

### 1. Estimate coral cover uncertainty
`Scripts/Use lookup tables/3.2_Results_Use_lookup_tables.R`

Uses disturbance-trajectory lookup tables to propagate coral cover uncertainty forward from each reef's last survey date.

### 2. Compute KSR score ranges
`Scripts/Key source reef uncertainty/01_recoveryreef2_test_for_min_and_max.R`

Runs the KSR algorithm under min, mean, and max coral cover scenarios to bound each reef's importance score.

### 3. Identify uncertain reefs
`Scripts/Key source reef uncertainty/02_results_get_uncertain_source_reefs.R`

Ranks reefs, computes net rank movement under uncertainty, and identifies reefs at risk of moving in or out of the KSR zone.

### 4. Attribute uncertainty to individual reefs
`Scripts/Key source reef uncertainty/03_results_flip_sink_values_ksmin.R`
`Scripts/Key source reef uncertainty/03_1_results_flip_sink_values_ksmax.R`

Iteratively flips coral cover of each reef in a source network to quantify its contribution to score uncertainty (from KSmin and KSmax baselines).

### 5. Summarise attribution results
`Scripts/Key source reef uncertainty/04_results_identify_reefs_of_interest.R`

Merges KSmin and KSmax flip results and computes mean uncertainty attribution per reef.

### 6. Optimise survey site selection
`Scripts/Key source reef uncertainty/05.2_optimiser_function.R`

Integer linear program to select K reefs whose surveying would resolve uncertainty in the greatest number of KSR networks.

## Data

Input data files are in `Data/`. Intermediate and final outputs are written to `Outputs/`.
