# Aim: Get mean percentage and absolute KS diff accounted for by each flipped reef
# in the source network (average of KSmin and KSmax flipping scenarios).

# Input CSVs 
ksmin_path <- file.path("Outputs", "Key source reef uncertainty", "o3_results_flip_KSmin.csv")
ksmax_path <- file.path("Outputs", "Key source reef uncertainty", "o3_1_results_flip_KSmax.csv")

ksmin <- read.csv(ksmin_path, stringsAsFactors = FALSE, check.names = FALSE)
ksmax <- read.csv(ksmax_path, stringsAsFactors = FALSE, check.names = FALSE)

# Add absolute KS diff accounted for by flip (out of 2), not percent
ksmin$Abs_KSdiff_accounted_by_flip <- ksmin$Source_KSdiff * (ksmin$Pc_KSdiff_accounted_by_flip / 100)
ksmax$Abs_KSdiff_accounted_by_flip <- ksmax$Source_KSdiff * (ksmax$Pc_KSdiff_accounted_by_flip / 100)

# Flipped_Reef: use Flipped_SinkReefID, but where NA use SourceReefID
is_na_ksmin <- is.na(ksmin$Flipped_SinkReefID)
is_na_ksmax <- is.na(ksmax$Flipped_SinkReefID)

ksmin$Flipped_Reef <- ksmin$Flipped_SinkReefID
ksmax$Flipped_Reef <- ksmax$Flipped_SinkReefID
ksmin$Flipped_Reef[is_na_ksmin] <- ksmin$SourceReefID[is_na_ksmin]
ksmax$Flipped_Reef[is_na_ksmax] <- ksmax$SourceReefID[is_na_ksmax]

# Type: 'Sink' by default, 'Source' where flipped sink reef id is NaN
ksmin$Type <- "Sink"
ksmax$Type <- "Sink"
ksmin$Type[is_na_ksmin] <- "Source"
ksmax$Type[is_na_ksmax] <- "Source"

# Select columns before outer join
left_df <- ksmin[, c("SourceReefID", "Flipped_Reef",
                     "Abs_KSdiff_accounted_by_flip", "Pc_KSdiff_accounted_by_flip",
                     "Type", "Source_KSdiff", "Source_opportunity", "Source_risk")]
right_df <- ksmax[, c("SourceReefID", "Flipped_Reef",
                      "Abs_KSdiff_accounted_by_flip", "Pc_KSdiff_accounted_by_flip",
                      "Type")]

merged <- merge(left_df, right_df,
                by = c("SourceReefID", "Flipped_Reef"),
                all = TRUE, suffixes = c("_left", "_right"), sort = TRUE)

# Rename columns
colnames(merged)[colnames(merged) == "Abs_KSdiff_accounted_by_flip_left"] <- "Abs_diff_ksmin"
colnames(merged)[colnames(merged) == "Abs_KSdiff_accounted_by_flip_right"] <- "Abs_diff_ksmax"
colnames(merged)[colnames(merged) == "Pc_KSdiff_accounted_by_flip_left"] <- "Pc_diff_ksmin"
colnames(merged)[colnames(merged) == "Pc_KSdiff_accounted_by_flip_right"] <- "Pc_diff_ksmax"
colnames(merged)[colnames(merged) == "Type_left"] <- "Type_ksmin"
colnames(merged)[colnames(merged) == "Type_right"] <- "Type_ksmax"

# Ensure numeric NA -> NaN in the 4 diff columns before mean
num_cols <- c("Abs_diff_ksmin","Abs_diff_ksmax","Pc_diff_ksmin","Pc_diff_ksmax")
for (cc in num_cols) {
  vals <- merged[[cc]]
  is_na <- is.na(vals)
  if (any(is_na)) merged[[cc]][is_na] <- NaN
}

# Compute means (row-wise); with NaN present, results remain NaN
merged$Abs_diff_mean <- (merged$Abs_diff_ksmin + merged$Abs_diff_ksmax) / 2
merged$Pc_diff_mean  <- (merged$Pc_diff_ksmin + merged$Pc_diff_ksmax) / 2

# Order columns
final_cols <- c(
  "SourceReefID","Flipped_Reef",
  "Abs_diff_ksmin","Pc_diff_ksmin","Type_ksmin",
  "Source_KSdiff","Source_opportunity","Source_risk",
  "Abs_diff_ksmax","Pc_diff_ksmax","Type_ksmax",
  "Abs_diff_mean","Pc_diff_mean"
)

# Some rows may be present only in ksmax -> Source_KSdiff/opportunity/risk become NA (to be written as NaN)
# Convert those NA to NaN
for (cc in c("Source_KSdiff","Source_opportunity","Source_risk")) {
  if (cc %in% names(merged)) {
    vals <- merged[[cc]]
    is_na <- is.na(vals)
    if (any(is_na)) merged[[cc]][is_na] <- NaN
  }
}

out <- merged[, final_cols]

# Sort tuple-wise by SourceReefID then Flipped_Reef
out <- out[order(out$SourceReefID, out$Flipped_Reef), ]

# Stats
n_distinct(out$SourceReefID) 
n_distinct(out$Flipped_Reef) 


# Write CSV with no quotes and NaN for missing numeric values; no row names
out_path <- file.path("Outputs", "Key source reef uncertainty", "o4_reefs_of_interest_summaries.csv")
write.csv(out, out_path, row.names = FALSE, quote = FALSE, na = "NaN")
