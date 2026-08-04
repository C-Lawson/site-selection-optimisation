## Read in data
score_table <- read.csv(file.path('Outputs','Key source reef uncertainty','o1_KSR_scores_range_all_reefs.csv'), stringsAsFactors = FALSE)

# Set colours for plots
if (!exists("col_risk")) stop("Run Scripts/00_MASTER.R instead of this script.")
if (!exists("col_opp")) stop("Run Scripts/00_MASTER.R instead of this script.")
if (!exists("col_default")) stop("Run Scripts/00_MASTER.R instead of this script.")

# # Set parameters
if (!exists("KSR_zone_threshold")) stop("Run Scripts/00_MASTER.R instead of this script.")
if (!exists("opp_number")) stop("Run Scripts/00_MASTER.R instead of this script.")
if (!exists("risk_number")) stop("Run Scripts/00_MASTER.R instead of this script.")

##  Get KSR Rankings for predicted (mean) and max and min.
## I.e. where a reef's position in list is based on its mean vs means of all
## reefs, its min vs the min of all reefs, and its max vs the max of all
## reefs. 
score_table$MeanRank <- rank(-score_table$Mean, ties.method = 'min')
score_table$MinRank  <- rank(-score_table$Min,  ties.method = 'min')
score_table$MaxRank  <- rank(-score_table$Max,  ties.method = 'min')

# check unique ranks (because tiedrank will give the same rank to ties)
print(length(unique(score_table$MeanRank)))

# Calculate 'net' movement to indicate likely direction of change given uncertainty.
## Compute both differences
score_table$diff_minrank <- score_table$MeanRank - score_table$MinRank
score_table$diff_maxrank <- score_table$MeanRank - score_table$MaxRank

## Preallocate Net
score_table$Net <- rep(0, nrow(score_table))

## Case 1: Opposite signs → use full sum
opposite_signs <- sign(score_table$diff_minrank) != sign(score_table$diff_maxrank)
score_table$Net[opposite_signs] <- score_table$diff_minrank[opposite_signs] + score_table$diff_maxrank[opposite_signs]

## Case 2: Same sign → take the one with greater absolute value (preserve sign)
same_signs <- !opposite_signs
abs_min <- abs(score_table$diff_minrank)
abs_max <- abs(score_table$diff_maxrank)

use_min <- abs_min > abs_max
use_max <- !use_min

## Apply for same_signs only
idx_min <- same_signs & use_min
idx_max <- same_signs & use_max

score_table$Net[idx_min] <- score_table$diff_minrank[idx_min]
score_table$Net[idx_max] <- score_table$diff_maxrank[idx_max]

## Plot Mean Rankings vs Net movement
# Calculate Xth percentile of reef rankings (simple)
n_total <- nrow(score_table)
KSRzone <- as.numeric(quantile(seq_len(n_total), probs = KSR_zone_threshold/100, type = 4)) + 0.5

xA <- KSRzone
xB <- KSRzone

# Save a main scatter plot
p_base <- ggplot(score_table, aes(x = MeanRank, y = Net)) +
  geom_point(aes(fill = "All reefs"), shape = 21, color = "black", stroke = 0.3, size = 2, alpha = 0.7) +
  geom_vline(aes(xintercept = KSRzone, color = "KSR Zone Threshold", linetype = "KSR Zone Threshold"), linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  scale_x_reverse() +
  scale_color_manual(values = c(
    "All reefs" = col_default,
    "KSR Zone Threshold" = "red",
    "Top Opportunity reefs" = col_opp,
    "Top Risk reefs" = col_risk
  )) +
  scale_fill_manual(values = c(
    "All reefs" = col_default,
    "Top Opportunity reefs" = col_opp,
    "Top Risk reefs" = col_risk
  )) +
  scale_linetype_manual(values = c("KSR Zone Threshold" = "dashed")) +
  labs(
    x = expression("Predicted KS Rank (" * bar(R) * ")"),
    y = "Net Movement with Uncertainty",
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  theme_bw()

ggsave(
  filename = file.path('Outputs','Key source reef uncertainty','plot_net_vs_predicted_KSR.png'),
  plot = p_base,
  width = 1600 / 150,
  height = 1200 / 150,
  dpi = 150
)

## ---------- GET HOW LIKELY EACH REEF IS TO CROSS THE KSR LINE ----------
idx_A <- score_table$MeanRank > xA
idx_B <- score_table$MeanRank < xB

score_table$Opportunity <- rep(NA_real_, nrow(score_table))  # Pre-fill with NA

score_table$Opportunity[idx_A] <- score_table$Net[idx_A] - (score_table$MeanRank[idx_A] - KSRzone)

## Find top Opportunity reefs - the number selected = opp_number
## Use stable ordering with Reef as secondary key and NA last
if (all(is.na(score_table$Opportunity))) {
  top_Opp <- integer(0)
} else {
  values <- score_table$Opportunity
  # Select top opp_number largest values (exclude NA)
  non_na_idx <- which(!is.na(values))
  if (length(non_na_idx) > opp_number) {
    # Order decreasing and take first opp_number
    ord <- non_na_idx[order(values[non_na_idx], decreasing = TRUE)]
    top_Opp <- ord[seq_len(opp_number)]
  } else {
    top_Opp <- non_na_idx
  }
}

## Save a plot with selected opportunity reefs overlayed (tagged)
p_opp <- ggplot(score_table, aes(x = MeanRank, y = Net)) +
  geom_point(aes(fill = "All reefs"), shape = 21, color = "black", stroke = 0.3, size = 2, alpha = 0.7) +
  geom_vline(aes(xintercept = KSRzone, color = "KSR Zone Threshold", linetype = "KSR Zone Threshold"), linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_point(
    data = score_table[top_Opp, , drop = FALSE],
    aes(fill = "Top Opportunity reefs"),
    shape = 21,
    color = "black",
    stroke = 0.3,
    size = 2,
    alpha = 0.7
  ) +
  scale_x_reverse() +
  scale_color_manual(values = c(
    "All reefs" = col_default,
    "KSR Zone Threshold" = "red",
    "Top Opportunity reefs" = col_opp,
    "Top Risk reefs" = col_risk
  )) +
  scale_fill_manual(values = c(
    "Top Opportunity reefs" = col_opp,
    "Top Risk reefs" = col_risk,
    "All reefs" = col_default
  )) +
  scale_linetype_manual(values = c("KSR Zone Threshold" = "dashed")) +
  labs(
    x = expression("Predicted KSR Rank (" * bar(R) * ")"),
    y = "Net Movement with Uncertainty",
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(
  filename = file.path('Outputs','Key source reef uncertainty','plot_net_vs_predicted_KSR_with_selection.png'),
  plot = p_opp,
  width = 1600 / 150,
  height = 1200 / 150,
  dpi = 150
)

## Same but for risk of dropping out of the KSR zone
score_table$Risk <- rep(NA_real_, nrow(score_table))  # Pre-fill with NA

score_table$Risk[idx_B] <- score_table$Net[idx_B] - (score_table$MeanRank[idx_B] - KSRzone)

## Find top Risk reefs - the number selected = risk_number.
## Find top Risk reefs - the number selected = risk_number
## Use stable ordering with Reef as secondary key and NA last
if (all(is.na(score_table$Risk))) {
  top_Risk <- integer(0)
} else {
  values_r <- score_table$Risk
  non_na_idx_r <- which(!is.na(values_r))
  # For risk we want the most negative (smallest) values
  if (length(non_na_idx_r) > risk_number) {
    ord_r <- non_na_idx_r[order(values_r[non_na_idx_r], decreasing = FALSE)]
    top_Risk <- ord_r[seq_len(risk_number)]
  } else {
    top_Risk <- non_na_idx_r
  }
}

## Plot them on top of your existing figure (separate file to avoid overwriting)
p_sel <- ggplot(score_table, aes(x = MeanRank, y = Net)) +
  geom_point(aes(fill = "All reefs"), shape = 21, color = "black", stroke = 0.3, size = 2, alpha = 0.7) +
  geom_vline(aes(xintercept = KSRzone, color = "KSR Zone Threshold", linetype = "KSR Zone Threshold"), linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_point(
    data = score_table[top_Opp, , drop = FALSE],
    aes(fill = "Top Opportunity reefs"),
    shape = 21,
    color = "black",
    stroke = 0.3,
    size = 2,
    alpha = 0.7
  ) +
  geom_point(
    data = score_table[top_Risk, , drop = FALSE],
    aes(fill = "Top Risk reefs"),
    shape = 21,
    color = "black",
    stroke = 0.3,
    size = 2,
    alpha = 0.7
  ) +
  annotate("rect", xmin = 0, xmax = 500, ymin = -1500, ymax = 500, fill = NA, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_reverse(limits = c(max(score_table$MeanRank) + 5, min(score_table$MeanRank) - 5), expand = c(0, 0)) +
  scale_color_manual(values = c(
    "KSR Zone Threshold" = "red"
  )) +
  scale_fill_manual(values = c(
    "Top Opportunity reefs" = col_opp,
    "Top Risk reefs" = col_risk,
    "All reefs" = col_default
  )) +
  scale_linetype_manual(values = c("KSR Zone Threshold" = "dashed")) +
  labs(
    x = expression("Predicted KSR Rank (" * bar(R) * ")"),
    y = "Net Movement with Uncertainty",
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  theme_bw() +
  theme(legend.position = "right") +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        color = "black",
        stroke = 0.3,
        alpha = 0.7
      )
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(
        linetype = "dashed",
        linewidth = 0.6
      )
    ),
    linetype = "none"
  )

ggsave(
  filename = file.path('Outputs','Key source reef uncertainty','plot_net_vs_predicted_KSR_selection.png'),
  plot = p_sel,
  width = 1600 / 225,
  height = 1200 / 300,
  dpi = 300
)

## Add a box showing the zoomed region (tagged image saved above)

# Create zoomed-in version, with reversed x-axis
p_zoom <- ggplot(score_table, aes(x = MeanRank, y = Net)) +
  geom_point(aes(fill = "All reefs"), shape = 21, color = "black", stroke = 0.3, size = 2, alpha = 0.7) +
  geom_vline(aes(xintercept = KSRzone, color = "KSR Zone Threshold", linetype = "KSR Zone Threshold"), linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_point(
    data = score_table[top_Opp, , drop = FALSE],
    aes(fill = "Top Opportunity reefs"),
    shape = 21,
    color = "black",
    stroke = 0.3,
    size = 2,
    alpha = 0.7
  ) +
  geom_point(
    data = score_table[top_Risk, , drop = FALSE],
    aes(fill = "Top Risk reefs"),
    shape = 21,
    color = "black",
    stroke = 0.3,
    size = 2,
    alpha = 0.7
  ) +
  scale_x_reverse(limits = c(500, 0), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1500, 500), expand = c(0, 0)) +
  scale_color_manual(values = c(
    "KSR Zone Threshold" = "red"
  )) +
  scale_fill_manual(values = c(
    "Top Opportunity reefs" = col_opp,
    "Top Risk reefs" = col_risk,
    "All reefs" = col_default
  )) +
  scale_linetype_manual(values = c("KSR Zone Threshold" = "dashed")) +
  labs(
    x = expression("Predicted KSR Rank (" * bar(R) * ")"),
    y = "Net Movement with Uncertainty",
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  theme_bw() +
  theme(legend.position = "right") +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        color = "black",
        stroke = 0.3,
        alpha = 0.7
      )
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(
        linetype = "dashed",
        linewidth = 0.6
      )
    ),
    linetype = "none"
  )

ggsave(
  filename = file.path('Outputs','Key source reef uncertainty','KSR_uncertainty_plot_zoom.png'),
  plot = p_zoom,
  width = 1600 / 225,
  height = 1200 / 300,
  dpi = 300
)



## ----------------------------------------------------------------------
## ---------- EXPORT THE TOP REEFS OF OPPORTUNITY AND RISK ----------
## Combine indices (and remove duplicates if needed)
selected_idx <- unique(c(top_Risk, top_Opp))

## Re-order selected reefs 
topIdx <- c(top_Risk, top_Opp)  

## Extract rows into new table
top20_table <- score_table[topIdx, ]

## save to file
write.csv(top20_table, file.path('Outputs','Key source reef uncertainty','o2_KSR_risky_and_opportunity_reefs.csv'), row.names = FALSE)
