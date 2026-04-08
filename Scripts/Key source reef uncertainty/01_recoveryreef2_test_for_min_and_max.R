# con is a connectivity matrix of probabilities. cover is a vector of the
# coral cover in each cell (%). Alpha is coefficient to decide importance
# to place on the number of reefs involved where 1 is max vs 0 is ignore
# this and just focus on connectivity enrichment. Critic is the lowest
# level of coral cover (%) that doesn't need rescue.
# Here assume all reefs same size but allow for variation (sizes).
# out
# 1 reef number
# 2 unweighted importance as larval source (ignores reef neediness)
# 3 weighted importance as larval source for reef neediness
# 4 neediness (0-1)
# 5 number of reefs source connects to
# 6 number of reefs that are in need that a source connects to

# note that input reefsizes and coral cover for GBR are vertical so need to
# switch them. Also I'm going to square root transform the reef area so it
# better approximates perimeter rather than area

library(foreach)
library(doSNOW)
library(Matrix)

load('Data/keysourcereefs2023.RData')
coral_covers <- read.csv('Outputs/Use lookup tables/3.2_lower_and_upper_coral_cover_estimates.csv')

cover_min <- matrix(coral_covers$CoverLower_estimate, nrow=1) # if cover is zero then algorithm says source reefs aren't connected to any sinks so make minimum possible value slightly higher than zero.
cover_max <- matrix(coral_covers$CoverUpper_estimate, nrow=1) # if cover is zero then algorithm says source reefs aren't connected to any sinks so make minimum possible value slightly higher than zero.
cover_mean <- matrix(coral_covers$Mean, nrow=1)
scenario <- c("Min_Source_Score", "Max_Source_Score", "Mean_Source_Score")
scores <- list()

con <- MEAN_CONNECT  # Keep as sparse matrix for efficiency
# con <- con - diag(diag(con))  # This removes self-seeding. Using self
# seeding in importance score calcs can on rare occasions mean the min score is
# higher than the max score. Comment out for now, but re-introduce it if
# it becomes a problem later. 
reefsizes <- reefsizes
alpha <- 1
critic <- 20
region <- region


## Create for loop that does both min and max
start_time <- Sys.time()

# Set up parallel processing - use all available cores minus 1
num_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(min(num_cores, length(scenario)))  # Don't use more cores than scenarios
registerDoSNOW(cl)

# Progress bar
pb <- txtProgressBar(max = length(scenario), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

scores <- foreach(k = 1:length(scenario), .options.snow = opts, 
                  .packages = c('Matrix')) %dopar% {
    
# Variables need to be re-declared in each parallel worker
con_local <- con
reefsizes_local <- reefsizes
region_local <- region
cover_min_local <- cover_min
cover_max_local <- cover_max
cover_mean_local <- cover_mean
critic_local <- critic
    
row <- nrow(con_local)
col <- ncol(con_local)
out <- matrix(0, nrow=row, ncol=13)
out[,12] <- region_local[,1]

# calculate weights of each reef based on current cover
need <- matrix(0, nrow=1, ncol=row)

    #### ADD ANOTHER NESTED IF STATEMENT THAT CHECKS IF THE REEF IN
    #### 'NEED' IS A SINK FOR REEF i - IF SO THEN IT GETS MAX COVER AND IF
    #### NOT IT GETS MEAN COVER.
    for (i in 1:col) { # Calculate score for each reef.
        cover <- cover_mean_local    # reset every time we start reef i
        need  <- matrix(0, nrow=1, ncol=col)  # (and probably reset need too)

        # Get sinks for this source (sparse matrix optimization)
        sink_indices <- which(con_local[i,] > 0)

        for (j in 1:col) { # Modify coral cover for each sink of reef i so the score can be calculated.
            if (scenario[k] == "Min_Source_Score") {
                if (i == j) {
                    cover[1, j] <- cover_min_local[1, j] # if 'sink' reef is the source (reef i) then it gets minimum cover value (to minimise importance score).
                } else if (j %in% sink_indices) { # if reef j is a sink to reef i
                    cover[1, j] <- cover_max_local[1, j] # then the sink gets maximum cover value to minimise score of reef i.
                } else {
                    cover[1, j] <- cover_mean_local[1, j]  # if reef j is not a sink of reef i then it just gets its normal value.
                }

            } else if (scenario[k] == "Max_Source_Score") {
                if (i == j) {
                    cover[1, j] <- cover_max_local[1, j]
                } else if (j %in% sink_indices) {
                    cover[1, j] <- cover_min_local[1, j]
                } else {
                    cover[1, j] <- cover_mean_local[1, j]
                }

            } else if (scenario[k] == "Mean_Source_Score") {
                cover[1, j] <- cover_mean_local[1, j] # All reefs get mean coral cover score for mean KSmean
            }

            if (cover[1,j] >= critic_local) { #then calculate needy score for each sink (j) of reef i
                need[1,j] <- 0
            } else {
                need[1,j] <- (critic_local - cover[1,j]) / critic_local
            }
        }

        # Now for reef i, we should have the appropriate coral covers for
        # all sinks to calculate the minimum score of reef i.

        # create weighted contributions to connectivity
        unweightcon <- matrix(0, nrow=1, ncol=row)
        # For sparse matrix: multiply each column by corresponding cover*reefsize value
        larvalinputs <- con_local * (as.vector(cover) * as.vector(reefsizes_local))
        larvalinput_reef <- Matrix::colSums(larvalinputs)
        larvalflux_reef <- matrix(larvalinput_reef / as.vector(reefsizes_local), nrow=1) # convert to a flux based on size

        enrichunweight <- 0 # calc enrichment done by each source
        enrichweight <- 0
        # Get non-zero connections for reef i (sparse matrix optimization)
        conn_indices <- which(con_local[i,] > 0)
        for (j in conn_indices) {
            # add an if to ensure connectivity

            inputflux <- larvalinputs[i,j] / reefsizes_local[1,j] # convert to a flux
            propinputflux <- inputflux / larvalflux_reef[1,j] # express input as proprotion of total flux to that reef
            enrichunweight <- enrichunweight + propinputflux #create total
            enrichweight <- enrichweight + (propinputflux * need[1,j]) # weight by need at the sink
        }
        out[i,1] <- i
        out[i,2] <- enrichunweight
        out[i,3] <- enrichweight
        out[i,4] <- need[1,i]


        # identify number of reefs each source contacts and then the weighted
        # number which ignores those that don't need it

        numreefs <- 0
        numreefsneeded <- 0
        # Use sparse matrix row for efficiency
        larvae_row <- larvalinputs[i,]
        nonzero_idx <- which(larvae_row > 0)
        for (j in nonzero_idx) {
            if (j != i) {
                numreefs <- numreefs + 1
                if (cover[1,j] <= critic_local) {
                    numreefsneeded <- numreefsneeded + 1
                }
            }
        }
        out[i,5] <- numreefs
        out[i,6] <- numreefsneeded

    }

## Store output as minimum or maximum importance scores for each reef.

out

}

close(pb)
stopCluster(cl)

# Convert list to named list
names(scores) <- scenario

elapsed_time <- Sys.time() - start_time
print(elapsed_time)

## Checks - mean results to check min, mean and max are relatively correct (i.e. min is lowest + max is highest). 
# 2 unweighted importance as larval source (ignores reef neediness)
# 3 weighted importance as larval source for reef neediness
# 4 neediness (o-1)
# 5 number of reefs source connects to
# 6 number of reefs that are in need that a source connects to
colMeans(scores$Min_Source_Score[, 2:6])
colMeans(scores$Mean_Source_Score[, 2:6])
colMeans(scores$Max_Source_Score[, 2:6])

# Save simplified version of data that just has minimum, maximum and range
# (difference) in possible importance scores. Just use column 3 for now as importance score (summed weighted
# importance as a larval source for reef neediness).

score_table <- data.frame(
    Reef = scores$Mean_Source_Score[,1],
    Number_of_sinks = scores$Mean_Source_Score[,5],
    Mean = scores$Mean_Source_Score[, 3],
    Min = scores$Min_Source_Score[, 3],
    Max = scores$Max_Source_Score[, 3],
    Difference = scores$Max_Source_Score[, 3] - scores$Min_Source_Score[, 3]
)

# save to file
write.csv(score_table, "Outputs/Key source reef uncertainty/o1_KSR_scores_range_all_reefs.csv", row.names = FALSE)
