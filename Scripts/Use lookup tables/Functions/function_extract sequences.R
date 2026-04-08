
# Get sequence of recovery periods (i.e. sequence of trajectories to sample).
extract_sequences <- function(regime) {
  years <- regime$Year
  cyclone <- regime$Cyclone
  dhw <- regime$DHW_cat
  
  cyclone_years <- which(cyclone != 0)
  dhw_years <- which(dhw != 0)
  sequences <- character()
  
  if (length(cyclone_years) == 0 & length(dhw_years) == 0) {
    # Entire period is undisturbed
    sequences <- paste0("Cyclone_Value_0_Recovery_", length(years)) # Use cyclone because only cyclone lookup table has no disturbance values.
  } else {
    # Add undisturbed period before first disturbance (if any)
    first_years <- c(
      if (length(cyclone_years) > 0) cyclone_years[1],
      if (length(dhw_years) > 0) dhw_years[1]
    )
    
    if (length(first_years) > 0 && min(first_years) > 1) {
      sequences <- c(sequences, paste0("Cyclone_Value_0_Recovery_", min(first_years) - 1))
    }
    
    # 1. merge & sort the two event‐year vectors
    event_years <- sort(c(cyclone_years, dhw_years))
    
    # 2. loop over them exactly as before
    for (i in seq_along(event_years)) {
      yr <- event_years[i]
      
      if (yr %in% cyclone_years) {
        typ   <- "Cyclone"
        value <- cyclone[yr]
      } else {
        typ   <- "DHW"
        value <- dhw[yr]
      }
      
      recovery <- if (i < length(event_years)) {
        event_years[i + 1] - yr - 1
      } else {
        length(years) - yr
      }
      
      sequences <- c(
        sequences,
        paste0(typ, "_Value_", value, "_Recovery_", recovery)
      )
    }
    
  }
  
  return(sequences)
}