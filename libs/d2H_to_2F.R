#' Delta 2H to 2F
#' @param d2H deuterium isotopic enrichment in delta notation relative to VSMOW
#' @return fractional abundance of deuterium in atom percent
#' 

d2H_to_2F <- function(d2H) {
  # Define VSMOW in ratio space
  R_VSMOW = 0.0001557643
  

  # Convert delta to ratio, remove permil units from d2H
  R_sample = ((d2H/1000) + 1) * R_VSMOW
  
  # Convert ratio to fractional abundance
  F_sample = R_sample/(1 + R_sample)
  
  # Return fractional abundance in atom percent (*100)
  # We do this because our IRMS reports fractional abundance in atom percent!
  return(F_sample*100)
}