#' Turnover Rate Calculation
#' @param a assimilation efficiency
#' @param FL 2F of label in at%
#' @param FT 2F of biomass at time t in at%
#' @param F0 2F of biomass at time 0 in at%
#' @param t incubation time
#' @return turnover rate, µ, in units of 't'


calculate_turnover <- function(a, FT, F0, FL, t) {
  
  # Calculate µ
  mu <- case_when(
    t == 0 ~ NA_real_,
    TRUE ~ - (1/t)*(log((FT - a*FL)/(F0 - a*FL)))
    )
  
  # Return µ if positive
  return(
    case_when(
      mu >= 0 ~ mu,
      mu < 0 ~ 0
      )
  )
}
