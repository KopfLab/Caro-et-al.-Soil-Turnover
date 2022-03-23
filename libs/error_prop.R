#' Propagation of uncertainty
#' @param a assimilation efficiency
#' @param F_L 2F of label at%
#' @param F_T 2F of biomass at time t at%
#' @param F_0 2F of biomass at time 0 at%
#' @param t incubation time
#' @param sF_0 uncertainty in F_0
#' @param sF_L uncertainty in F_L
#' @param sF_T uncertainty in F_T
#' @param sa uncertainty in a
#' @return the uncertainty in µ, turnover rate 'su'

calculate_sigma_mu <- function(a, F_L, F_T, F_0, t, sF_0, sF_L, sF_T, sa) {
  
  # the big uncertainty propagation equation :)
  
 
  # numerator
  num = sqrt(
    ((a*F_L - F_T)^2 * sF_0^2) +
      
    ((a*F_L - F_0)^2 * sF_T^2) +
      
    (a^2 * (F_0 - F_T)^2 * sF_L^2) +
      
    (F_L^2 * (F_0 - F_T)^2 * sa^2)
  )
  
  
  # denominator
  denom = t * (a*F_L - F_0) * (a*F_L - F_T)
  
  su = num/denom
  return(su)
}