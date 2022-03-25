#' Turnover rate to generation time
#' @param mu a rate of turnover in units turnover per unit time
#' @return generation time/doubling time in units of t
#' 

turnover_to_gen <- function(mu) {
  return(
    case_when(
      mu < 0 ~ 0,
      mu >= 0 ~ log(2) / mu
    )
  )
}
