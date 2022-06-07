#' Rename fatty acids 
#' @param df a data frame containing a "compound" column containing FA names
#' @return data frame with renamed FAs

rename_FAs <- function(df) {
  return(
    df %>% 
      mutate(compound = case_when(
        compound == "16:1 cis 9" ~ "16:1ω7c",
        compound == "16:1 trans 9" ~ "16:1:ω5c",
        compound == "18:2 9,12" ~ "18:2",
        compound == "18:1 cis 9" ~ "18:1ω9",
        compound == "18:1 trans 9" ~ "18:1ω7",
        TRUE ~ compound
      ))
  )
}