
#' read a chromleon export file
#' @param file_path path to the excel file to be read
#' @return a list with two data frames for injection_details and integration_results
read_chromeleon_export_peaks <- function(file_path) {
  
  # read raw data
  raw_data <- suppressMessages(readxl::read_excel(file_path, sheet = "Peak Analysis"))
  names(raw_data) <- paste0("x", 1:ncol(raw_data))
  
  
  injection_details <- raw_data %>% 
    # remove rows that have no information
    filter(!is.na(x1)) %>% 
    # focus on rows between Injection Details and Chromatogram sections
    filter(row_number() > which(x1 == "Injection Details")[1] & row_number() < which(x1 == "Chromatogram")[1]) %>% 
    # remove columns that have no information in them at all
    { .[map_lgl(., ~!all(is.na(.x)))] }
  injection_details %>% knitr::kable()
  
  all_results <- raw_data %>% 
    filter(row_number() > which(x1 == "Peak Analysis")[1]) 
  
  integration_results <- 
    setNames(all_results[-c(1:3),], t(all_results[1,])[,1]) %>%
    mutate_at(vars(-`Peak Name`), function(x) suppressWarnings(as.numeric(x)))
  
  list(injection_details = injection_details, peak_analysis = peak_analysis)
}