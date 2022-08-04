# Get dayflow from EXPORTS
# Adapted from https://github.com/goertler/inundation/blob/main/R/get_dayflow.R

get_dayflow <- function(){
  
  # get metadata
  m <- jsonlite::fromJSON("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536.jsonld")
  
  file_table <- m$`@graph`
  
  file_table <- subset(file_table, `dct:format` == "CSV")
  
  urls <- grep("results", file_table$`dcat:accessURL`$`@id`, value = TRUE)
  
  
  # read in the data
  col_types <- readr::cols(.default = readr::col_character())
  dat <- lapply(urls, readr::read_csv, col_types=col_types, show_col_types = FALSE, progress = FALSE)
  suppressWarnings(dat <- lapply(dat, function(x){
    if (is.null(x$EXPORTS)){
      x$EXPORTS <- NA
    }
    return(x[, c("Date", "EXPORTS")])
  }))
  
  
  # bind data
  dayflow <- do.call(rbind, dat)
  
  # rename columns
  dayflow$Date <- lubridate::parse_date_time(dayflow$Date, orders = c("mdy", "ymd", "dmy"))
  dayflow$EXPORTS <- as.numeric(dayflow$EXPORTS)
  
  # clean names
  dayflow <- janitor::clean_names(dayflow)
  
  # remove duplicates
  i <- which(!duplicated(dayflow))
  dayflow <- dayflow[i, ]
  
  
  # write out
  #utils::write.csv(dayflow, file.path(rappdirs::user_cache_dir("data_raw"), "exports_dayflow.csv"), row.names = FALSE)
  write.csv(dayflow, file.path("data_raw", "exports_dayflow.csv"), row.names = FALSE)
  return(dayflow)
  
}
