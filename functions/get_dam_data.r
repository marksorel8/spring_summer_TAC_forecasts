#function to get DART dam counts
get_dam_data<-function(yr_start,yr_end,dam,species){
  for(yr in yr_start : yr_end){
    if(yr==yr_start){
      series<-tibble(NULL)
    }
    series<-series%>%
      bind_rows(read_csv(paste0(
        "https://www.cbr.washington.edu/dart/cs/php/rpt/adult_daily.php?sc=1&outputFormat=csv&year=",
        yr,
        "&proj=",
        dam,
        "&span=no&startdate=",
        1,
        "%2F",
        1,
        "&enddate=",
        12,
        "%2F",
        31,
        "&run=&syear=",
        yr,
        "&eyear=",
        yr
      )
      )%>%
        dplyr::select(Date,as.name(species))%>%
        filter(!is.na(Date)))
  }
  series<-series%>%
    dplyr::rename(abundance=as.name(species), date = Date)%>%
    mutate(year = year(date),
           abundance = ifelse(is.na(abundance)|abundance < 0,0,abundance),
           species = species
    )
}
