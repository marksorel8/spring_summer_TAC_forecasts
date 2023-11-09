make_dat<-function(file_path=NULL,dat1,redo=TRUE){
  
  
  if(!is.null(file_path)&file.exists(file_path)&!redo){
    return(read_csv(file_path))
    
  }else{
    covariates<-  c(#new
      "lag1_log_Jack"
      ,"lag1_log_age4"
      ,"lag4_log_adults"
      ,"lag5_log_adults"
      ,"lag1_log_SAR1" #snake river sockeye
      ,"lag1_log_SAR2" #snake river sockeye
      ,"lag1_NPGO"
      ,"lag1_PDO"
      ,"lag2_NPGO"
      ,"lag2_PDO"
      ,"WSST_A"
      #NOAA Ocean Indicators
      ,"lag2_PC1"
      ,"lag2_PC2"
      ,"lag2_sp_phys_trans"
      ,"pink_ind"
    ) 
    
    
    Yrlist<-data.frame(year=c(min(dat1$year):(max(dat1$year)+1)))
    
    dat1%<>%
      right_join(Yrlist)%>%
      arrange(year)
    #=========================================================
    #get PDO data
    #=========================================================
    PDO<-read_table("https://psl.noaa.gov/pdo/data/pdo.timeseries.ersstv5.csv",skip=1,col_names=F,comment="#")%>%
      dplyr::rename(Date=X1,PDO=X2)%>%
      filter(!PDO < -99)%>%
      mutate(Date=as.Date(Date),Month=month(Date),Year=as.integer(year(Date)))%>%
      group_by(Year)%>%
      add_tally()%>%
      #filter(!Month>6)%>% #use only spring (Jan-June) NPGO
      #filter(!n < 12)%>% #use only complete years
      group_by(Year)%>%
      dplyr::rename(year=Year)%>%
      dplyr::summarise(PDO=mean(PDO))%>%
      
      dplyr::select(year,PDO)
    #=========================================================
    #get NPGO data
    #=========================================================
    NPGO<-read_table("http://www.o3d.org/npgo/npgo.php",skip=29,col_names=F,comment="#")%>%
      filter(!is.na(X2))%>%
      dplyr::rename(Year=X1,Month=X2,NPGO=X3)%>%
      mutate(Year=as.integer(as.numeric(Year)))%>%
      group_by(Year)%>%
      add_tally()%>%
      #filter(!Month>6)%>% #use only spring (Jan-June) NPGO
      #filter(!n < 12)%>% #use only complete years
      group_by(Year)%>%
      dplyr::summarise(NPGO=mean(NPGO))%>%
      dplyr::rename(year=Year)%>%
      
      dplyr::select(year,NPGO)
    #=========================================================
    #get NOAA indicator data, wrangle into usable format, plot
    #=========================================================
    indicators<-read_csv("https://www.fisheries.noaa.gov/s3//2022-12/OEI-spotlight-cvs-2022-NWFSC.csv",skip=1)%>%
      filter(!is.na(`Ecosystem Indicators`))%>%
      pivot_longer(names_to = "Year",
                   cols=c(starts_with("1"),starts_with("2")),
                   values_to = "value")%>%
      pivot_wider(names_from=`Ecosystem Indicators`,values_from=value)%>%
      mutate(year=as.integer(Year))%>%
      dplyr::select(-Year)
    #=========================================
    # Get ENSO index
    #=========================================
    enso<-read_table("https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/nino34.long.anom.data",skip=1,col_names = F)%>%
      as_tibble()%>%
      setNames(c("year",1:12))%>%
      # filter(year%in%as.character(yr_start:yr_end))%>%
      mutate(across(everything(),~as.numeric(.)))%>%
      pivot_longer(names_to = "month",values_to = "Nino3.4",cols=c(!year))%>%
      filter(Nino3.4>-10 & Nino3.4 <10)%>%
      filter(month%in%c(10:12))%>%
      group_by(year)%>%
      summarise(fall_Nino3.4=mean(Nino3.4))%>%
      right_join(Yrlist)%>%
      mutate(lag2_fall_Nino3.4=lag(fall_Nino3.4,2))%>%
      dplyr::select(-fall_Nino3.4)
    
    #===========================================================================================================
    #Get ERSST data: get_ersst_v5_data takes A LONG TIME (1 hr) vs. get_ersst_v5_data_V2 which is much quicker!
    #==========================================================================================================
    sstdat<-get_ersst_v5_data_V2(years=c(min(Yrlist$year):max(Yrlist$year)),
                                 data.dir=#"https://stateofwa.sharepoint.com/:u:/r/sites/DFW-TeamWDFWOPIModelReview/Shared%20Documents/General/sst.mnmean.nc?csf=1&web=1&e=mhx1bc",
                                   "C:\\Users\\sorelmhs\\Washington State Executive Branch Agencies\\DFW-Team WDFW OPI Model Review - General" ,
                                 ncfilename="sst.mnmean.nc",
                                 latrange=c(44,50),
                                 lonrange=c(-125,-120)
    )
    
    ssta<-sstdat%>%
      dplyr::select(year,month,resid)%>%
      mutate(month=as.numeric(month))%>%
      mutate(year=ifelse(month>9,year+1,year))%>%
      group_by(year)%>%
      summarise(WSST_A = mean(resid[month>9|month<11]),
                SSST_A = mean(resid[month<=9|month>=4]),
      )
  
    PIT<-read_csv("https://www.cbr.washington.edu/dart/cs/php/rpt/pit_sar_esu.php?queryType=year&proj=BON&esu_type=SR_Sock&rt=A&age=all&grouptype=basin&csvOnly=1") %>% 
      mutate(across(year,as.numeric)) %>% 
      filter(year!=year(Sys.Date()))%>%
      mutate(OutmigrationYear=year,Year=OutmigrationYear+2)
    
    PIT<-PIT%>%bind_cols(data.frame(SAR1=gam(cbind(ocean1Count,juvCount-ocean1Count)~s(OutmigrationYear,k=(dim(PIT)[1]),m=1,bs="ps"),family=binomial,data=PIT)$fitted))%>%
      bind_cols(data.frame(SAR2=c(gam(cbind(ocean2Count,juvCount-ocean2Count)~s(OutmigrationYear,k=(dim(PIT)[1]-1),m=1,bs="ps"),family=binomial,data=PIT)$fitted,NA)))%>%
      mutate(lag1_log_SAR1 = log(SAR1),lag1_log_SAR2=lag(log(SAR2),1))%>%
      dplyr::select(year=Year,lag1_log_SAR1,lag1_log_SAR2)
    #================================================================
    dat<-dat1%>%
      ungroup %>% 
      # left_join(Yrlist)%>%
      left_join(PDO)%>%
      left_join(NPGO)%>%
      left_join(indicators)%>%
      left_join(enso)%>%
      left_join(ssta)%>%
      left_join(PIT)%>%
      # left_join(OCN)%>%
      # dplyr::rename(lag1_JackOPI = lagJackOPI,
      #               lag1_SmAdj = lagSmAdj
      #               )%>%
      mutate(species = "fish",
             period = 1,
             lag1_log_Jack = lag(log(Jack)),
             lag1_log_age4=lag(log(Age4)),
             lag4_log_adults=lag(log(abundance),4),
             lag5_log_adults=lag(log(abundance),5),
             lag2_sp_phys_trans = lag(`Physical Spring Trans.\nUI based (day of year)`,2),
             lag2_PC1 = lag(scale(`Principal Component scores (PC1)`),2),
             lag2_PC2 = lag(scale(`Principal Component scores (PC2)`),2),
             lag1_NPGO = lag(c(scale(NPGO))),
             lag2_NPGO = lag(c(scale(NPGO)),2),
             lag1_PDO = lag(c(scale(PDO))),
             lag2_PDO = lag(c(scale(PDO)),2),
             pink_ind = ifelse(year>1999 & year%%2==0,0,1)
      )%>%
      ungroup()%>%
      dplyr::select(year,species,period,abundance,all_of(unique(unlist(covariates))))%>%
      mutate(across(all_of(unique(unlist(covariates))),\(x) c(scale(x))),
             across(all_of(unique(unlist(covariates))),\(x) replace_na(x,0))) %>% 
      filter(
        !is.na(abundance)|year==yr_end
        # across(
        #   .cols = everything(),
        #   .fns = ~ !is.na(.x)
        # )
      )#%>%mutate(across(!c(year,species,period,abundance),boxcox_scale(.)))
    
    
    
    write_csv(dat,file=file_path)  
    return(dat)
  }
  
  
}


