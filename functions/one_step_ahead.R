#function to evaluate performance of ARIMA model (produces season total forecasts only)
one_step_ahead<-function(series,
                            leave_yrs,
                            TY_ensemble,
                            covariates,
                            first_forecast_period
                          
                            ){
  
  start<-Sys.time()

  
  series<-series%>%
    ungroup()%>%
    dplyr::select(year,species,period,abundance,all_of(unique(unlist(covariates))))%>%
    filter(
      across(
        .cols = all_of(unique(unlist(covariates))),
        .fns = ~ !is.na(.x)
      )
    )
  exists1 =ifelse(is.na(series%>%dplyr::select(abundance)%>%tail(n=1)%>%pull()),1,0)
  
  library(doParallel)
  cl <- makeCluster(parallel::detectCores()-3)
  registerDoParallel(cl)
  forecasts_out<- foreach::foreach(i=1:leave_yrs,.combine =  'rbind',.packages=c("tidyverse","forecast")) %dopar% {
  source("functions/arima_forecast.r")
  for(c in 1:length(covariates)){
    
    
    # for(i in 1:leave_yrs){
      last_train_yr = max(series$year) - (leave_yrs-i+exists1)
      tdat<-series%>%
        filter(year <= (last_train_yr + 1))%>%
        mutate(train_test = ifelse(year > last_train_yr & period >= first_forecast_period, 1, 0),
        )

      xreg<-tdat%>%
        filter(train_test==0)%>%
        ungroup%>%
        dplyr::select(all_of(covariates[[c]]))%>%
        as.matrix()
      
      xreg_pred<-tdat%>%
        filter(train_test==1)%>%
        ungroup%>%
        dplyr::select(all_of(covariates[[c]]))%>%
        as.matrix()
      
      temp<-list(
        pred=NA,
       CI=tibble(`Lo 50`=NA,
                 `Hi 50`=NA,
                 `Lo 95`=NA,
                 `Hi 95`=NA,
                 year=last_train_yr+1,
                 period=1),
       arma=NA,
       aicc=NA
      )
      try(temp<-arima_forecast(tdat,xreg,xreg_pred,last_train_yr,first_forecast_period))
      pred<-temp$pred %>% tail(1)
      CI<-temp$CI


        
      tdat<-tdat%>% filter(train_test==1) %>% dplyr::select( c("year","period")) %>% 
                    bind_cols(data.frame(
                      predicted_abundance=pred,
                      arma=temp$arma,
                      aicc=temp$aicc))%>%
                    left_join(CI, by = c("year","period"))%>%
                    # dplyr::rename(predicted_abundance = pred)%>%
        mutate(model=as.character(c))
      
      if(c==1){forecasts = tdat
      }else{forecasts = forecasts %>% bind_rows(tdat)}
      
     
    }
    # forecasts<-forecasts%>%
    #   mutate(error = predicted_abundance-abundance,
    #          pct_error=scales::percent(error/abundance),
    #          model = as.character(c)
    #   )
    # if(c==1){
    #   tdat2 <- forecasts
    # }else{
    #   if(sum(is.na(forecasts$predicted_abundance))==0){
    #     tdat2 <- tdat2%>%
    #       bind_rows(forecasts)
    #   }
    # }
    return(forecasts )
  }
  # forecasts<-tdat2
  
  #do model averaging and stacking and calculate performance metrics
  # forecast_eval<-evaluate_forecasts_with_ensembles2(forecasts=forecasts,series=series,TY_ensemble=TY_ensemble,k=k,leave_yrs=leave_yrs)
  # 
  # return(forecast_eval)
  stopCluster(cl)
  
  
  print((Sys.time()-start))
  return(  forecasts_out)
}

