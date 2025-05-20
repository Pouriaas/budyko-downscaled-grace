rm(list=ls())
library(ncdf4)
library(raster)
library(sf)
library(stars)
library(ggplot2)
library(hydroGOF)
library(GA)
library(lubridate)
library(gdata)
library(xlsx)

###we read the function in here
source(paste0("C:/R_for_Machine learning/bodyko/functions_codes/",
              "budyko_optimization_functions/with_deltaS/with_corroloatory.R"))

load(paste0("C:/R_for_Machine learning/bodyko/",
        "outputs/comparison_dataframe/with_Deltas/Grace/",
        "GRACE_REC_v03_JPL_ERA5_monthly_ensemble_mean/",
        "comparison_data_with_Era_land_ppt_pet_aet/",
        "comparison_yearly_bodyko_corrolatory.RData"))



  
  ###I selected  stations with the most bias to be removed

#omitted_station_codes<-c("35-009","21-223","21-105","41-907","41-111",
 #                        "44-035","41-083","21-261","12-027","21-131",
  #                       "16-003","44-105","41-149","41-143","21-273",
   #                      "34-011","41-117","41-039","41-165","16-051",
    #                     "21-275","22-045","18-061","41-115","41-159",
                         #          #STREAMFLOW BASED
        #                 "16-021","12-043","12-085","15-015","46-011",
         #                "41-109","41-161","41-101","41-041","18-089",
                         #worst KGE
           #              "13-005","31-013")




###I selected 33 stations with the most bias to be removed

####based on KGE of each statiion
omitted_station_codes<-c("35-009","21-223","21-105","41-907","41-111",
                         "44-035","41-083","21-261","16-003","12-027",
                         "21-131","21-273","34-011","44-105","41-039",
                         "18-061","22-045","41-149","18-089","41-143",
                         "41-117","21-275","33-015","41-165","15-015",
                         #STREAMFLOW BASED
                         "16-051","16-021","12-043","12-085","16-023",
                         "13-005","34-003","41-115","46-011","41-159",

      #                      "41-109","41-161","41-101","41-041","18-089",
                         #worst KGE
                         "47-029","31-013","41-099")



omitted_station_rows<-c()
  for (i in 1:length(omitted_station_codes)) {
    
    omitted_station_rows<-c(omitted_station_rows,
                            which(comparison_datafram_yearly$Code==omitted_station_codes[i]))
  }
  
 data_comparison_dataframe_yearly<-comparison_datafram_yearly[-c(omitted_station_rows),]
 
# write.xlsx(data_comparison_dataframe_yearly,
  #    "C:/R_for_Machine learning/bodyko/outputs/data_for_Dr/excel_files/GRACE_ERA5/Era_land_Pre_aev_pot/comparison.xlsx")

data_comparison_dataframe_yearly<-comparison_datafram_yearly
  
  dec_related_variable<-c("NDVI","DEM","curve number","lon","lat")
  colnum_min_max<-length(dec_related_variable)
  
  max_comparison<-data.frame(matrix(nrow=1,ncol=colnum_min_max))
  min_comparison<-data.frame(matrix(nrow=1,ncol=colnum_min_max))
  
  for (i in 1:colnum_min_max) {
    colnames(max_comparison)[i]<-dec_related_variable[i]
    colnames(min_comparison)[i]<-dec_related_variable[i]
    
    colnum_min_max<-which(colnames(data_comparison_dataframe_yearly)==dec_related_variable[i])
    
    max_comparison[i]<-max(data_comparison_dataframe_yearly[,colnum_min_max])
    min_comparison[i]<-min(data_comparison_dataframe_yearly[,colnum_min_max])
    
    data_comparison_dataframe_yearly[,colnum_min_max]<-
      TSPred::minmax(data_comparison_dataframe_yearly[,colnum_min_max])
    
  }
  
  
  train_percent<<-0.7
  test_percent<<-1-train_percent
  
  shomarande_train<-sample.int(nrow(data_comparison_dataframe_yearly),
                               floor(train_percent*nrow(data_comparison_dataframe_yearly)))
  
  shomarande_test<-(1:nrow(data_comparison_dataframe_yearly))[-c(shomarande_train)]
  
  ##train and test data are gathered in here
  train_data<-data_comparison_dataframe_yearly[shomarande_train,]
  test_data<-data_comparison_dataframe_yearly[shomarande_test,]
  
  ###relative importance
  relative_importance<<-3
  
  
  kz<<-1
  
  dec_numbers<-6
  dec_ob_number<-dec_numbers+1
  Result_GA<<-data.frame(matrix(nrow=10^6,ncol=dec_ob_number))
  
  NDVI<<-train_data$NDVI
  Long<<-train_data$lon
  Lat<<-train_data$lat
  Alt<<-train_data$DEM
  curve_number<<-train_data$`curve number`
  
  #####delta_storage is entered here as effective rainfall
  Rainfall<<-train_data$Rain_prod_vol-train_data$delta_storage
  pot_Evaporation<<-train_data$Potential_Eva_prod_vol
  actual_Evaporation<<-train_data$Eva_prod_vol
  Streamflow<<-train_data$hydr_vol
  #deltaS<<-train_data$delta_storage
  
  opt<-function(dec)
  {
    ###dec1=alpha
    ###dec2=v
    
    ##parameters:alhpha
    ##A,B,C,D,E
    #  alpha<-dec[1]
    shomareshgar<-min(length(Rainfall),length(NDVI),length(Lat),
                      length(Long),length(curve_number),length(Alt),
                      length(pot_Evaporation))
    
    simulated_Evapotranspiration<-c()
    simulated_Streamflow<-c()
    for (i in 1:shomareshgar) {
      
      
      simulated_Evapotranspiration[i]<-ifelse(Rainfall[i]<=0,0,Fu(dec[1],dec[2],dec[3],dec[4],dec[5],
                                                                     dec[6],
                                                                     NDVI[i],Long[i],Lat[i],Alt[i],curve_number[i],
                                                                     Rainfall[i],pot_Evaporation[i])) 
      
      simulated_Streamflow[i]<-Rainfall[i]-ifelse(Rainfall[i]<=0,0,Fu(dec[1],dec[2],dec[3],dec[4],dec[5],
                                                                      dec[6],
                                                                      NDVI[i],Long[i],Lat[i],Alt[i],curve_number[i],
                                                                      Rainfall[i],pot_Evaporation[i]))
      
    }
    
  #  objective<-(-1*sum((Rainfall-simulated_Evapotranspiration-Streamflow)^2))
    
  #  objective<-1/(sum(abs(Streamflow-simulated_Streamflow))*relative_importance+
   #               sum(abs(actual_Evaporation-simulated_Evapotranspiration)))
    
    objective<-(hydroGOF::KGE(simulated_Streamflow,Streamflow)-1)*relative_importance+
      (hydroGOF::KGE(actual_Evaporation,simulated_Evapotranspiration)-1)
    
    Result_GA[kz,]<<-c(dec,objective)
    kz<<-kz+1
    
    return(objective)
    #  Result_GA[kz,]<-c(dec,objective)
  }
  
  lower <- rep(-20,dec_numbers)
  upper <- c(rep(20,(dec_numbers)))
  
  iteration<-200
  GA <- ga(type = "real-valued", 
           fitness =opt,
           lower = lower, upper = upper, 
           popSize = 200, 
           maxiter = iteration, run =iteration,pmutation = 0.15)
  
 # save.image("C:/R_for_Machine learning/bodyko/outputs/Simulation/with_deltas/yearly/Era_land_P_ERA_pot/Fu.Rdata")

  
  #####finding dec varaibles for each generation
  
  Result_GA<-Result_GA[complete.cases(Result_GA),]
  
  shomarande_solgen<-matrix(nrow=iteration,ncol=1)
  
  dec_gen<-matrix(nrow=iteration,ncol=2)
  
  for (i in 1:iteration) {
    
    shomarande_solgen[i,1]<-which.min(abs(GA@summary[i,1]-Result_GA[,ncol(Result_GA)]))
    
  }
  dec_gen<-Result_GA[shomarande_solgen,1:(dec_numbers)]
  
  ##function for test
  
  NDVI<<-test_data$NDVI
  Long<<-test_data$lon
  Lat<<-test_data$lat
  Alt<<-test_data$DEM
  curve_number<<-test_data$`curve number`
  #####delta_storage is entered here as effective rainfall
  Rainfall<<-test_data$Rain_prod_vol-test_data$delta_storage
  pot_Evaporation<<-test_data$Potential_Eva_prod_vol
  actual_Evaporation<<-test_data$Eva_prod_vol
  Streamflow<<-test_data$hydr_vol
  
  opt<-function(dec)
  {
    
    shomareshgar<-min(length(Rainfall),length(NDVI),length(Lat),
                      length(Long),length(curve_number),length(Alt),
                      length(pot_Evaporation))
    
    simulated_Evapotranspiration<-c()
    simulated_Streamflow<-c()
    for (i in 1:shomareshgar) {
      
      simulated_Evapotranspiration[i]<-ifelse(Rainfall[i]<=0,0,Fu(dec[1],dec[2],dec[3],dec[4],dec[5],
                                                                     dec[6],
                                                                     NDVI[i],Long[i],Lat[i],Alt[i],curve_number[i],
                                                                     Rainfall[i],pot_Evaporation[i])) 
      
      simulated_Streamflow[i]<-Rainfall[i]-ifelse(Rainfall[i]<=0,0,Fu(dec[1],dec[2],dec[3],dec[4],dec[5],
                                                                      dec[6],
                                                                      NDVI[i],Long[i],Lat[i],Alt[i],curve_number[i],
                                                                      Rainfall[i],pot_Evaporation[i]))
    }
    
    objective<-(hydroGOF::KGE(simulated_Streamflow,Streamflow)-1)*relative_importance+
     (hydroGOF::KGE(actual_Evaporation,simulated_Evapotranspiration)-1)
    
 #   objective<-1/(sum(abs(Streamflow-simulated_Streamflow))*relative_importance+
  #                  sum(abs(actual_Evaporation-simulated_Evapotranspiration)))
    
    return(objective)
    #  Result_GA[kz,]<-c(dec,objective)
  }
  
  ####finding test results
  
  test_result<-data.frame(matrix(nrow = length(shomarande_solgen),ncol=1))
  for (i in 1:length(shomarande_solgen)) {
    dec<-c()
    for (j in 1:dec_numbers) {
      dec[j]<-dec_gen[i,j]
    }
    test_result[i,1]<-opt(dec)
  }
  
  ####
  ###finding best answer
  Result_train_test<-data.frame(matrix(nrow = length(shomarande_solgen),ncol=3))
  #Result_train_test<-data.frame(matrix(nrow = length(shomarande_solgen),ncol=2))
  
  colnames(Result_train_test)<-c("generation","train","test")
  #olnames(Result_train_test)<-c("generation","train")
  
  Result_train_test[,1]<-1:iteration
  
  Result_train_test[,2]<-GA@summary[,1]
  
  Result_train_test[,3]<-test_result[,1]
  
  ####
  
  #####finally we find decision variables here
  dec<-as.numeric(dec_gen[which.max(Result_train_test[,3]),])
  
  NDVI<<-data_comparison_dataframe_yearly$NDVI
  Long<<-data_comparison_dataframe_yearly$lon
  Lat<<-data_comparison_dataframe_yearly$lat
  Alt<<-data_comparison_dataframe_yearly$DEM
  curve_number<<-data_comparison_dataframe_yearly$`curve number`
  
  Rainfall<<-data_comparison_dataframe_yearly$Rain_prod_vol-
  data_comparison_dataframe_yearly$delta_storage
  pot_Evaporation<<-data_comparison_dataframe_yearly$Potential_Eva_prod_vol
  actual_Evaporation<<-data_comparison_dataframe_yearly$Eva_prod_vol
  #deltaS<<-data_comparison_dataframe_yearly$delta_storage
  opt_final<-function(dec)
  {
    #  alpha<-dec[1]
    
    shomareshgar<-min(length(Rainfall),length(NDVI),length(Lat),
                      length(Long),length(curve_number),length(Alt),
                      length(pot_Evaporation))
    
    simulated_streamflow<-c()
    simulated_Evapotranspiration<-c()
    for (i in 1:shomareshgar) {
      
      simulated_streamflow[i]<-Rainfall[i]-ifelse(Rainfall[i]<=0,0,Fu(dec[1],dec[2],dec[3],dec[4],dec[5],
                                                                         dec[6],
                                                                         NDVI[i],Long[i],Lat[i],Alt[i],curve_number[i],
                                                                         Rainfall[i],pot_Evaporation[i]))
      
      simulated_Evapotranspiration[i]<-ifelse(Rainfall[i]<=0,0,Fu(dec[1],dec[2],dec[3],dec[4],dec[5],
                                                                  dec[6],
                                                                  NDVI[i],Long[i],Lat[i],Alt[i],curve_number[i],
                                                                  Rainfall[i],pot_Evaporation[i])) 
    }
    
    return(cbind(simulated_streamflow,simulated_Evapotranspiration))
    #  Result_GA[kz,]<-c(dec,objective)
  }
  ###simulated streamflow and evaporation are gathered here
  
  simulated_result<-data.frame(opt_final(dec))
  
  
  ####finding simulated streamflow
  Streamflow_final<-data.frame(matrix(ncol=4,nrow=nrow(data_comparison_dataframe_yearly)))
  row.names(Streamflow_final)<-row.names(data_comparison_dataframe_yearly)
  
  colnames(Streamflow_final)<-c("Code","year","simulated","observed")
  
  Streamflow_final$Code<-data_comparison_dataframe_yearly$Code
  Streamflow_final$year<-data_comparison_dataframe_yearly$year
  ###simulated and observed
  Streamflow_final$simulated<-simulated_result$simulated_streamflow
  
  Streamflow_final$observed<-data_comparison_dataframe_yearly$hydr_vol
  
  KGE_streamflow=hydroGOF::KGE(Streamflow_final[,"simulated"],Streamflow_final[,"observed"])
  
  print(KGE_streamflow)
  
  ####finding simulated evapotranspiration
  
  evapotranspiration_final<-data.frame(matrix(ncol=4,nrow=nrow(data_comparison_dataframe_yearly)))
  row.names(evapotranspiration_final)<-row.names(data_comparison_dataframe_yearly)
  
  colnames(evapotranspiration_final)<-c("Code","year","simulated","observed")
  
  evapotranspiration_final$Code<-data_comparison_dataframe_yearly$Code
  evapotranspiration_final$year<-data_comparison_dataframe_yearly$year
  ###simulated and observed
  evapotranspiration_final$simulated<-simulated_result$simulated_Evapotranspiration
  
  evapotranspiration_final$observed<-data_comparison_dataframe_yearly$Eva_prod_vol
  
  KGE_evapotranspiration=hydroGOF::KGE(evapotranspiration_final[,"simulated"],evapotranspiration_final[,"observed"])
  
  print(KGE_evapotranspiration)  
  
######order
  delta_fitness_list<-data.frame(matrix(nrow=length(unique(Streamflow_final$Code)),ncol=2))
  delta_fitness_list[,1]<-unique(Streamflow_final$Code)
  row.names(delta_fitness_list)<-unique(Streamflow_final$Code)
  colnames(delta_fitness_list)<-c("Code","delta_fitness")
  
  for (i in 1:length(unique(Streamflow_final$Code))) {
    #####streamflow
    shomarande_KGE_omitted_streamflow<-
      which(Streamflow_final$Code==unique(Streamflow_final$Code)[i])
    
    Streamflow_final_KGE_omitted<-Streamflow_final[-c(shomarande_KGE_omitted_streamflow),]
    
    KGE_streamflow_omitted<-hydroGOF::KGE(Streamflow_final_KGE_omitted[,"simulated"],
                                          Streamflow_final_KGE_omitted[,"observed"])
    #####evapotranspiration
    shomarande_KGE_omitted_evapotranspiration<-
      which(evapotranspiration_final$Code==unique(evapotranspiration_final$Code)[i])
    
    evapotranspiration_final_KGE_omitted<-evapotranspiration_final[-c(shomarande_KGE_omitted_evapotranspiration),]
    
    KGE_evapotranspiration_omitted<-hydroGOF::KGE(evapotranspiration_final_KGE_omitted[,"simulated"],
                                                  evapotranspiration_final_KGE_omitted[,"observed"])
    
    
    ####findingfitness
#   delta_fitness_list[i,2]<-((KGE_evapotranspiration_omitted-1)+relative_importance*(KGE_streamflow_omitted-1))-
  #    ((KGE_evapotranspiration-1)+relative_importance*(KGE_streamflow-1))
    
    delta_fitness_list[i,2]<-(KGE_streamflow_omitted-1)-(KGE_streamflow-1)
    
  }
  
  delta_fitness_list<-delta_fitness_list[order(-delta_fitness_list$delta_fitness),]
  
  KGE_vector<-data.frame(matrix(ncol=2,nrow=length(unique(Streamflow_final$Code))))
  for (i in 1:length(unique(Streamflow_final$Code))){
    row_KGE<-which(Streamflow_final$Code==
                     unique(Streamflow_final$Code)[i])
    
    KGE_vector[i,2]<-  hydroGOF::KGE(Streamflow_final$simulated[row_KGE],
                                     Streamflow_final$observed[row_KGE])
  }
  KGE_vector[,1]<-unique(Streamflow_final$Code)
  KGE_vector<-KGE_vector[order(KGE_vector[,2]),]
  
#  save(list = c("max_comparison",
  #              "min_comparison",
  #               "dec"),
   #   file = "C:/R_for_Machine learning/bodyko/outputs/bodyko_application/dec/GRACE/GRACE_REC_v03_JPL_ERA5_monthly_ensemble_mean/dec_max_min_Era_land_ppt_pev_aev_Fu2.Rdata")
  
# load("C:/R_for_Machine learning/bodyko/outputs/bodyko_application/dec/GRACE/GRACE_REC_v03_JPL_ERA5_monthly_ensemble_mean/dec_max_min_Era_land_ppt_pev_aev_Fu.Rdata")
  
# load("C:/R_for_Machine learning/bodyko/outputs/bodyko_application/dec/GRACE/GRACE_REC_v03_JPL_ERA5_monthly_ensemble_mean/dec_max_min_Era_land_ppt_pev_aev_Fu.Rdata")
 
# load("C:/R_for_Machine learning/bodyko/outputs/Simulation/with_deltas/GRACE/GRACE_REC_v03_JPL_ERA5_monthly_ensemble_mean/Era_land/Fu2.Rdata")
  # save.image("C:/R_for_Machine learning/bodyko/outputs/Simulation/with_deltas/GRACE/GRACE_REC_v03_JPL_ERA5_monthly_ensemble_mean/Era_land/Fu2.Rdata")
  
 
 
 # load("C:/R_for_Machine learning/bodyko/outputs/Simulation/with_deltas/GRACE/GRACE_REC_v03_JPL_MSWEP_monthly_ensemble_mean/Era_land/Fu.Rdata")
  
 # load("C:/R_for_Machine learning/bodyko/outputs/Simulation/with_deltas/yearly/Era_land_P_ERA_land_pot/Fu.Rdata")
  
  
  