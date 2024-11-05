rm(list=ls())

# Ordi perso ou Cluster Crest
path <- "C:/Users/milen/Desktop/DisentanglingAP/"
path <- "//abra/current/CRT_PRJ_AIR_POLLUTION/Disentangling"
setwd(path)

path_code <- "W:/Disentangling"

library(stargazer)
library(lubridate)
library(lfe) 
library(haven)
library(hdm)
library(ggplot2)
library(data.table)
library(missMDA)

if(!"Input Data" %in% dir()){dir.create("Input Data")} 
if(!"Data" %in% dir()){dir.create("Data")}
if(!"Sorties" %in% dir()){dir.create("Sorties")}

# Overwrites hdm package functions with own routines, and additional routines. 
source(paste0(path_code,"/R/lambdaCalculationClusteredOption.R"))
source(paste0(path_code,"/R/selected_instruments.R"))
source(paste0(path_code,"/R/rlasso_unpenalized_controls_ClusteredOption.R"))
source(paste0(path_code,"/R/rlasso_reintroduce_controls_ClusteredOption.R"))
source(paste0(path_code,"/R/twoSLS_on_selected.R"))
source(paste0(path_code,"/R/scale_penalty.R"))
source(paste0(path_code,"/R/reg_FE_clustered.R"))

source(paste0(path_code,"/R/purge_FE.R"))
source(paste0(path_code,"/R/interact_city_continuous.R"))
source(paste0(path_code,"/R/get_angle_cat.R"))


source(paste0(path_code,"/R/MHT_corrections.R"))
source(paste0(path_code,"/R/randomized_inference.R"))

####################################################################################################################################
########## =========================== Estimation Sample For Multi Pollutant Analysis  ============================== ##############
####################################################################################################################################

load("data/base_IV_purged_FE.Rdata")

setDT(data)

data$month_year<-factor(substr(data$date,1,7))
data$month_year_ville<-factor(paste0(data$month_year,data$ville)) 
data$FEP_ff_2 <-data$FEP_ff^2
data$FEP_temperature_2 <-data$FEP_temperature^2
data$FEP_precipitations_2 <-data$FEP_precipitations^2
data$FEP_humidite_2 <-data$FEP_humidite^2
data$humidite_2 <-data$humidite^2
data[,all5in:=!is.na(PM2_5_mean_city) & !is.na(CO_mean_city) & !is.na(SO2_mean_city) & !is.na(NO2_mean_city) & !is.na(O3_mean_city)]
data$FEP_temperature_min_2 <- data$FEP_temperature_min^2
data$FEP_temperature_max_2 <- data$FEP_temperature_max^2
data$FEP_precipitations_min_2 <- data$FEP_precipitations_min^2
data$FEP_precipitations_max_2 <- data$FEP_precipitations_max^2
data$FEP_ensoleil_glo_2 <- data$FEP_ensoleil_glo^2
data[,c('scaled_ipblh','scaled_invth') := lapply(.SD, FUN = scale), .SDcols = c('inv_hcl_mean', 'invth_simple_mean')]
data[,FEP_PM2_5_CO_interacted:= FEP_PM2_5_mean_city*FEP_CO_mean_city]
data[,FEP_O3_CO_interacted:= FEP_O3_mean_city*FEP_CO_mean_city]
data[,FEP_O3_SO2_interacted:= FEP_O3_mean_city*FEP_SO2_mean_city]
data[,FEP_PM2_5_SO2_interacted:= FEP_PM2_5_mean_city*FEP_SO2_mean_city]

data5<-data[all5in==T]

# Weather -- R?vision ici
m<-c("FEP_neige_sol","FEP_precipitations","FEP_precipitations_2",
     "FEP_temperature","FEP_temperature_2","FEP_ff","FEP_ff_2","FEP_humidite","FEP_humidite_2","FEP_ensoleil_glo")

# Pollutants 
pollutants <- c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_O3_mean_city","FEP_NO2_mean_city","FEP_SO2_mean_city")
pollutants_max <- paste0(pollutants, "_max24h")
pollutants_min <- paste0(pollutants, "_min24h")
pollutants6 <- c('FEP_PM10_mean_city',pollutants)
pollutant_interaction <- c('FEP_PM2_5_CO_interacted', 'FEP_O3_CO_interacted', "FEP_O3_SO2_interacted", "FEP_PM2_5_SO2_interacted")

# Outcomes
outcomes<-c("FEP_rate_adm_C_10_Full_all","FEP_rate_adm_C_09_Full_all","FEP_rate_adm_C_11_Full_all","FEP_n_J_in_CertD","FEP_n_I_in_CertD","FEP_n_K_in_CertD")
totaloutcomes <- c("FEP_rate_emergencies_tot","FEP_rate_dec_all")
  
####################################################################################################################################
####################### ====================== Lasso selection ============================= #######################################
####################################################################################################################################


# Instruments.
# 3 output variables from LMDZ among the variables not indexed by pressure level (3)
# 4 output variables indexed by pressure levels, + wind norm, at 10-79 layers (69*5)
# PBLH, INVPBLH, INVTH, INVTH_strength: on average + on average at 6 moments of the day.
# IPBLH*moments-of-the-day*city

### Choice of the instrument set.

instruments <- c(names(data)[grepl('FEP_s_pblh4hours', names(data))],
                 names(data)[grepl('FEP_s_pblh_mean$', names(data))],
                 
                 names(data)[grepl('FEP_invth_simple_[0-5]$',names(data))],   
                 names(data)[grepl('FEP_invth_strength_[0-5]$',names(data))],  
                 names(data)[grepl('FEP_invth_simple_mean$',names(data))],
                 names(data)[grepl('FEP_invth_strength_mean$',names(data))],
                 
                 names(data)[grepl('FEP_inv_hcl_[0-5]$',names(data))],    
                 names(data)[grepl('FEP_intcity_inv_hcl_[0-5]_[A-Z]',names(data))],
                 names(data)[grepl('FEP_inv_hcl_mean$',names(data))],    
                 
                 names(data)[grepl('FEP_vitu_[2-9][0-9]_mean$', names(data))],
                 names(data)[grepl('FEP_vitv_[2-9][0-9]_mean$', names(data))],
                 names(data)[grepl('FEP_norm_v_[2-9][0-9]_mean$', names(data))],
                 names(data)[grepl('FEP_zfull_[2-9][0-9]_mean$', names(data))]) 

# Robustness (1): hourly measures for IPBLH, PBLH, INVTH.
# Robustness (2): interacted terms for invth as well
# Robustness (3): no interacted terms.
instruments1 <- c(instruments,
                  names(data)[grepl('FEP_s_pblh_[0-9]+$', names(data))], 
                  names(data)[grepl('FEP_ipblh_',names(data))],
                  names(data)[grepl('FEP_hourly_invth',names(data))])

instruments2 <- c(instruments,names(data)[grepl('FEP_intcity_invth',names(data))])

instruments3 <- setdiff(instruments,
                        names(data)[grepl('FEP_intcity_inv_hcl_[0-5]_[A-Z]',names(data))])

instruments4 <- c(instruments, 
                  names(data)[grepl('FEP_ovap_[2-9][0-9]_mean$', names(data))])

instruments5 <- setdiff(instruments, 
                  names(data)[grepl('FEP_zfull_[2-9][0-9]_mean$', names(data))])

instruments6 <- c(names(data)[grepl('FEP_s_pblh4hours', names(data))],
                 names(data)[grepl('FEP_s_pblh_mean$', names(data))],
                 
                 names(data)[grepl('FEP_invth_simple_[0-5]$',names(data))],   
                 names(data)[grepl('FEP_invth_strength_[0-5]$',names(data))],  
                 names(data)[grepl('FEP_invth_simple_mean$',names(data))],
                 names(data)[grepl('FEP_invth_strength_mean$',names(data))],
                 
                 names(data)[grepl('FEP_inv_hcl_[0-5]$',names(data))],    
                 names(data)[grepl('FEP_intcity_inv_hcl_[0-5]_[A-Z]',names(data))],
                 names(data)[grepl('FEP_inv_hcl_mean$',names(data))],    
                 
                 names(data)[grepl('FEP_vitu_[1-9][0-9]_mean$', names(data))],
                 names(data)[grepl('FEP_vitv_[1-9][0-9]_mean$', names(data))],
                 names(data)[grepl('FEP_norm_v_[1-9][0-9]_mean$', names(data))],
                 names(data)[grepl('FEP_zfull_[1-9][0-9]_mean$', names(data))]) 

instruments7 <- c(names(data)[grepl('FEP_s_pblh4hours', names(data))],
                  names(data)[grepl('FEP_s_pblh_mean$', names(data))],
                  
                  names(data)[grepl('FEP_invth_simple_[0-5]$',names(data))],   
                  names(data)[grepl('FEP_invth_strength_[0-5]$',names(data))],  
                  names(data)[grepl('FEP_invth_simple_mean$',names(data))],
                  names(data)[grepl('FEP_invth_strength_mean$',names(data))],
                  
                  names(data)[grepl('FEP_inv_hcl_[0-5]$',names(data))],    
                  names(data)[grepl('FEP_intcity_inv_hcl_[0-5]_[A-Z]',names(data))],
                  names(data)[grepl('FEP_inv_hcl_mean$',names(data))],    
                  
                  names(data)[grepl('FEP_vitu_[3-9][0-9]_mean$', names(data))],
                  names(data)[grepl('FEP_vitv_[3-9][0-9]_mean$', names(data))],
                  names(data)[grepl('FEP_norm_v_[3-9][0-9]_mean$', names(data))],
                  names(data)[grepl('FEP_zfull_[3-9][0-9]_mean$', names(data))]) 




instruments_set_list <- list(instruments, instruments1, instruments2, instruments3,instruments4, instruments5, instruments6, instruments7)

res <-  lapply(instruments_set_list, FUN = function(inst) lapply(pollutants, FUN = function(x)
  selected_instruments(data5, m, inst, x, clustervar = 'month_year_ville')))
res_max <- lapply(pollutants_max, FUN = function(x)
  selected_instruments(data5, m, instruments_set_list[[1]], x, clustervar = 'month_year_ville',baseline_only = T))
res_min <- lapply(pollutants_min, FUN = function(x)
  selected_instruments(data5, m, instruments_set_list[[1]], x, clustervar = 'month_year_ville',baseline_only = T))
res_sampleA <- lapply(pollutants, FUN = function(x)
  selected_instruments(data, m, instruments_set_list[[1]], x, clustervar = 'month_year_ville',baseline_only = T))
res_interaction <- lapply(pollutant_interaction, FUN = function(x)
  selected_instruments(data, m, instruments_set_list[[1]], x, clustervar = 'month_year_ville',baseline_only = T))

save(instruments_set_list, res, res_max, res_min, res_sampleA, res_interaction, file = "data/selection_step_rev.Rdata")

### Robustness Lasso Selection to Weather Controls.

m1 <-  c(setdiff(m, c('FEP_temperature','FEP_temperature_2','FEP_precipitations','FEP_precipitations_2')), 
         'FEP_temperature_min', 'FEP_temperature_min_2', 'FEP_temperature_max', 'FEP_temperature_max_2', 
         'FEP_precipitations_min', 'FEP_precipitations_min_2', 'FEP_precipitations_max', 'FEP_precipitations_max_2')

m2 <- c(setdiff(m, c('FEP_temperature','FEP_temperature_2')), names(data)[grepl('^FEP_temperature_bin_[0-9]+$',names(data))])

m3 <- c(setdiff(m, c('FEP_temperature','FEP_temperature_2','FEP_ff','FEP_ff_2','FEP_humidite','FEP_humidite_2','FEP_precipitations','FEP_precipitations_2')),   
        names(data)[grepl('^FEP_temperature_bin_[0-9]+',names(data))],
        names(data)[grepl('(^FEP_ff_bin_set3)|(^FEP_humidite_bin_set3)|(^FEP_precipitations_bin_)',names(data))])

m4 <- c(setdiff(m, c('FEP_temperature','FEP_temperature_2','FEP_ff','FEP_ff_2','FEP_humidite','FEP_humidite_2','FEP_precipitations','FEP_precipitations_2')),   
        names(data)[grepl('^FEP_temperature_bin_set3',names(data))],
        names(data)[grepl('(^FEP_ff_bin_set3)|(^FEP_humidite_bin_set3)|(^FEP_precipitations_bin_)',names(data))])

m5 <- c(setdiff(m, c('FEP_temperature','FEP_temperature_2')), names(data)[grepl('^FEP_temperature_bin_breaks',names(data))])

m6 <- c(m, names(data)[grepl('^FEP_dd_binwd',names(data))])

m7 <- c(setdiff(m, c('FEP_temperature','FEP_temperature_2','FEP_ff','FEP_ff_2','FEP_humidite','FEP_humidite_2','FEP_precipitations','FEP_precipitations_2')),   
        names(data)[grepl('^FEP_temperature_bin_decile_[0-9]+$',names(data))],
        names(data)[grepl('(^FEP_ff_bin_set3)|(^FEP_humidite_bin_set3)|(^FEP_precipitations_bin_)',names(data))],
        names(data)[grepl('^FEP_dd_binwd',names(data))])

m8 <- c(setdiff(m, c('FEP_temperature','FEP_temperature_2')), names(data)[grepl('^FEP_temperature_bin_decile_',names(data))])

weather_list_rob <- list(m1,m2,m3,m4, m5, m6, m7, m8)

res_robW <- lapply(weather_list_rob, FUN = function(mw)  lapply(pollutants, FUN = function(x) selected_instruments(data5, mw, instruments, x, clustervar = 'month_year_ville', baseline_only = TRUE)))
res_robW_largersample  <- lapply(weather_list_rob, FUN = function(mw)  lapply(pollutants, FUN = function(x) selected_instruments(data, mw, instruments, x, clustervar = 'month_year_ville', baseline_only = TRUE)))

save(weather_list_rob, res_robW, res_robW_largersample, file = "data/selection_step_variant_W_rev.Rdata")




### Selected Sets.

load("data/selection_step_rev.Rdata")
baseline <- lapply(1:length(res), FUN = function(i) (lapply(res[[i]], FUN=function(x) setdiff(x$selected_robustunpC, m))))
baselineNC <- lapply(1:length(res), FUN = function(i) (lapply(res[[i]], FUN=function(x) setdiff(x$selected_robust, m))))
baseline_max <- unique(unlist(res_max))
baseline_min <- unique(unlist(res_min))
baseline_sampleA <- lapply(1:length(res_sampleA), FUN = function(i) unlist(res_sampleA[[i]]))
baseline_interactions <- lapply(1:length(res_interaction), FUN = function(i) unlist(res_interaction[[i]]))

lapply(baseline[[1]], length)
length(baseline_max)
lapply(baseline_sampleA, length)
lapply(baseline_interactions, length)

selected_set <- unique(unlist(baseline[[1]]))
setdiff(baseline_max,selected_set)
setdiff(selected_set,baseline_max)

#selected_set<-lapply(pollutants, FUN = function(x) selected_instruments(data5, m, instruments, x, clustervar = 'month_year_ville', baseline_only = TRUE))

length(selected_set)
selected_set_ordered <- c("FEP_invth_simple_0","FEP_invth_simple_2","FEP_invth_simple_5","FEP_invth_strength_0", "FEP_invth_strength_1",
                          "FEP_s_pblh4hours_0","FEP_s_pblh4hours_2","FEP_s_pblh4hours_3",
                          "FEP_inv_hcl_mean","FEP_inv_hcl_0", "FEP_inv_hcl_2","FEP_inv_hcl_1","FEP_inv_hcl_4",
                          "FEP_intcity_inv_hcl_0_LILLE","FEP_intcity_inv_hcl_0_NICE", "FEP_intcity_inv_hcl_0_PARIS",
                          "FEP_intcity_inv_hcl_1_NICE","FEP_intcity_inv_hcl_1_LILLE",
                          "FEP_intcity_inv_hcl_2_LYON","FEP_intcity_inv_hcl_2_NICE", "FEP_intcity_inv_hcl_2_NANTES", 
                          "FEP_intcity_inv_hcl_3_NANTES",
                          "FEP_intcity_inv_hcl_4_NANTES","FEP_intcity_inv_hcl_5_MARSEILLE","FEP_intcity_inv_hcl_5_STRASBOURG",
                           "FEP_vitu_20_mean","FEP_vitu_40_mean","FEP_vitv_32_mean",
                           "FEP_norm_v_38_mean","FEP_norm_v_39_mean","FEP_norm_v_45_mean","FEP_norm_v_52_mean",
                           "FEP_zfull_25_mean","FEP_zfull_46_mean","FEP_zfull_78_mean")

# What about prediction?
pred_qualit <- function(pol, select_inst_list, polname, quant_seq = seq(0.04,0.96,0.02), insample = FALSE, inst_name = c("PM2.5","CO","O3","NO2","SO2")){
  lms <- lapply(select_inst_list, FUN = function(x) lm(as.formula(paste(pol," ~ ",paste0(x, collapse = "+")) ), data = data5))
  haty <- lapply(lms, FUN = function(x) predict(x, data[all5in == insample]))
  y <- data[all5in== insample,get(pol)]
  print(lapply(haty, FUN =function(x) summary(lm(x ~ y))))
  df <- data.table(y = rep(y, length(select_inst_list)), haty = unlist(haty), inst = rep(inst_name, each = length(y)))
  df[, quantY:= cut(get('y'), quantile(y, quant_seq, na.rm = T)), by = inst]
  print(df[,.(sum(!is.na(get('y')))),by = inst])
  df <- df[,.(haty = median(haty, na.rm=T), y = median(y, na.rm=T)), by = .(quantY,inst)]
  df[, isMain:= get('inst') ==  polname]
  ggplot(df[!is.na(quantY)],aes(y,haty, col = as.factor(inst), group = as.factor(inst), alpha = isMain)) + 
    geom_point() + geom_line()  + ggtitle(label = ifelse(insample, "Prediction In Sample", "Prediction Out-of-Sample")) + 
    xlab(paste("Observed ", polname)) + ylab(paste("Predicted ",polname)) + theme(legend.title = element_blank())  + theme_bw() + scale_alpha_manual(values = c(0.3,1))
}
pred_qualit('FEP_PM2_5_mean_city', baseline[[1]], "PM2.5", quant_seq = seq(0.05,0.95,0.05), insample = TRUE)
pred_qualit('FEP_PM2_5_mean_city', baseline[[1]], "PM2.5", quant_seq = seq(0.05,0.95,0.05), insample = FALSE)
pred_qualit('FEP_CO_mean_city', baseline[[1]], "CO", insample = TRUE)
pred_qualit('FEP_CO_mean_city', baseline[[1]], "CO", insample = FALSE)
pred_qualit('FEP_O3_mean_city', baseline[[1]], "O3", insample = TRUE)
pred_qualit('FEP_O3_mean_city', baseline[[1]], "O3", insample = FALSE)
pred_qualit('FEP_NO2_mean_city', baseline[[1]], "NO2", insample = TRUE)
pred_qualit('FEP_NO2_mean_city', baseline[[1]], "NO2", insample = FALSE)
pred_qualit('FEP_SO2_mean_city', baseline[[1]], "SO2", insample = TRUE)
pred_qualit('FEP_SO2_mean_city', baseline[[1]], "SO2", insample = FALSE)



pred_qualit_v0 <- function(pol, select_inst, polname){
  pm <- lm(as.formula(paste(pol," ~ ",paste0(select_inst, collapse = "+"))), data = data5)
  haty <- predict(pm, data[all5in == F])
  y <- data[all5in== F,get(pol)]
  print(summary(lm(y~haty)))
  print(summary(lm(haty~y)))
  plot(y,haty, xlab = paste(polname," Observed Out Of Sample"), ylab = paste(polname," Predicted Out Of Sample")) + abline(0,1, col = "blue") 
}
pred_qualit_v0('FEP_CO_mean_city', baseline[[1]][[1]], "CO")
pred_qualit_v0('FEP_CO_mean_city', baseline[[1]][[2]], "CO")

###################################################################################################################
################# A Stata Sample For Wild Bootstrap, First Stage etc, F-Stat etc... ###############################
###################################################################################################################

pm10_selection <-  selected_instruments(data5, m, instruments, "FEP_PM10_mean_city", clustervar = 'month_year_ville', baseline_only = TRUE)
selected_6 <- union(selected_set, unlist(pm10_selection))
x<-c(outcomes,pollutants6,m,selected_6)
y <- gsub('FEP_','', x)
sd <- data[, lapply(.SD,FUN =function(x) sd(x, na.rm=T))  , .SDcols= y]
tot <- c("ville","date","month_year_ville","month_year","day","all5in","pollution_index_pca_imputed_F1","invth_simple_mean","FEP_invth_simple_mean",x,y) 
tot <- unique(tot)
stata <- data[,..tot]
lapply(1:length(x), FUN = function(i) stata[,paste0(gsub("_intcity","",x[i]),"_sc") := get(x[i])/sd[[i]]]) 
write_dta(stata,path="data/base_IV_for_stata.dta")


####################################################################
################# On the First Stage ###############################
####################################################################

# Full First Stage (Compare With Stata) - scaled first stage Stata
firststage <- lapply(pollutants,  FUN = function(i) reg_FE_clustered(data5, x2 = m, y = i, x = selected_set_ordered, FE = '0'))
stargazer(firststage,omit.stat=c("ser","F","adj.rsq"), digits = 2,   omit = c(m,'Constant'),
          dep.var.labels.include = F, column.labels = c("PM2.5","CO","O3","NO2","SO2"), column.separate = c(1,1,1,1,1,1),
          add.lines =list(c("# Selected", unlist(lapply(baseline[[1]], length)))))

# Intersection between selected sets (why do we exclude PM10)
intersection <- lapply(pollutants6, FUN = function(x) selected_instruments(data5, m, instruments_set_list[[1]], x, clustervar = 'month_year_ville', baseline_only = TRUE))
save(intersection, file = "data/selection_6_pollutants_baseline.Rdata")
load("data/selection_6_pollutants_baseline.Rdata")
cross <- sapply(intersection, FUN = function(y) sapply(intersection, FUN = function(x) length(setdiff(intersect(x$selected_robustunpC, y$selected_robustunpC),m))/length(setdiff(x$selected_robustunpC,m))))
cross
min(cross[cross!=0])
max(cross[cross!=1])
mean(cross[cross!=1 & cross!=0])
crossnopm10 <- cross[c(1,3:6),c(1,3:6)]
mean(crossnopm10[crossnopm10!=1 & crossnopm10!=0])


FSV <- lapply(pollutants6, 
              FUN = function(x) scale_penalty(data5, m, instruments_set_list[[1]], x, clustervar = 'month_year_ville', grid = seq(1,3,0.2)))
save(FSV, file = 'data/selection_with_higher_penalty_rev.Rdata')
load('data/selection_with_higher_penalty_rev.Rdata')
lapply(FSV[[1]]$selected_around_robustC,FUN = function(x) length(setdiff(x,m)))
lapply(FSV[[2]]$selected_around_robustC,FUN = function(x) length(setdiff(x,m)))
lapply(FSV[[3]]$selected_around_robustC,FUN = function(x) length(setdiff(x,m)))
lapply(FSV[[4]]$selected_around_robustC,FUN = function(x) length(setdiff(x,m)))
lapply(FSV[[5]]$selected_around_robustC,FUN = function(x) length(setdiff(x,m)))
lapply(FSV[[6]]$selected_around_robustC,FUN = function(x) length(setdiff(x,m)))

setdiff(FSV[[1]]$selected_around_robustC[[6]], m) #5 
setdiff(FSV[[2]]$selected_around_robustC[[11]], m) #5
setdiff(FSV[[4]]$selected_around_robustC[[4]], m) #5
setdiff(FSV[[5]]$selected_around_robustC[[7]], m) #4
setdiff(FSV[[6]]$selected_around_robustC[[1]], m) #3

pm10 <- setdiff(FSV[[1]]$selected_around_robustC[[6]], m)
pm25 <- setdiff(FSV[[2]]$selected_around_robustC[[11]], m)
co <-  setdiff(FSV[[3]]$selected_around_robustC[[1]], m)
o3 <- setdiff(FSV[[4]]$selected_around_robustC[[4]], m)
no2 <- setdiff(FSV[[5]]$selected_around_robustC[[7]], m)
so2 <- setdiff(FSV[[6]]$selected_around_robustC[[1]], m)

setdiff(pm25, baseline[[1]][[1]])
setdiff(co, baseline[[1]][[2]])
setdiff(o3, baseline[[1]][[3]])
setdiff(no2, baseline[[1]][[4]])
setdiff(so2, baseline[[1]][[5]])

small <- list(pm10, pm25, co, o3, no2, so2)


## Table Lasso
firststage <- lapply(2:6,  FUN = function(i) reg_FE_clustered(data5, x2 = m, y = pollutants6[i], x = small[[i]], FE = '0'))
stargazer(firststage,omit.stat=c("ser","F","adj.rsq"), digits = 2,   omit = c(m,'Constant'),  type = "text",
          dep.var.labels.include = F, column.labels = c("PM2.5","CO","O3","NO2","SO2"), column.separate = c(1,1,1,1,1,1),
          add.lines =list(c("# Selected (Optimal Constraint)", unlist(lapply(baseline[[1]], length)))))
# With normalization.
y <- c(pollutants6, unique(unlist(small)))
sd <- data[, lapply(.SD,FUN =function(x) sd(x, na.rm=T))  , .SDcols= gsub("FEP_","",y)]
tot <- c(y, m, "month_year_ville")
scaled <- data5[,..tot]
lapply(1:length(y), FUN = function(i) scaled[,paste0(y[i],"") := get(y[i])/sd[[i]]])
firststage <- lapply(2:6,  FUN = function(i) reg_FE_clustered(scaled, x2 = m, y = pollutants6[i], x = small[[i]], FE = '0'))
stargazer(firststage,omit.stat=c("ser","F","adj.rsq"), digits = 2,   omit = c(m,'Constant'), 
          dep.var.labels.include = F, column.labels = c("PM2.5","CO","O3","NO2","SO2"), column.separate = c(1,1,1,1,1,1),
          add.lines =list(c("# Selected (Optimal Constraint)", unlist(lapply(baseline[[1]], length)))))

## Table Lasso

firststage <- lapply(2:6,  FUN = function(i) reg_FE_clustered(data5, x2 = m, y = pollutants6[i], x = baseline[[1]][[i-1]], FE = '0'))
stargazer(firststage,omit.stat=c("ser","F","adj.rsq"), digits = 2,   omit = c(m,'Constant'),  type = "text",
          dep.var.labels.include = F, column.labels = c("PM2.5","CO","O3","NO2","SO2"), column.separate = c(1,1,1,1,1,1),
          add.lines =list(c("# Selected (Optimal Constraint)", unlist(lapply(baseline[[1]], length)))))
# With normalization.
y <- c(pollutants6, unique(unlist(baseline[[1]])))
sd <- data[, lapply(.SD,FUN =function(x) sd(x, na.rm=T))  , .SDcols= gsub("FEP_","",y)]
tot <- c(y, m, "month_year_ville")
scaled <- data5[,..tot]
lapply(1:length(y), FUN = function(i) scaled[,paste0(y[i],"") := get(y[i])/sd[[i]]])
firststage <- lapply(2:6,  FUN = function(i) reg_FE_clustered(scaled, x2 = m, y = pollutants6[i], x = baseline[[1]][[i-1]], FE = '0'))
stargazer(firststage,omit.stat=c("ser","F","adj.rsq"), digits = 2,   omit = c(m,'Constant'), 
          dep.var.labels.include = F, column.labels = c("PM2.5","CO","O3","NO2","SO2"), column.separate = c(1,1,1,1,1,1),
          add.lines =list(c("# Selected (Optimal Constraint)", unlist(lapply(baseline[[1]], length)))))



## Ilustrate in an appendix
pm25_selection_6P_sample <-  selected_instruments(data5[!is.na(FEP_PM10_mean_city)], m, instruments, "FEP_PM2_5_mean_city", clustervar = 'month_year_ville', baseline_only = TRUE)
comp <- list(unlist(pm10_selection), unlist(pm25_selection_6P_sample))
firststage <- lapply(1:2,  FUN = function(i) reg_FE_clustered(data5[!is.na(FEP_PM10_mean_city)], x2 = m, y = pollutants6[i], x = comp[[i]], FE = '0'))
stargazer(firststage,omit.stat=c("ser","F","adj.rsq"), digits = 2,   omit = c(m,'Constant'), type = "text",
          dep.var.labels.include = F, column.labels = c("PM10","PM2.5","CO","O3","NO2","SO2"), column.separate = c(1,1,1,1,1,1),
          add.lines =list(c("# Selected (Optimal Constraint)", unlist(lapply(baseline[[1]], length)))))

## Largest Sample for each pollutant

pm10 <- scale_penalty(data, m, instruments_set_list[[1]], "FEP_PM10_mean_city", clustervar = 'month_year_ville', grid = seq(3,5,0.2))
pm25 <- scale_penalty(data, m, instruments_set_list[[1]], "FEP_PM2_5_mean_city", clustervar = 'month_year_ville', grid = seq(4,6,0.2))
co <- scale_penalty(data, m, instruments_set_list[[1]], "FEP_CO_mean_city", clustervar = 'month_year_ville', grid = seq(2,5,0.2))
no2 <- scale_penalty(data, m, instruments_set_list[[1]], "FEP_NO2_mean_city", clustervar = 'month_year_ville', grid = seq(3,5,0.2))
so2 <- scale_penalty(data, m, instruments_set_list[[1]], "FEP_SO2_mean_city", clustervar = 'month_year_ville', grid = seq(1,2,0.2))
o3 <- scale_penalty(data, m, instruments_set_list[[1]], "FEP_O3_mean_city", clustervar = 'month_year_ville', grid = seq(2,3,0.2))
save(pm10,pm25,co,no2,so2,o3, file = 'data/selection_with_higher_penalty_large_sample_rev.Rdata')
load('data/selection_with_higher_penalty_large_sample_rev.Rdata')

lapply(pm10[[1]],length)
lapply(pm25[[1]],length)
lapply(co[[1]],length)
lapply(no2[[1]],length)
lapply(so2[[1]],length)
lapply(o3[[1]],length)

setdiff(pm25[[1]][[3]], baseline[[1]][[1]])
pm25[[1]][[3]]
setdiff(co[[1]][[4]], baseline[[1]][[2]])
co[[1]][[4]]
setdiff(o3[[1]][[6]], baseline[[1]][[3]])
o3[[1]][[6]]
setdiff(no2[[1]][[3]], baseline[[1]][[4]])
no2[[1]][[3]]
setdiff(so2[[1]][[4]], baseline[[1]][[5]])
so2[[1]][[4]]
small <- list(pm10[[1]][[1]], pm25[[1]][[3]], co[[1]][[4]], o3[[1]][[6]], no2[[1]][[3]], so2[[1]][[4]])

firststage <- lapply(1:6,  FUN = function(i) reg_FE_clustered(data, x2 = m, y = pollutants6[i], x = small[[i]], FE = '0'))
firststage2 <- lapply(1:6,  FUN = function(i) reg_FE_clustered(data5, x2 = m, y = pollutants6[i], x = small[[i]], FE = '0'))
stargazer(firststage,firststage2,omit.stat=c("ser","F","adj.rsq"), digits = 2,   omit = c(m,'Constant'), type = "text",
          dep.var.labels.include = F, column.labels = c("PM10","PM2.5","CO","O3","NO2","SO2"), column.separate = c(1,1,1,1,1,1))


################################################################################################################################################### 
################### =============== IV LASSO - Health Outcomes, Five Pollutant models  ===================== ######################################
###################################################################################################################################################
# Panel A, Table 6
a <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, c("-1",m), pollutants, selected_set))

stargazer(a, omit.stat=c("ser","adj.rsq","rsq"), 
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F, digits = 4, type = "text")

# 5% coldest and 5% hottest
a <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5[temperature <23.3 & temperature >1.5], vd, c("-1",m), pollutants, selected_set))
stargazer(a, omit.stat=c("ser","adj.rsq","rsq"), type = "text",
          dep.var.labels.include = F, digits = 4)
a <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5[humidite <95.3 & humidite >48.7], vd, c("-1",m), pollutants, selected_set))
stargazer(a, omit.stat=c("ser","adj.rsq","rsq"), type = "text",
          dep.var.labels.include = F, digits = 4)


# To compare with wild bootstrap p values
stargazer(a[[1]],a[[2]],a[[4]],a[[5]], omit.stat=c("ser","adj.rsq","rsq"),   covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"), 
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol","tsv"),
          dep.var.labels.include = F, digits = 4, report = "vc*p")

tot <- lapply(totaloutcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants,selected_set))
stargazer(tot, omit.stat=c("ser","adj.rsq","rsq"),   covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"), 
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol","tsv"),
          dep.var.labels.include = F, digits = 4)

## Alternative sets of instruments
load("data/selection_step_rev.Rdata")
alt_inst <- lapply(outcomes, FUN =function(vd) lapply(c(1,3:8), FUN = function(i) 
  twoSLS_on_selected(data5, vd, m, pollutants, 
    setdiff(unique(c(res[[i]][[1]]$selected_robustunpC,
                     res[[i]][[2]]$selected_robustunpC,
                     res[[i]][[3]]$selected_robustunpC,
                     res[[i]][[4]]$selected_robustunpC,
                     res[[i]][[5]]$selected_robustunpC)),m) )))

stargazer(alt_inst[[1]], omit.stat=c("ser","adj.rsq","rsq"), column.labels = c("Ref", "+Interactions InvTH", "- Interactions","+Humidity","-Altitude Pres"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report = "vc*p",
          dep.var.labels.include = F, digits = 4, type = "text")
stargazer(alt_inst[[2]], omit.stat=c("ser","adj.rsq","rsq"), column.labels = c("Ref", "+Interactions InvTH", "- Interactions","+Humidity","-Altitude Pres"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report = "vc*p",
          dep.var.labels.include = F, digits = 4, type = "text")
stargazer(alt_inst[[3]], omit.stat=c("ser","adj.rsq","rsq"), column.labels = c("Ref", "+Interactions InvTH", "- Interactions","+Humidity","-Altitude Pres"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report = "vc*p",
          dep.var.labels.include = F, digits = 4, type = "text")
stargazer(alt_inst[[4]], omit.stat=c("ser","adj.rsq","rsq"), column.labels = c("Ref", "+Interactions InvTH", "- Interactions","+Humidity","-Altitude Pres"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report = "vc*p",
          dep.var.labels.include = F, digits = 4, type = "text")
stargazer(alt_inst[[5]], omit.stat=c("ser","adj.rsq","rsq"), column.labels = c("Ref", "+Interactions InvTH", "- Interactions","+Humidity","-Altitude Pres"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report = "vc*p",
          dep.var.labels.include = F, digits = 4, type = "text")
stargazer(alt_inst[[6]], omit.stat=c("ser","adj.rsq","rsq"), column.labels = c("Ref", "+Interactions InvTH", "- Interactions","+Humidity","-Altitude Pres"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report = "vc*p",
          dep.var.labels.include = F, digits = 4, type = "text")

## For Online Appendix.

n_inst<-unlist(lapply(1:6, FUN = function(i) length(setdiff(unique(c(res[[i]][[1]]$selected_robustunpC,
                                                                     res[[i]][[2]]$selected_robustunpC,
                                                                     res[[i]][[3]]$selected_robustunpC,
                                                                     res[[i]][[4]]$selected_robustunpC,
                                                                     res[[i]][[5]]$selected_robustunpC)),m))))

n_inst_init<-unlist(lapply(instruments_set_list,length))


stargazer(alt_inst[[1]],alt_inst[[2]], omit.stat=c("ser","adj.rsq","rsq"), column.labels = rep(c("Baseline", "(+ City x TI)",
                                                                                             "(- City x IPBLH) ","(+ Humidity)","(- Altitude Pres)"),2),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          covariate.labels = c("PM2.5","CO","O3","NO2","SO2"), type = "text",
          dep.var.labels.include = F, digits = 4, 
          add.lines = list(c("Instrument Set Size",n_inst_init[c(1,3:6)],n_inst_init[c(1,3:6)]),
                           c("Selected Instrument Set Size",n_inst[c(1,3:6)],n_inst[c(1,3:6)])))
stargazer(alt_inst[[4]],alt_inst[[5]], omit.stat=c("ser","adj.rsq","rsq"), column.labels = rep(c("Baseline", "(+ City x TI)",
                                                                                             "(- City x IPBLH) ","(+ Humidity)","(- Altitude Pres)"),2),
          add.lines = list(c("Instrument Set Size",n_inst_init[c(1,3:6)],n_inst_init[c(1,3:6)]),
                           c("Selected Instrument Set Size",n_inst[c(1,3:6)],n_inst[c(1,3:6)])),

                    covariate.labels = c("PM2.5","CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F, digits = 4)




# Randomized Inference For Standard Errors As Suggested by Referee 1
r <- randomization_estimates_ville(100, data, a)
# Estimates.
comp_sd <- r[,.(pvalless5 = round(mean(pval <0.05),2), mean_est = round(mean(est),4), est_ref = round(first(est_ref),4), sd = round(sd(est),5), sd_ref = round(first(se_ref),5)) , by = .(outcome, pol)]
comp_sd[,.(pol, outcome,mean_est,est_ref, sd, sd_ref)]

# Referee Suggestion: add atmospheric pressure
# From LMDZ: psol Surface Pressure Pa
ar <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, c(m,"FEP_pstat", "I(FEP_pstat^2)"), pollutants, selected_set))

stargazer(a[[1]],ar[[1]],a[[2]],ar[[2]],a[[4]],ar[[4]],a[[5]],ar[[5]],omit.stat=c("ser","adj.rsq","rsq"),  type = "text", 
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol","pstat"),
          dep.var.labels.include = F, digits = 4)

# Referee Suggestion: Report in Table 6, Rate of emergencies admissions. Comment on substitution Effects 
s0<-(twoSLS_on_selected(data5, 'FEP_rate_emergencies_tot', m, pollutants, selected_set))
s1<-(twoSLS_on_selected(data5, 'FEP_rate_emergencies_9_10', m, pollutants, selected_set))
s2<-(twoSLS_on_selected(data5, 'FEP_rate_emergencies_not_9_10', m, pollutants, selected_set))
# 
stargazer(s0,s1,s2,omit.stat=c("ser","adj.rsq","rsq"),  type = "text",
          dep.var.labels.include = F,digits=4, omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"))

# When only an interaction term is selected for a given city, complete with the same interaction term with all cities as a robustness
int <- paste0('FEP_intcity_inv_hcl', unique(substr(selected_set[grepl("intcity",selected_set)],20,22)))
selected_set_extended <- unique(c(selected_set, unlist(lapply(int, FUN = function(i) names(data)[grepl(i,names(data))]))))
ar <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants, selected_set_extended))

stargazer(a[[1]],ar[[1]],a[[2]],ar[[2]],a[[4]],ar[[4]],a[[5]],ar[[5]], omit.stat=c("ser","adj.rsq","rsq"), 
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F, digits = 4, type = "text")


#######################
####### MHT ###########
#######################

MHT1 <- MHT_corrections(alpha = 0.05, grouping = "outcome", outc =  c("FEP_rate_adm_C_10_Full_all","FEP_rate_adm_C_09_Full_all","FEP_rate_dec_all","FEP_rate_adm_C_11_Full_all","FEP_n_I","FEP_n_J","FEP_n_I_in_CertD","FEP_n_J_in_CertD"))
MHT1[ pval <0.05 ]
MHT1 <- MHT_corrections(alpha = 0.1, grouping = "outcome", outc =  c("FEP_rate_adm_C_10_Full_all","FEP_rate_adm_C_09_Full_all","FEP_rate_dec_all","FEP_rate_adm_C_11_Full_all","FEP_n_I","FEP_n_J","FEP_n_I_in_CertD","FEP_n_J_in_CertD"))
MHT1[ pval < 0.1 ]
MHT2 <- MHT_corrections(alpha = 0.1, grouping = "", outc =  c("FEP_rate_adm_C_10_Full_all","FEP_rate_adm_C_09_Full_all","FEP_rate_adm_C_11_Full_all"))
MHT2[FWEReject == TRUE | FDRReject == TRUE]
MHT3 <- MHT_corrections(alpha = 0.1, grouping = "", outc =  c("FEP_n_I_in_CertD","FEP_n_J_in_CertD","FEP_n_K_in_CertD"))
MHT3[FWEReject == TRUE | FDRReject == TRUE]
# Possibilit? de BKY :
MHT3 <- MHT_corrections(alpha = 0.1/(0.1+1), grouping = "", outc =  c("FEP_rate_adm_C_10_Full_all","FEP_rate_adm_C_09_Full_all","FEP_n_I_in_CertD","FEP_n_J_in_CertD"))
MHT3[FWEReject == TRUE | FDRReject == TRUE]
r<-MHT3[FDRReject==T][,.N]
MHT3 <- MHT_corrections(alpha = 0.1*(dim(MHT3)[1])/((dim(MHT3)[1]-r)*1.1), grouping = "", outc =  c("FEP_rate_adm_C_10_Full_all","FEP_rate_adm_C_09_Full_all","FEP_n_I_in_CertD","FEP_n_J_in_CertD"))
MHT3[FWEReject == TRUE | FDRReject == TRUE]
MHT4 <- MHT_corrections(alpha = 0.05, grouping = "", outc =  c("FEP_rate_adm_C_10_Full_all","FEP_rate_adm_C_09_Full_all","FEP_n_I_in_CertD","FEP_n_J_in_CertD"))
MHT4[FWEReject == TRUE | FDRReject == TRUE]


################################################################################################################################################### 
################### ============ IV LASSO -  Health Outcomes, Single Pollutant models  ===================== ######################################
###################################################################################################################################################

tab <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(1:5, FUN = function(i) twoSLS_on_selected(data5, vd, m, pollutants[i],baseline[[1]][[i]])))
tab2 <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(pollutants, FUN = function(x) twoSLS_on_selected(data5, vd, m, x,selected_set)))

stargazer(tab, omit.stat=c("ser","adj.rsq","rsq"),   covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4, type = "text", title = "Only pollutant-specific instrument for each pollutant")
stargazer(tab2, omit.stat=c("ser","adj.rsq","rsq"),   covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4, type = "text", title = "All selected instruments used for each pollutant")

# Report in Table 7.
tab.s.c <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data5, vd, m, pollutants[i],baseline[[1]][[i]]))$coefficients[,1]; u[grepl('mean_city',names(u))]} ))
tab.s.se <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(1:5, FUN = function(i){ u<- summary(twoSLS_on_selected(data5, vd, m, pollutants[i],baseline[[1]][[i]]))$coefficients[,2]; u[grepl('mean_city',names(u))]} ))
tab.s.p <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(1:5, FUN = function(i){ u<- summary(twoSLS_on_selected(data5, vd, m, pollutants[i],baseline[[1]][[i]]))$coefficients[,4]; u[grepl('mean_city',names(u))]} ))
a.s <- lapply(a, summary)
stargazer(a[[1]],a[[1]],a[[2]],a[[2]],
          # coefficients
          coef = list(a.s[[1]]$coefficients[grepl('mean_city',rownames(a.s[[1]]$coefficients)),1],
                      unlist(tab.s.c[[1]]),
                      a.s[[2]]$coefficients[grepl('mean_city',rownames(a.s[[2]]$coefficients)),1],
                      unlist(tab.s.c[[2]])),
          # standard errors
          se = list(a.s[[1]]$coefficients[grepl('mean_city',rownames(a.s[[1]]$coefficients)),2],
                   unlist(tab.s.se[[1]]),
                    a.s[[2]]$coefficients[grepl('mean_city',rownames(a.s[[2]]$coefficients)),2],
                    unlist(tab.s.se[[2]])),
          p = list(a.s[[1]]$coefficients[grepl('mean_city',rownames(a.s[[1]]$coefficients)),4],
                    unlist(tab.s.p[[1]]),
                   a.s[[2]]$coefficients[grepl('mean_city',rownames(a.s[[2]]$coefficients)),4],
                   unlist(tab.s.p[[2]])),omit.stat=c("ser","adj.rsq","rsq"),   covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4,  
          add.lines = list(c('Instruments',length(selected_set),paste0(lapply(baseline[[1]],length),collapse = ";"),length(selected_set),paste0(lapply(baseline[[1]],length),collapse = ";"))))


stargazer(a[[4]],a[[4]],a[[5]],a[[5]],
          # coefficients
          coef = list(a.s[[4]]$coefficients[grepl('mean_city',rownames(a.s[[4]]$coefficients)),1],
                      unlist(tab.s.c[[3]]),
                      a.s[[5]]$coefficients[grepl('mean_city',rownames(a.s[[5]]$coefficients)),1],
                      unlist(tab.s.c[[4]])),
          # standard errors
          se = list(a.s[[4]]$coefficients[grepl('mean_city',rownames(a.s[[1]]$coefficients)),2],
                    unlist(tab.s.se[[3]]),
                    a.s[[5]]$coefficients[grepl('mean_city',rownames(a.s[[2]]$coefficients)),2],
                    unlist(tab.s.se[[4]])),
          p = list(a.s[[4]]$coefficients[grepl('mean_city',rownames(a.s[[1]]$coefficients)),4],
                   unlist(tab.s.p[[3]]),
                   a.s[[5]]$coefficients[grepl('mean_city',rownames(a.s[[2]]$coefficients)),4],
                   unlist(tab.s.p[[4]])),omit.stat=c("ser","adj.rsq","rsq"),   covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4, add.lines = list(c('Instruments',length(selected_set),paste0(lapply(baseline[[1]],length),collapse = ";"),length(selected_set),paste0(lapply(baseline[[1]],length),collapse = ";"))))

length(selected_set)
lapply(baseline[[1]],length)

## MHT.
lapply(tab,FWER_models)


tab <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(1:5, FUN = function(i) twoSLS_on_selected(data5, vd, m, pollutants[i],c("FEP_invth_simple_mean","FEP_inv_hcl_mean" ))))

stargazer(tab, omit.stat=c("ser","adj.rsq","rsq"),   covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4, type = "text", title = "Only pollutant-specific instrument for each pollutant")

# Report in Table 7.
tab.s.c <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data5, vd, m, pollutants[i],c("FEP_invth_simple_mean","FEP_inv_hcl_mean" )))$coefficients[,1]; u[grepl('mean_city',names(u))]} ))
tab.s.se <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(1:5, FUN = function(i){ u<- summary(twoSLS_on_selected(data5, vd, m, pollutants[i],c("FEP_invth_simple_mean","FEP_inv_hcl_mean" )))$coefficients[,2]; u[grepl('mean_city',names(u))]} ))
tab.s.p <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(1:5, FUN = function(i){ u<- summary(twoSLS_on_selected(data5, vd, m, pollutants[i],c("FEP_invth_simple_mean","FEP_inv_hcl_mean" )))$coefficients[,4]; u[grepl('mean_city',names(u))]} ))
a.s <- lapply(a, summary)
stargazer(a[[1]],a[[1]],a[[2]],a[[2]],
          # coefficients
          coef = list(a.s[[1]]$coefficients[grepl('mean_city',rownames(a.s[[1]]$coefficients)),1],
                      unlist(tab.s.c[[1]]),
                      a.s[[2]]$coefficients[grepl('mean_city',rownames(a.s[[2]]$coefficients)),1],
                      unlist(tab.s.c[[2]])),
          # standard errors
          se = list(a.s[[1]]$coefficients[grepl('mean_city',rownames(a.s[[1]]$coefficients)),2],
                    unlist(tab.s.se[[1]]),
                    a.s[[2]]$coefficients[grepl('mean_city',rownames(a.s[[2]]$coefficients)),2],
                    unlist(tab.s.se[[2]])),
          p = list(a.s[[1]]$coefficients[grepl('mean_city',rownames(a.s[[1]]$coefficients)),4],
                   unlist(tab.s.p[[1]]),
                   a.s[[2]]$coefficients[grepl('mean_city',rownames(a.s[[2]]$coefficients)),4],
                   unlist(tab.s.p[[2]])),omit.stat=c("ser","adj.rsq","rsq"),   covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4, 
          add.lines = list(c('Instruments',length(selected_set),2,length(selected_set),2)))


stargazer(a[[4]],a[[4]],a[[5]],a[[5]],
          # coefficients
          coef = list(a.s[[4]]$coefficients[grepl('mean_city',rownames(a.s[[4]]$coefficients)),1],
                      unlist(tab.s.c[[3]]),
                      a.s[[5]]$coefficients[grepl('mean_city',rownames(a.s[[5]]$coefficients)),1],
                      unlist(tab.s.c[[4]])),
          # standard errors
          se = list(a.s[[4]]$coefficients[grepl('mean_city',rownames(a.s[[1]]$coefficients)),2],
                    unlist(tab.s.se[[3]]),
                    a.s[[5]]$coefficients[grepl('mean_city',rownames(a.s[[2]]$coefficients)),2],
                    unlist(tab.s.se[[4]])),
          p = list(a.s[[4]]$coefficients[grepl('mean_city',rownames(a.s[[1]]$coefficients)),4],
                   unlist(tab.s.p[[3]]),
                   a.s[[5]]$coefficients[grepl('mean_city',rownames(a.s[[2]]$coefficients)),4],
                   unlist(tab.s.p[[4]])),omit.stat=c("ser","adj.rsq","rsq"),   covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4,  
          add.lines = list(c('Instruments',length(selected_set),2,length(selected_set),2)))

length(selected_set)
lapply(baseline[[1]],length)


###################################################################################################################
############################## Robustness Exercices : Timing, Weathers, Min/Max etc.  #############################
###################################################################################################################


## Timing ##
# Leads/Lags
timing <- data[order(ville,date)]
cols <- c(pollutants, instruments_set_list[[1]])
leadscols <- paste0(cols, "_lead")
lagscols <- paste0(cols, "_lag")
lagscols2 <- paste0(cols, "_lag2")
timing[, (leadscols):= shift(.SD, -1 , fill = NA), .SDcols = cols, by = ville]
timing[, (lagscols):= shift(.SD, 1 , fill = NA), .SDcols = cols, by = ville]
timing[, (lagscols2):= shift(.SD, 2 , fill = NA), .SDcols = cols, by = ville]

### Respiratory Diseases. 
resp_pollutants <- c('FEP_O3_mean_city','FEP_CO_mean_city','FEP_SO2_mean_city','FEP_NO2_mean_city')
endog <- c(resp_pollutants,paste0(resp_pollutants,'_lag'),paste0(resp_pollutants, '_lag2'))
endog <- endog[order(endog)]
t0 <- twoSLS_on_selected(timing, 'FEP_rate_adm_C_10_Full_all', 
                         m, endog, 
                         c(selected_set,paste0(selected_set, '_lag'),paste0(selected_set,"_lag2")))
resp_pollutants <- c('FEP_O3_mean_city','FEP_SO2_mean_city')
endog <- c(paste0(resp_pollutants,'_lead'),c(resp_pollutants),paste0(resp_pollutants,'_lag'))
endog <- endog[order(endog)]
t1 <- twoSLS_on_selected(timing, 'FEP_rate_adm_C_10_Full_all', 
                         m, endog,c(selected_set,paste0(selected_set, '_lag'),paste0(pollutants, '_lag2')))

## Cardiovascular Diseases.
cardio_pollutants <- c("FEP_CO_mean_city")
endog <- c(cardio_pollutants,paste0(cardio_pollutants,'_lag'),paste0(cardio_pollutants, '_lag2'))
endog <- endog[order(endog)]
c0 <- twoSLS_on_selected(timing, 'FEP_rate_adm_C_09_Full_all', 
                         m, endog,
                         c(selected_set,paste0(selected_set, '_lag'),paste0(selected_set,"_lag2")))
endog <- c(paste0(cardio_pollutants,'_lead'),cardio_pollutants,paste0(cardio_pollutants,'_lag'))
endog <- endog[order(endog)]

c1 <- twoSLS_on_selected(timing, 'FEP_rate_adm_C_09_Full_all', 
                         m, endog, 
                         c(selected_set,paste0(selected_set, '_lag'),paste0(selected_set,"_lag2")))
## Deaths.
deaths_pollutants <- c('FEP_PM2_5_mean_city')
endog <- c(deaths_pollutants,paste0(deaths_pollutants,'_lag'),paste0(deaths_pollutants, '_lag2'))
endog <- endog[order(endog)]

d0 <- twoSLS_on_selected(timing, 'FEP_n_I_in_CertD',
                         m, endog, 
                         c(selected_set,paste0(selected_set, '_lag'),paste0(selected_set,"_lag2")))

endog <- c(paste0(deaths_pollutants,'_lead'),deaths_pollutants,paste0(deaths_pollutants,'_lag'),paste0(deaths_pollutants, '_lag2'))
endog <- endog[order(endog)]
d1 <- twoSLS_on_selected(timing, 'FEP_n_I_in_CertD', 
                         m, endog,
                         c(selected_set,paste0(selected_set, '_lag'),paste0(selected_set,"_lag2")))
deaths_pollutants <- c('FEP_SO2_mean_city')
endog <- c(deaths_pollutants,paste0(deaths_pollutants,'_lag'),paste0(deaths_pollutants, '_lag2'))
endog <- endog[order(endog)]
h0 <- twoSLS_on_selected(timing, 'FEP_n_J_in_CertD', 
                         m, endog, 
                         c(selected_set,paste0(selected_set, '_lag'),paste0(selected_set,"_lag2")))
endog <- c(paste0(deaths_pollutants,'_lead'),deaths_pollutants,paste0(deaths_pollutants,'_lag'),paste0(deaths_pollutants, '_lag2'))
endog <- endog[order(endog)]
h1 <- twoSLS_on_selected(timing, 'FEP_n_J_in_CertD', 
                         m, endog,
                         c(selected_set,paste0(selected_set, '_lag'),paste0(selected_set,"_lag2")))
stargazer(t0,t1, c0,c1, d0,d1,h0,h1, omit.stat=c("ser","adj.rsq","rsq"), 
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)

# Measurement in pollutant data. mean? Min/max #
b1 <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants_max, baseline_max))
b2 <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants_min, baseline_min))

stargazer(b1,omit.stat=c("ser","adj.rsq","rsq"), covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"), 
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)

stargazer(b2,omit.stat=c("ser","adj.rsq","rsq"), covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)
length(baseline_max)
length(baseline_min)

pollutants_nam <- gsub('FEP_','',pollutants)
pollutants_max_nam <- gsub('FEP_','',pollutants_max)
pollutants_min_nam <- gsub('FEP_','',pollutants_min)
apply(data[,..pollutants_nam], MARGIN = 2, FUN = function(x) sd(x, na.rm=T))
apply(data[,..pollutants_max_nam], MARGIN = 2, FUN = function(x) sd(x, na.rm=T))
apply(data[,..pollutants_min_nam], MARGIN = 2, FUN = function(x) sd(x, na.rm=T))

# PM2.5 vs PM10
d <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, setdiff(pollutants6, 'FEP_PM2_5_mean_city'), 
                                                            unique(c(baseline[[1]][[2]], baseline[[1]][[3]], 
                                                                     baseline[[1]][[4]], baseline[[1]][[5]],unlist(pm10_selection)))))
stargazer(a,d, type ='text', omit.stat=c("ser","adj.rsq","rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)
# Other set of instruments.
e <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants, unique(setdiff(unique(unlist(baseline[[2]])),m))))
f <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants, unique(setdiff(unique(unlist(baseline[[3]])),m))))
g <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants, unique(setdiff(unique(unlist(baseline[[4]])),m))))
h <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants, unique(setdiff(unique(unlist(baseline[[5]])),m))))

stargazer(a[[1]],e[[1]],f[[1]],g[[1]],h[[1]],type ='text', omit.stat=c("ser","adj.rsq","rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4, add.lines = list(c('Instruments', length(selected_set),length(unique(unlist(baseline[[2]]))), length(unique(unlist(baseline[[3]]))), length(unique(unlist(baseline[[4]]))), length(unique(unlist(baseline[[5]]))))))
stargazer(a[[2]],e[[2]],f[[2]],g[[2]],h[[2]],type ='text', omit.stat=c("ser","adj.rsq","rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4, add.lines = list(c('Instruments', length(selected_set),length(unique(unlist(baseline[[2]]))), length(unique(unlist(baseline[[3]]))), length(unique(unlist(baseline[[4]]))), length(unique(unlist(baseline[[5]]))))))
stargazer(a[[4]],e[[3]],f[[3]],g[[3]],h[[3]],type ='text', omit.stat=c("ser","adj.rsq","rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4, add.lines = list(c('Instruments', length(selected_set),length(unique(unlist(baseline[[2]]))), length(unique(unlist(baseline[[3]]))), length(unique(unlist(baseline[[4]]))), length(unique(unlist(baseline[[5]]))))))
stargazer(a[[5]],e[[4]],f[[4]],g[[4]],h[[4]],type ='text', omit.stat=c("ser","adj.rsq","rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4, add.lines = list(c('Instruments', length(selected_set),length(unique(unlist(baseline[[2]]))), length(unique(unlist(baseline[[3]]))), length(unique(unlist(baseline[[4]]))), length(unique(unlist(baseline[[5]]))))))


# Predict on large sample, run estimation on reduced sample.
h <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants, unique(unlist(baseline_sampleA))))

stargazer(a,h, type ='text', omit.stat=c("ser","adj.rsq","rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)

## Interacted terms (Response to Sebastien Roux).
pollutant_interaction <- c('FEP_PM2_5_CO_interacted', 'FEP_O3_CO_interacted', "FEP_O3_SO2_interacted","FEP_PM2_5_SO2_interacted")


## Mortality and CO. small interaction on n_J .. (not n_I) 

aint <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, c(pollutants, pollutant_interaction), c(selected_set, unique(unlist(baseline_interactions[[1]])))))
aint2 <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_SO2_mean_city", 'FEP_PM2_5_CO_interacted'), c(baseline[[1]][[1]],baseline[[1]][[2]],baseline[[1]][[5]], baseline_interactions[[1]][[1]])))
aint3 <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, c("FEP_PM2_5_mean_city","FEP_CO_mean_city", 'FEP_PM2_5_CO_interacted'), c(baseline[[1]][[1]],baseline[[1]][[2]], baseline_interactions[[1]][[1]])))
stargazer(aint[[4]],aint2[[4]],aint3[[5]],type ='text', omit.stat=c("ser","adj.rsq","rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)

aint3 <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, c(pollutants,  'FEP_O3_CO_interacted'), c(selected_set,  baseline_interactions[[1]][[2]] )))
aint4 <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, c(pollutants,  'FEP_O3_SO2_interacted'), c(selected_set,  baseline_interactions[[1]][[3]] )))
aint5 <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, c("FEP_O3_mean_city","FEP_SO2_mean_city","FEP_CO_mean_city",  'FEP_O3_CO_interacted'), c(baseline[[1]][[2]],baseline[[1]][[3]], baseline[[1]][[5]],  baseline_interactions[[1]][[2]] )))
aint6 <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, c("FEP_O3_mean_city","FEP_SO2_mean_city", 'FEP_O3_SO2_interacted'), c(baseline[[1]][[3]], baseline[[1]][[5]],  baseline_interactions[[1]][[3]] )))
aint7 <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, c("FEP_PM2_5_mean_city","FEP_SO2_mean_city", 'FEP_PM2_5_SO2_interacted'), c(baseline[[1]][[3]], baseline[[1]][[5]],  baseline_interactions[[1]][[4]] )))
stargazer(aint7,type ='text', omit.stat=c("ser","adj.rsq","rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)
stargazer(aint[[1]], aint3[[1]], aint4[[1]], aint5[[1]], aint6[[1]],type ='text', omit.stat=c("ser","adj.rsq","rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)



## ----  4 - Pollutants Models. ---- ## 
modelList <- list(c("FEP_CO_mean_city","FEP_O3_mean_city","FEP_SO2_mean_city","FEP_NO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_O3_mean_city","FEP_SO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_O3_mean_city","FEP_NO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_NO2_mean_city","FEP_SO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_SO2_mean_city","FEP_O3_mean_city","FEP_NO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_SO2_mean_city","FEP_O3_mean_city","FEP_NO2_mean_city","FEP_CO_mean_city"))

# Same set of pooled instruments as baseline five pollutant model.
for(vd in outcomes){
  
  c <- lapply(1:length(modelList), FUN = function(s) twoSLS_on_selected(data, vd, m, modelList[[s]], selected_set))
  stargazer(c, omit.stat=c("ser","adj.rsq","rsq"), type ='text',
            omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report=('vc*p'),
            dep.var.labels.include = F,digits=4)
  
  

}


## ---- Increasing the # of pollutants in models ---- ## 
modelList <- list(c('FEP_PM2_5_mean_city','FEP_SO2_mean_city'), c('FEP_CO_mean_city','FEP_O3_mean_city'), 
                  c('FEP_PM2_5_mean_city','FEP_O3_mean_city',"FEP_CO_mean_city"), 
                   c("FEP_CO_mean_city","FEP_O3_mean_city","FEP_SO2_mean_city","FEP_NO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_NO2_mean_city","FEP_O3_mean_city","FEP_SO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_O3_mean_city","FEP_NO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_O3_mean_city","FEP_NO2_mean_city","FEP_SO2_mean_city"))

# Same set of pooled instruments as baseline five pollutant model.

for(vd in outcomes[c(1,2,4,5)]){
  
  tab.s.c <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data, vd, m, pollutants[i],selected_set))$coefficients[,1]; u[grepl('mean_city',names(u))]} )
  tab.s.se <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data, vd, m, pollutants[i],selected_set))$coefficients[,2]; u[grepl('mean_city',names(u))]} )
  tab.s.p <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data, vd, m, pollutants[i],selected_set))$coefficients[,4]; u[grepl('mean_city',names(u))]} )
  
  c <- lapply(1:length(modelList), FUN = function(s) twoSLS_on_selected(data, vd, m, modelList[[s]], selected_set))
  
  stargazer(c[[length(modelList)]],c, omit.stat=c("ser","adj.rsq","rsq"), 
            # coefficients
            coef = list(unlist(tab.s.c)),
            # standard errors
            se = list(unlist(tab.s.se)),
            p = list(unlist(tab.s.p)),
            covariate.labels = c("PM2.5","CO","O3","NO2","SO2"),
            omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
            dep.var.labels.include = F,digits=4)
}




for(vd in c('FEP_n_I','FEP_n_J','FEP_n_K')){
  
  c <- lapply(1:length(modelList), FUN = function(s) twoSLS_on_selected(data, vd, m, modelList[[s]], selected_set))
  stargazer(c, omit.stat=c("ser","adj.rsq","rsq"), type ='text',
            omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report=('vc*p'),
            dep.var.labels.include = F,digits=4)
}

# Selection on same sample as estimation, only instruments from pollutants entering the models.

modelList <- list(c("FEP_O3_mean_city","FEP_SO2_mean_city"),
                  c("FEP_NO2_mean_city","FEP_O3_mean_city"),
                  c('FEP_PM2_5_mean_city','FEP_SO2_mean_city'),
                  c('FEP_PM2_5_mean_city','FEP_CO_mean_city'),
                  c('FEP_CO_mean_city','FEP_O3_mean_city'),
                  c("FEP_CO_mean_city","FEP_O3_mean_city","FEP_SO2_mean_city","FEP_NO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_O3_mean_city","FEP_SO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_O3_mean_city","FEP_NO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_NO2_mean_city","FEP_SO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_SO2_mean_city","FEP_O3_mean_city","FEP_NO2_mean_city"))



monotomulti <- NULL
i <- 1
for(mod in modelList){
  
  sample <- complete.cases(data[,..mod])
  selec <- lapply(mod, FUN = function(x) selected_instruments(data[sample], m, instruments_set_list[[1]], x, clustervar = 'month_year_ville', baseline_only = TRUE))
  selec <- unique(unlist(lapply(selec, FUN = function(x) setdiff(x[[1]],m))))
  
  monotomulti[[i]] <-  lapply(outcomes,FUN = function(vd) twoSLS_on_selected(data, vd, m, mod, selec))
  stargazer(monotomulti[[i]],omit.stat=c("ser","adj.rsq","rsq"), type ='text',
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report=('vc*p'),
          dep.var.labels.include = F,digits=4)
  i <- i+1
}
  

stargazer(lapply(monotomulti, FUN = function(x) x[[1]]),a[[1]], omit.stat=c("ser","adj.rsq","rsq"), type ='text',
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report=('vc*p'),
          dep.var.labels.include = F,digits=4)


stargazer(lapply(monotomulti, FUN = function(x) x[[2]]),a[[2]], omit.stat=c("ser","adj.rsq","rsq"), type ='text',
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report=('vc*p'),
          dep.var.labels.include = F,digits=4)

stargazer(lapply(monotomulti, FUN = function(x) x[[3]]),a[[3]], omit.stat=c("ser","adj.rsq","rsq"), type ='text',
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report=('vc*p'),
          dep.var.labels.include = F,digits=4)

stargazer(lapply(monotomulti, FUN = function(x) x[[4]]),a[[4]], omit.stat=c("ser","adj.rsq","rsq"), type ='text',
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), report=('vc*p'),
          dep.var.labels.include = F,digits=4)


## --------- Table 16-18: sample for estimation and selection: do they matter ? ---- ##

a <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants, selected_set))
b <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants, unique(unlist(baseline_sampleA))))
stargazer(b, omit.stat=c("ser","adj.rsq","rsq"),covariate.labels = c('PM2.5','CO','O3',"NO2","SO2"), type ='text',
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"), #report=('vc*p'),
          dep.var.labels.include = F,digits=4)

tab <- lapply(outcomes, FUN = function(vd) list( lapply(pollutants, FUN = function(x) twoSLS_on_selected(data, vd, m, x, unique(unlist(baseline_sampleA)))),
                                                 lapply(pollutants, FUN = function(x) twoSLS_on_selected(data5, vd, m, x, unique(unlist(baseline_sampleA)))),
                                                 lapply(pollutants, FUN = function(x) twoSLS_on_selected(data5, vd, m, x, unique(unlist(baseline[[1]]))))))
paste0(sapply(baseline_sampleA,length), collapse = ";")
paste0(sapply(baseline[[1]],length), collapse = ";")

for(vd in outcomes){
  ref <-  twoSLS_on_selected(data5, vd, m, pollutants, selected_set)
  tab1.s.c <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data, vd, m, pollutants[i],baseline_sampleA[[i]]))$coefficients[,1]; u[grepl('mean_city',names(u))]} )
  tab1.s.se <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data, vd, m, pollutants[i],baseline_sampleA[[i]]))$coefficients[,2]; u[grepl('mean_city',names(u))]} )
  tab1.s.p <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data, vd, m, pollutants[i],baseline_sampleA[[i]]))$coefficients[,4]; u[grepl('mean_city',names(u))]} )
  
  tab2.s.c <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data5, vd, m, pollutants[i],baseline_sampleA[[i]]))$coefficients[,1]; u[grepl('mean_city',names(u))]} )
  tab2.s.se <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data5, vd, m, pollutants[i],baseline_sampleA[[i]]))$coefficients[,2]; u[grepl('mean_city',names(u))]} )
  tab2.s.p <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data5, vd, m, pollutants[i],baseline_sampleA[[i]]))$coefficients[,4]; u[grepl('mean_city',names(u))]} )
  
  tab3.s.c <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data5, vd, m, pollutants[i],baseline[[1]][[i]]))$coefficients[,1]; u[grepl('mean_city',names(u))]} )
  tab3.s.se <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data5, vd, m, pollutants[i],baseline[[1]][[i]]))$coefficients[,2]; u[grepl('mean_city',names(u))]} )
  tab3.s.p <- lapply(1:5, FUN = function(i){ u<-summary(twoSLS_on_selected(data5, vd, m, pollutants[i],baseline[[1]][[i]]))$coefficients[,4]; u[grepl('mean_city',names(u))]} )
  
  stargazer(ref,ref,ref,ref, omit.stat=c("ser","adj.rsq","rsq"), 
            # coefficients
            coef = list(unlist(tab1.s.c),unlist(tab2.s.c),unlist(tab3.s.c)),
            # standard errors
            se = list(unlist(tab1.s.se),unlist(tab2.s.se),unlist(tab3.s.se)),
            p = list(unlist(tab1.s.p),unlist(tab2.s.p),unlist(tab3.s.p)),
            covariate.labels = c("PM2.5","CO","O3","NO2","SO2"),
            omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
            dep.var.labels.include = F,digits=4)
  
  
}


## 
load("data/selection_step_variant_W.Rdata")
robust <- lapply(1:length(res_robW), FUN = function(i) setdiff(unique(unlist(lapply(res_robW[[i]], FUN=function(x) x$selected_robustunpC))), weather_list_rob[[i]]))
robust2 <- lapply(1:length(res_robW_largersample), FUN = function(i) setdiff(unique(unlist(lapply(res_robW_largersample[[i]], FUN=function(x) x$selected_robustunpC))), weather_list_rob[[i]]))
sapply(robust, FUN = function(x) length(setdiff(x,m)))
sapply(robust2, FUN = function(x) length(setdiff(x,m)))

sample <- apply(data5[,..m][,.(na = !is.na(.SD))], MARGIN = 1 , all)

a <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants, selected_set))
aw <- lapply(outcomes, FUN = function(vd) lapply(c(1:8), FUN = function(i) twoSLS_on_selected(data5[sample], vd, weather_list_rob[[i]], pollutants, setdiff(robust[[i]],m))))

stargazer(a[[1]],aw[[1]][[1]],aw[[1]][[8]], aw[[1]][[4]], 
          a[[2]],aw[[2]][[1]],aw[[2]][[8]], aw[[2]][[4]],
          a[[4]],aw[[4]][[1]],aw[[4]][[8]], aw[[4]][[4]],
          a[[5]],aw[[5]][[1]],aw[[5]][[8]], aw[[5]][[4]], covariate.labels = c('PM2.5','CO','O3',"NO2","SO2"),
          omit.stat=c("ser","adj.rsq","rsq"), type = "text",
          dep.var.labels.include = F,digits=4)


stargazer(a[[1]],aw[[1]][[1]],aw[[1]][[8]], aw[[1]][[4]],a[[2]],aw[[2]][[1]],aw[[2]][[8]], aw[[2]][[4]], covariate.labels = c('PM2.5','CO','O3',"NO2","SO2"),
          omit.stat=c("ser","adj.rsq","rsq"), omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4,
          add.lines = list(c('', length(selected_set),length(robust[[1]]), length(robust[[8]]), length(robust[[4]]), length(selected_set),length(robust[[1]]), length(robust[[8]]), length(robust[[4]]))))
stargazer(a[[4]],aw[[4]][[1]],aw[[4]][[8]], aw[[4]][[4]],
  a[[5]],aw[[5]][[1]],aw[[5]][[8]], aw[[5]][[4]], covariate.labels = c('PM2.5','CO','O3',"NO2","SO2"),
          omit.stat=c("ser","adj.rsq","rsq"), omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4,
          add.lines = list(c('', length(selected_set),length(robust[[1]]), length(robust[[8]]), length(robust[[4]]), length(selected_set),length(robust[[1]]), length(robust[[8]]), length(robust[[4]]))))


stargazer(a[[1]],aw[[1]][[1]],aw[[1]][[8]], aw[[1]][[4]], 
          omit.stat=c("ser","adj.rsq","rsq"), type = "text",
          dep.var.labels.include = F,digits=4)
stargazer(a[[2]],aw[[2]][[1]],aw[[2]][[8]], aw[[2]][[4]],
          omit.stat=c("ser","adj.rsq","rsq"), type = "text",
          dep.var.labels.include = F,digits=4)
stargazer(a[[4]],aw[[4]][[1]],aw[[4]][[8]], aw[[4]][[4]],
          omit.stat=c("ser","adj.rsq","rsq"), type = "text",
          dep.var.labels.include = F,digits=4)
stargazer(a[[5]],aw[[5]][[1]],aw[[5]][[8]], aw[[5]][[4]],
          omit.stat=c("ser","adj.rsq","rsq"), type = "text",
          dep.var.labels.include = F,digits=4)

## Justify not to bins to much (in referee responses):
stargazer(a[[5]], aw[[5]][[4]], a[[4]],aw[[4]][[4]],a[[2]], aw[[2]][[4]],a[[1]], aw[[1]][[4]],
          omit.stat=c("ser","adj.rsq","rsq"), type = "text",
          dep.var.labels.include = F,digits=4)


rob1 <- twoSLS_on_selected(data5[sample], "FEP_n_I_in_CertD", weather_list_rob[[4]], c('FEP_PM2_5_mean_city','FEP_SO2_mean_city'), robust[[4]])
rob2 <- twoSLS_on_selected(data5[sample], "FEP_n_J_in_CertD", weather_list_rob[[4]], c('FEP_PM2_5_mean_city','FEP_SO2_mean_city'), robust[[4]])
rob3 <- twoSLS_on_selected(data, "FEP_n_I_in_CertD", weather_list_rob[[4]], c('FEP_PM2_5_mean_city','FEP_SO2_mean_city'), unlist(res_robW_largersample[[4]][[1]]))
rob4 <- twoSLS_on_selected(data, "FEP_n_J_in_CertD", weather_list_rob[[4]], c('FEP_PM2_5_mean_city','FEP_SO2_mean_city'), unlist(res_robW_largersample[[4]][[1]]))
stargazer(a[[4]], rob2, rob4,a[[5]], rob1, rob3, omit.stat=c("ser","adj.rsq","rsq"), type ="text",report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)

## Add one by one weather controls. 
w <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(list(c("FEP_humidite","FEP_temperature" ),
                                                                 c("FEP_ensoleil_glo","FEP_humidite","FEP_temperature" ),
                                                                 c("FEP_ensoleil_glo","FEP_humidite","FEP_temperature","FEP_precipitations"),
                                                                 c("FEP_ensoleil_glo","FEP_humidite","FEP_temperature","FEP_precipitations","FEP_ff"),m ), FUN = function(n) twoSLS_on_selected(data5, vd, n, pollutants, selected_set)))
stargazer(w[[1]],w[[2]], omit.stat=c("ser","adj.rsq","rsq"), 
          covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)

stargazer(w[[3]],w[[4]], omit.stat=c("ser","adj.rsq","rsq"), 
          covariate.labels = c('PM2.5',"CO","O3","NO2","SO2"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F,digits=4)

w <- lapply(outcomes[c(1,2,4,5)], FUN = function(vd) lapply(list(names(data)[grepl("bin",names(data))],m ), FUN = function(n) twoSLS_on_selected(data5, vd, n, pollutants, selected_set)))



### ===================== By Age ================================== ##

resp <- c('FEP_rate_adm_C_10_Full_0-4','FEP_rate_adm_C_10_Full_5-14',
          'FEP_rate_adm_C_10_Full_15-59','FEP_rate_adm_C_10_Full_60-79','FEP_rate_adm_C_10_Full_80-PP')
res5 <- lapply(resp, FUN = function(x)  twoSLS_on_selected(data5, x, m, pollutants,  selected_set))
res2 <- lapply(resp, FUN = function(x)  twoSLS_on_selected(data, x, m, c('FEP_O3_mean_city','FEP_SO2_mean_city'),  unique(c(baseline_sampleA[[3]],baseline_sampleA[[5]]))))

length(unique(c(baseline_sampleA[[3]],baseline_sampleA[[5]])))

plot_res <- function(n_pol, pol, res_list, outc, level = 0.95, labels = c("0-4","5-14","15-59","60-79",">=80")){
  ci <- lapply(res_list, FUN = function(res) sapply(res, FUN = function(x) confint(x, level = level)[grepl(pol, rownames(confint(x)))] ))
  pe <- lapply(res_list, FUN = function(res) sapply(res, FUN = function(x) x$coefficients[grepl(pol, rownames(x$coefficients))]))
  ci <- Reduce(rbind, lapply(ci, FUN = function(x) {dt <-  data.frame(x); names(dt) <- outc; dt$y <- c('min','max'); dt}))
  ci$N_Pollutants <- rep(n_pol, each = 2)
  pe <- Reduce(rbind, lapply(pe, FUN = function(x) {dt <-  data.frame(t(x)); names(dt) <- outc; dt$y <- c('estimate'); dt}))
  pe$N_Pollutants <- n_pol
  toplot <- melt(rbind(ci,pe), id.var = c('y','N_Pollutants'))
  toplot <- dcast(toplot, N_Pollutants + variable ~ y)
  toplot$N_Pollutants <- as.factor(toplot$N_Pollutants)
  toplot$aux <- as.vector(matrix(seq(1,length(outc)), nrow= length(outc), ncol = length(n_pol))+ t(matrix(seq(0,(length(n_pol)-1)*0.1,0.1), ncol= length(outc), nrow = length(n_pol)))) 
  print(ggplot(toplot, aes(aux, estimate, ymin = min, ymax = max, col = N_Pollutants, shape = N_Pollutants)) + geom_pointrange() + ggtitle(pol) +
    scale_x_continuous("Age Group",breaks=seq(1.1,length(outc) +0.1,1), labels= labels) + ylab('Causal Estimate') +theme_bw() + 
      theme(legend.title = element_blank(), legend.position = 'bottom') + scale_color_manual(values = c('blue','#2980b9','black')))
  
}

pdf('O3age.pdf',3,4)
plot_res(c(5,2), 'O3', list(res5,res2), resp, level = 0.9)
dev.off()
pdf('SO2age.pdf',3,4)
plot_res(c(5,2), 'SO2', list(res5,res2), resp, level = 0.9)
dev.off()
stargazer(res5,res2, omit.stat=c("ser","adj.rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F)


cardio <-  c('FEP_rate_adm_C_09_Full_0-14','FEP_rate_adm_C_09_Full_15-59','FEP_rate_adm_C_09_Full_60-79','FEP_rate_adm_C_09_Full_80-PP')

cardiores5 <- lapply(cardio, FUN = function(x)  twoSLS_on_selected(data5, x, m, pollutants,  selected_set))
cardiores1 <- lapply(cardio, FUN = function(x)  twoSLS_on_selected(data, x, m, 'FEP_CO_mean_city',  baseline_sampleA[[2]]))

pdf('COage.pdf',2.5,3.5)
plot_res(c(5,1), 'CO', list(cardiores5,cardiores1), cardio, level = 0.9, labels = c("0-14","15-59","60-79",">=80"))
dev.off()

stargazer(cardiores5,cardiores1, omit.stat=c("ser","adj.rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F)
length(baseline_sampleA[[2]])

death <- c('FEP_rate_dec_0-14','FEP_rate_dec_15-59','FEP_rate_dec_60-79','FEP_rate_dec_80-PP')
deathres <- lapply(death, FUN = function(x)  twoSLS_on_selected(data5, x, m, pollutants,  selected_set))
deathres2 <- lapply(death, FUN = function(x)  twoSLS_on_selected(data, x, m, c('FEP_PM2_5_mean_city','FEP_SO2_mean_city'),  unique(c(baseline_sampleA[[1]], baseline_sampleA[[5]]))))

length(unique(c(baseline_sampleA[[1]],baseline_sampleA[[5]])))

stargazer(deathres, deathres2, omit.stat=c("ser","adj.rsq"),  report=('vc*p'),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F)


pdf('PM2_5age_death.pdf',2.5,3.5)
plot_res(c(5,2), 'PM2_5', list(deathres,deathres2), death, level = 0.9, labels = c('0-14',"15-59","60-79",">=80"))
dev.off()

pdf('SO2age_death.pdf',2.5,3.5)
plot_res(c(5,2), 'SO2', list(deathres,deathres2), death, level = 0.9, labels = c('0-14',"15-59","60-79",">=80"))
dev.off()

resp_age <- gsub('FEP_','',resp)
stargazer(apply(data[,..resp_age], MARGIN = 2, FUN= function(x) round(mean(x),1)), summary = FALSE , digits = 1)
stargazer(res2, omit.stat=c("ser","adj.rsq"), type = 'text',
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F)
cardio_age <- gsub('FEP_','',cardio)
stargazer(apply(data[,..cardio_age], MARGIN = 2, FUN= function(x) round(mean(x),1)), summary = FALSE, digits = 1)
stargazer(cardiores1, omit.stat=c("ser","adj.rsq"),
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F)
death_age <- gsub('FEP_','',death)
stargazer(apply(data[,..death_age], MARGIN = 2, FUN= function(x) round(mean(x),1)), summary = FALSE, digits = 1)
stargazer(deathres2, omit.stat=c("ser","adj.rsq"), 
          omit=c("Cons","neige","vergl","precip","temp","ff","humid","ensol"),
          dep.var.labels.include = F)



### ===================== Simple IV =============================== ##

instruments_IV <- list('inv_hcl_mean', 'invth_simple_mean', c('invth_simple_mean','inv_hcl_mean'))
m<- gsub('FEP_','',m)
unid_pol <- 'pollution_index_pca_imputed_F1'
FE <-  "as.factor(day)*as.factor(ville)+as.factor(month_year_ville)"
# Table 5 : First Stage
fs0 <- lapply(instruments_IV, FUN = function(inst) lapply(list(data,data5), FUN = function(sample) reg_FE_clustered(sample, x2 =  m, y = unid_pol , x = inst, FE =FE)))
fs <- lapply(instruments_IV, FUN = function(inst) lapply(c('PM2_5_mean_city',
               'PM10_mean_city',
               'NO2_mean_city',
               'O3_mean_city',
               'CO_mean_city',
               'SO2_mean_city'), FUN = function(pol) reg_FE_clustered(data, x2 = m, y = pol , x = inst, FE = FE)))


stargazer(fs0[[3]],fs[[3]],omit.stat=c("ser","F","adj.rsq"), digits = 1,
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))


fs0 <- lapply(list(data,data5), FUN = function(sample) reg_FE_clustered(sample, x2 =  m, y = unid_pol , x = c('scaled_ipblh','scaled_invth'), FE = FE))
fs <- lapply(c('PM2_5_mean_city_scaled',
               'PM10_mean_city_scaled',
               'NO2_mean_city_scaled',
               'O3_mean_city_scaled',
               'CO_mean_city_scaled',
               'SO2_mean_city_scaled'), FUN = function(pol) reg_FE_clustered(data, x2 = m, y = pol , x =c('scaled_ipblh','scaled_invth'), FE = FE))

stargazer(fs0,fs,omit.stat=c("ser","F","adj.rsq"), digits = 2, 
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))


# Table 7 : Reduced form.
rf0 <-  lapply(instruments_IV, FUN = function(inst) lapply(outcomes,  FUN = function(outcome)  reg_FE_clustered(data, x2 = m , y = outcome, x = inst, FE = FE)))
rf <-   lapply(instruments_IV, FUN = function(inst) lapply(outcomes, FUN = function(outcome) reg_FE_clustered(data5, x2 = m , y = outcome, x =  inst, FE = FE)))

stargazer(rf0[[1]],omit.stat=c("ser","F","adj.rsq"), digits = 2,
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))
stargazer(rf0[[2]],omit.stat=c("ser","F","adj.rsq"), digits = 2,
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))
stargazer(rf[[1]],omit.stat=c("ser","F","adj.rsq"), digits = 2, type = "text",
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))
stargazer(rf[[2]],omit.stat=c("ser","F","adj.rsq"), digits = 2, type = "text",
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))

# Show totals in Reduced Form
rf0 <-  lapply(instruments_IV, FUN = function(inst) lapply(c("FEP_rate_emergencies_tot","rate_dec_all"),  FUN = function(outcome)  reg_FE_clustered(data, x2 = m , y = outcome, x = inst, FE = FE)))
rf <-   lapply(instruments_IV, FUN = function(inst) lapply(c("FEP_rate_emergencies_tot","rate_dec_all"), FUN = function(outcome) reg_FE_clustered(data5, x2 = m , y = outcome, x =  inst, FE = FE)))
stargazer(rf0,omit.stat=c("ser","F","adj.rsq"), digits = 2, type = "text",
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))
stargazer(rf[[1]],omit.stat=c("ser","F","adj.rsq"), digits = 2, covariate.labels = "IPBLH",
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))


rf0 <-  lapply(instruments_IV, FUN = function(inst) lapply(c("FEP_n_I","FEP_n_J","FEP_n_K"),  FUN = function(outcome)  reg_FE_clustered(data, x2 = m , y = outcome, x = inst, FE = FE)))
rf <-   lapply(instruments_IV, FUN = function(inst) lapply(c("FEP_n_I","FEP_n_J","FEP_n_K"), FUN = function(outcome) reg_FE_clustered(data5, x2 = m , y = outcome, x =  inst, FE = FE)))
stargazer(rf0[[1]],omit.stat=c("ser","F","adj.rsq"), digits = 2, type = "text",
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))
stargazer(rf[[1]],omit.stat=c("ser","F","adj.rsq"), digits = 2, type = "text",
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))

## Weather: introducing FE?
rf <-   lapply(outcomes, FUN = function(outcome) lapply(c(FE,"0"), 
                                                        FUN = function(FE) reg_FE_clustered(data5, x2 = m[2:10] , y = outcome, x =  m[1], FE = FE)))
stargazer(rf,omit.stat=c("ser","F","adj.rsq"), digits = 2, type = "text",
          dep.var.labels.include = F)

# All Right

# IV (not shown)
iv <- lapply(c("FEP_rate_emergencies_tot",outcomes[1:3],"rate_dec_all",outcomes[4:6]), FUN = function(vd) lapply(list(data,data5), FUN = function(sample) twoSLS_on_selected(sample, vd, m, unid_pol, instruments_IV[[3]] , FE = FE)))
stargazer(iv,omit.stat=c("ser","F","adj.rsq"), digits = 2, type = 'text',
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))
iv <- lapply(c(  "FEP_rate_emergencies_tot","FEP_rate_emergencies_9_10","FEP_rate_emergencies_not_9_10","rate_adm_C_10_Full_all","rate_adm_C_09_Full_all","rate_adm_C_11_Full_all"), FUN = function(vd) lapply(list(data,data5), FUN = function(sample) twoSLS_on_selected(sample, vd, m, unid_pol, instruments_IV[[3]] , FE = FE)))
stargazer(iv,omit.stat=c("ser","F","adj.rsq"), digits = 2, type = 'text',
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))
# May justify that there is a sizeable effect on total emergencies - if anything, the effect is stronger. See with single p models?
# Issue: A stronger effect on "other emergencies" (not looked at)
iv <- lapply(c("FEP_rate_emergencies_tot",outcomes[1:2],"rate_dec_all",outcomes[4:5]), FUN = function(vd) twoSLS_on_selected(data5, vd, m, unid_pol, instruments_IV[[3]] , FE = FE))
stargazer(iv[[1]],iv[[4]],omit.stat=c("ser","F","adj.rsq"), digits = 2,
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))



# Table 13 : OLS estimates
ols <- lapply(outcomes, FUN = function(vd) reg_FE_clustered(data, x2 = m, y = vd, x = unid_pol))
stargazer(ols,omit.stat=c("ser","F","adj.rsq"), digits = 2, type = "text",
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))

# Table 14: 
ivMonoP0 <- lapply(outcomes, FUN = function(vd) lapply(pollutants6, FUN = function(pol) twoSLS_on_selected(data, vd, m, pol,  instruments_IV[[3]], FE = FE)))
stargazer(ivMonoP0[[1]],ivMonoP0[[2]],omit.stat=c("ser","F","adj.rsq"), digits = 2, covariate.labels = c('PM10',"PM2.5","CO","O3","NO2","SO2"), type = "text",
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))
stargazer(ivMonoP0[[4]],ivMonoP0[[5]],omit.stat=c("ser","F","adj.rsq"), digits = 2, covariate.labels = c('PM10',"PM2.5","CO","O3","NO2","SO2"), type = "text",
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))


### Other set of fixed effects (and weather controls?)
iv <- lapply(gsub('FEP_','',outcomes), FUN = function(vd) twoSLS_on_selected(data5, vd,gsub('FEP_','',m), gsub('FEP_','',pollutants), gsub('FEP_','',selected_set) , FE = "as.factor(day)*as.factor(ville)+as.factor(month_year_ville)"))
stargazer(iv,omit.stat=c("ser","F","adj.rsq"), digits = 2, type = 'text',
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))

data5$month <- month(data5$date)
data5$year<-year(data5$date)
iv <- lapply(gsub('FEP_','',outcomes), FUN = function(vd) twoSLS_on_selected(data5, vd,gsub('FEP_','',m), gsub('FEP_','',pollutants), gsub('FEP_','',selected_set) , FE = "as.factor(day)+as.factor(ville)+as.factor(month_year)"))
stargazer(iv,omit.stat=c("ser","F","adj.rsq"), digits = 2, type = 'text',
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))

iv <- lapply(gsub('FEP_','',outcomes), FUN = function(vd) twoSLS_on_selected(data5, vd,gsub('FEP_','',m), gsub('FEP_','',pollutants), gsub('FEP_','',selected_set) , FE = "as.factor(day)+as.factor(ville)+as.factor(month)+as.factor(year)"))
stargazer(iv,omit.stat=c("ser","F","adj.rsq"), digits = 2, type = 'text',
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))

iv <- lapply(gsub('FEP_','',outcomes), FUN = function(vd) twoSLS_on_selected(data5, vd,gsub('FEP_','',m), gsub('FEP_','',pollutants), gsub('FEP_','',selected_set) , FE = "as.factor(day)+as.factor(ville)+as.factor(month)+as.factor(year)"))
stargazer(iv,omit.stat=c("ser","F","adj.rsq"), digits = 2, type = 'text',
          dep.var.labels.include = F, omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))

