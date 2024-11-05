load("data/base_IV_purged_FE.Rdata")

setDT(data)

data$month_year<-factor(substr(data$date,1,7))
data$month_year_ville<-factor(paste0(data$month_year,data$ville)) 
data$FEP_ff_2 <-data$FEP_ff^2
data$FEP_temperature_2 <-data$FEP_temperature^2
data$FEP_precipitations_2 <-data$FEP_precipitations^2
data[,all5in:=!is.na(PM2_5_mean_city) & !is.na(CO_mean_city) & !is.na(SO2_mean_city) & !is.na(NO2_mean_city) & !is.na(O3_mean_city)]

data5<-data[all5in==T]

# Weather
m<-c("FEP_neige_sol","FEP_verglas_gel_sol","FEP_neige_couvrant_tout_sol","FEP_precipitations","FEP_precipitations_2",
     "FEP_temperature","FEP_temperature_2","FEP_ff","FEP_ff_2","FEP_humidite","FEP_ensoleil_glo")

## ====================== Lasso selection ============================= ##

load("data/selection_step.Rdata")
baseline <- lapply(1:length(res), FUN = function(i) unlist(lapply(res[[i]], FUN=function(x) x$selected_robustunpC)))

ci_df <- function(regs, pollutants){
  ci <- data.frame(Reduce(rbind,lapply(regs,FUN = function(x) {cint <- confint(x); cint <- cint[grepl('mean_city',rownames(cint)),]; rownames(cint) <- gsub('.FEP_','',gsub('_mean_city.fit..','',rownames(cint))); cint})))
  names(ci) <- c('cim','ciM')
  ci[,'age'] <- rep(gsub('FEP_rate_decAGE_','',gsub('FEP_rate_adm_C_.._Full_','',labels)), each = dim(ci)[1]/length(labels))
  ci[,'polluant'] <- rep(gsub('.?FEP_','',gsub('_mean_city(.fit..)?','',pollutants)), length(labels))
  ci[,'estimate'] <- unlist(lapply(regs,FUN = function(x) {coef <- coefficients(x); coef[grepl('mean_city',names(coef))] }))
  ci[,'N Pollutants'] <- as.character(length(pollutants))
  rownames(ci) <- NULL
  ci
}

## ----- By Age ------- ##
labels <- c("FEP_rate_adm_C_09_Full_all", "FEP_rate_adm_C_09_Full_0-4",  "FEP_rate_adm_C_09_Full_5-14",
            "FEP_rate_adm_C_09_Full_15-59",  "FEP_rate_adm_C_09_Full_60-64", "FEP_rate_adm_C_09_Full_65-69" ,"FEP_rate_adm_C_09_Full_70-74", "FEP_rate_adm_C_09_Full_75-79", "FEP_rate_adm_C_09_Full_80-PP")

modelList <- list("FEP_CO_mean_city",
                  c("FEP_CO_mean_city","FEP_PM2_5_mean_city"),
                  c("FEP_CO_mean_city","FEP_PM2_5_mean_city","FEP_O3_mean_city"),  
                  c("FEP_CO_mean_city","FEP_PM2_5_mean_city","FEP_O3_mean_city","FEP_SO2_mean_city"), 
                  c("FEP_CO_mean_city","FEP_PM2_5_mean_city","FEP_O3_mean_city","FEP_SO2_mean_city","FEP_NO2_mean_city"))

c09 <- lapply(modelList, FUN = function(poll) lapply(labels, FUN = function(vd) twoSLS_on_selected(data, vd, m, poll, baseline[[1]])))

ci <- Reduce(rbind,lapply(2:length(modelList), FUN = function(i) ci_df(c09[[i]], modelList[[i]]))) 
ci$age <- ifelse(ci$age == '80-PP', '>= 80',ci$age) 
ci$age <- factor(ci$age, levels = c('0-4','5-14','15-59','60-64','65-69','70-74','75-79','>= 80','all')) 
ggplot(ci[ci$polluant == "CO",],aes(age,estimate, ymin = cim, ymax = ciM, col = `N Pollutants`))+geom_pointrange(position = position_dodge(width = 0.3)) + facet_wrap(~polluant, scales = "free_y") + scale_color_manual(values = c("blue","#3A4F8F","grey","black")) + theme_bw() + ylab('Estimate') + xlab("Age Group")


labels <- c("FEP_rate_adm_C_10_Full_all", "FEP_rate_adm_C_10_Full_nn",  "FEP_rate_adm_C_10_Full_29j-1",  "FEP_rate_adm_C_10_Full_0-4",  "FEP_rate_adm_C_10_Full_5-14",
            "FEP_rate_adm_C_10_Full_15-59",  "FEP_rate_adm_C_10_Full_60-64", "FEP_rate_adm_C_10_Full_65-69" ,"FEP_rate_adm_C_10_Full_70-74", "FEP_rate_adm_C_10_Full_75-79", "FEP_rate_adm_C_10_Full_80-PP")

modelList <- list("FEP_O3_mean_city",
                  c("FEP_SO2_mean_city","FEP_O3_mean_city"),
                  c("FEP_CO_mean_city","FEP_SO2_mean_city","FEP_O3_mean_city"),  
                  c("FEP_CO_mean_city","FEP_PM2_5_mean_city","FEP_O3_mean_city","FEP_SO2_mean_city"), 
                  c("FEP_CO_mean_city","FEP_PM2_5_mean_city","FEP_O3_mean_city","FEP_SO2_mean_city","FEP_NO2_mean_city"))

c10 <- lapply(modelList, FUN = function(poll) lapply(labels, FUN = function(vd) twoSLS_on_selected(data, vd, m, poll, baseline[[1]])))
ci <- Reduce(rbind,lapply(2:length(modelList), FUN = function(i) ci_df(c10[[i]], modelList[[i]]))) 
ci$age <- factor(ci$age, levels = c('nn','29j-1','0-4','5-14','15-59','60-64','65-69','70-74','75-79','80-PP','all')) 
ggplot(ci[ci$polluant!= "NO2" & ci$polluant!="PM2_5" & ci$age!= "nn" & ci$age!="29j-1",],aes(age,estimate, ymin = cim, ymax = ciM, col = `N Pollutants`))+geom_pointrange(position = position_dodge(width = 0.3)) + facet_wrap(~polluant, scales = "free_y") + scale_color_manual(values = c("blue","#3A4F8F","grey","black")) + theme_bw() + ylab('Estimate') + xlab("Age Group")
ggplot(ci[ci$polluant!= "NO2" & ci$polluant!="PM2_5" & ci$age %in% c("nn","29j-1","0-4","all"),],aes(age,estimate, ymin = cim, ymax = ciM, col = `N Pollutants`))+geom_pointrange(position = position_dodge(width = 0.3)) + facet_wrap(~polluant, scales = "free_y") + scale_color_manual(values = c("blue","#3A4F8F","grey","black")) + theme_bw() + ylab('Estimate') + xlab("Age Group")


labels <- c( "FEP_rate_decAGE_all" , "FEP_rate_decAGE_m4ans" , "FEP_rate_decAGE_m65ans" ,"FEP_rate_decAGE_6574","FEP_rate_decAGE_o75"  )

modelList <- list("FEP_PM2_5_mean_city",
                  c("FEP_PM2_5_mean_city","FEP_SO2_mean_city"),
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_SO2_mean_city"),  
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_SO2_mean_city","FEP_O3_mean_city"), 
                  c("FEP_PM2_5_mean_city","FEP_CO_mean_city","FEP_SO2_mean_city","FEP_O3_mean_city","FEP_NO2_mean_city"))


dec <- lapply(modelList, FUN = function(poll) lapply(labels, FUN = function(vd) twoSLS_on_selected(data, vd, m, poll, baseline[[1]])))
ci <- Reduce(rbind,lapply(2:length(modelList), FUN = function(i) ci_df(dec[[i]], modelList[[i]]))) 
ci$age <- ifelse(ci$age == 'm4ans', '0-4',ifelse(ci$age == '6574', '65-74', ifelse(ci$age == 'o75','>= 75',ifelse(ci$age == 'm65ans','5-64','all')))) 
ci$age <- factor(ci$age, levels = c('0-4','5-64','65-74','>= 75','all')) 
ggplot(ci[ci$polluant!= "NO2" & ci$polluant!="O3" & ci$polluant!="CO",],aes(age,estimate, ymin = cim, ymax = ciM, col = `N Pollutants`))+geom_pointrange(position = position_dodge(width = 0.3)) + facet_wrap(~polluant, scales = "free_y") + scale_color_manual(values = c("blue","#3A4F8F","grey","black")) + theme_bw() + ylab('Estimate') + xlab("Age Group")
