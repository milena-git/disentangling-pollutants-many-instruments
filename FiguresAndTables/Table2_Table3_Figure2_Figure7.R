rm(list=ls())

load("data/base_IV_purged_FE.Rdata")

# ==================================== Table with Descriptive Statistics (Table 2) ============================================ #

library(dplyr)
data <- as.data.frame(data)

var<-c( names(data)[grepl("mean_city$",names(data)) & !grepl('^FEP_',names(data))], "pollution_index_pca_imputed_F1",
        "rate_emergencies_tot",
   names(data)[grepl("^rate_adm_C_[0-9]+_Full_all",names(data)) & !grepl("rate_adm_C_15_Full_all",names(data)) ],
   names(data)[grepl("^rate_dec.*_all$",names(data)) | grepl("^n_[I|J|K]_in",names(data)) ],"precipitations","temperature","ff","humidite","ensoleil_glo"  )
var2<-c("neige_sol","verglas_gel_sol",
        "neige_couvrant_tout_sol")
desc<-data[,var]
tab<-cbind(t(desc %>% summarise_all(funs(round(quantile(.,0.1,na.rm=T),1)))),
      t(desc %>% summarise_all(funs(round(mean(.,na.rm=T),1)))),
      t(desc %>% summarise_all(funs(round(quantile(.,0.9,na.rm=T),1)))),
      t(desc %>% summarise_all(funs(sum(is.na(.))))))
desc2<-data[,var2]
tab2<-cbind(t(desc2 %>% summarise_all(funs(round(quantile(.,0.1,na.rm=T),3)))),
           t(desc2 %>% summarise_all(funs(round(mean(.,na.rm=T),3)))),
           t(desc2 %>% summarise_all(funs(round(quantile(.,0.9,na.rm=T),3)))),
           t(desc2 %>% summarise_all(funs(sum(is.na(.))))))
nobs<-dim(data)[1]
rownames(tab)<-c("PM2.5","PM10","NO2","O3","CO","SO2",
                 "Air Pollutant Index",
                 "Emergency admissions",
                 "Cardiovascular Diseases","Respiratory Diseases","Digestive Diseases",
                 "Mortality Rate",
                 "Mortality Rate from at least one Cardiovascular Diseases Cause",
                 "Mortality Rate from at least one Respiratory Diseases Cause",
                 "Mortality Rate from at least one Digestive Diseases Cause",
                 "Precipitations","Temperature","Wind Strength","Humidity","Sun Light")
stargazer(tab,digits=1)
rownames(tab2)<-c("Snow (Dummy)","Ice (Dummy)","Snow covering Ground (Dummy)")
stargazer(tab2)
stargazer(c("Total Observations","","","",nobs))

## Descriptive statistics of Main Instruments.
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
vars <- c('s_pblh_mean','inv_hcl_mean',"invth_simple_mean","invth_strength_mean",
          "vitu_20_mean","vitu_40_mean","vitu_60_mean","vitv_20_mean","vitv_40_mean","vitv_60_mean","norm_v_20_mean","norm_v_40_mean","norm_v_60_mean",
          "zfull_20_mean","zfull_40_mean","zfull_60_mean")
data$inv_hcl_mean <- 1000*data$inv_hcl_mean
desc<-data[,..vars]
tab<-cbind(t(desc %>% summarise_all(funs(round(quantile(.,0.1,na.rm=T),1)))),
           t(desc %>% summarise_all(funs(round(mean(.,na.rm=T),1)))),
           t(desc %>% summarise_all(funs(round(quantile(.,0.9,na.rm=T),1)))),
           t(desc %>% summarise_all(funs(sum(is.na(.))))))
rownames(tab)<-c("PBLH","IPBLH (x 1,000)","Thermal Inversions (# Hours during the day)","Thermal Inversion Strength (T_up-T_0)",
                 "Zonal Wind (Layer 20)","Zonal Wind (Layer 40)","Zonal Wind (Layer 60)",
                 "Meridional Wind (Layer 20)","Meridional Wind (Layer 40)","Meridional Wind (Layer 60)",
                 "Wind Strength (Layer 20)","Wind Strength (Layer 40)","Wind Strength (Layer 60)",
                 "Altitude of Layer 20","Altitude of Layer 40","Altitude of Layer 60")
stargazer(tab,digits=1)




# =================================== Table 3 ======================================== #

stargazer(round(cor(data[,c("pollution_index_pca_imputed_F1","PM2_5_mean_city","PM10_mean_city","NO2_mean_city","O3_mean_city","CO_mean_city","SO2_mean_city")],use="pairwise.complete.obs"),2), digits = 2)


# ==================================== Figures: Inverse of PBLH and Pollutant Concentrations ============================================ #

# Seasonal Fixed-Effects #
EFT<-"ville+day_ville+month_year_ville"
data$day_ville<-factor(paste0(data$day,data$ville))
data$month_year_ville<-factor(paste0(substr(data$date,1,7),data$ville))

# Ground-level weather # 
m<-"neige_sol+verglas_gel_sol+neige_couvrant_tout_sol+precipitations+I(precipitations^2)+temperature+I(temperature^2)+ff+I(ff^2)+humidite+ensoleil_glo"

# Because of 1-day lag variables: 
data <- data %>% filter(date!="2010-01-01")
# Non-missing weather and instruments
data <- data %>% filter(!is.na(invth_simple_mean) & !is.na(inv_hcl_mean) & !is.na(neige_sol) & !is.na(verglas_gel_sol) & !is.na(neige_couvrant_tout_sol) & !is.na(precipitations) & !is.na(temperature) & !is.na(ff) & !is.na(dd) & !is.na(humidite) & !is.na(pstat) & !is.na(tsv) & !is.na(ensoleil_glo) )


# Sample is pollutant-specific #
pol <- 'NO2_mean_city'
pol_title <- 'NO2'

plot_one_corr <-function(pol = "PM2_5_mean_city", pol_title = "PM2.5", pol_fig_title = "pm25_invhcl.pdf"){
   data_temp <- data[!is.na(data[,pol]),]
   
   r<-lm(as.formula(paste0(pol,"~ ",m,"+",EFT)),data=data_temp)
   s<-lm(as.formula(paste0("inv_hcl_mean~ ",m,"+",EFT)),data=data_temp)
   data_temp$res_pol<-residuals(r)
   data_temp$res_inv_hcl_mean<-residuals(s)
   
   data_temp$inv_hcl_mean_bins<-cut(data_temp$res_inv_hcl_mean,quantile(data_temp$res_inv_hcl_mean,seq(0,1,0.05)),include.lowest = T,right=T)
   
   v<-data_temp %>% group_by(inv_hcl_mean_bins) %>% summarise(mean=mean(res_pol),
                                                              mean_inv_hcl = mean(res_inv_hcl_mean),
                                                              q25=quantile(res_pol,0.25),
                                                              q75=quantile(res_pol,0.75))
   xlim <- quantile(data_temp$res_inv_hcl_mean, c(0.01,0.99))
   ylim <- quantile(data_temp$res_pol, c(0.2,0.8))
   
   pdf(pol_fig_title,width=5,height=4)
   
   print(ggplot()+ggtitle(pol_title)+geom_smooth(data=data_temp,aes(res_inv_hcl_mean,res_pol),method='lm') +
      coord_cartesian(xlim = xlim, ylim = ylim) +
      geom_point(data=v,aes( mean_inv_hcl,mean,ymin=q25,ymax=q75,group=1))+theme_bw()+ylab(paste0(pol_title," deviation from prediction"))+
      xlab("Inverse of boundary layer deviation from prediction"))
   
   dev.off()

}

plot_one_corr()
plot_one_corr('NO2_mean_city',"NO2","no2_invhcl.pdf")
plot_one_corr('O3_mean_city',"O3","o3_invhcl.pdf")
plot_one_corr('SO2_mean_city',"SO2","so2_invhcl.pdf")
plot_one_corr('CO_mean_city',"CO","co_invhcl.pdf")
plot_one_corr('PM10_mean_city',"PM10","pm10_invhcl.pdf")


plot_one_corr_uncond <-function(pol = "PM2_5_mean_city", pol_title = "PM2.5", pol_fig_title = "pm25_invhcl_uncond.pdf"){
   data_temp <- data[!is.na(data[,pol]),]
   
   eval(parse(text = paste0("data_temp$pol <- data_temp$", pol)))
   
   data_temp$inv_hcl_mean_bins<-cut(data_temp$inv_hcl_mean, quantile(data_temp$inv_hcl_mean,seq(0,1,0.05)),include.lowest = T,right=T)
   
   v<-data_temp %>% group_by(inv_hcl_mean_bins) %>% summarise(mean=mean(pol),
                                                                      mean_inv_hcl = mean(inv_hcl_mean),
                                                                      q25=quantile(pol,0.25),
                                                                      q75=quantile(pol,0.75))
   xlim <- quantile(data_temp$inv_hcl_mean, c(0.05,0.95))
   ylim <- quantile(data_temp$pol, c(0.1,0.9))
           
   pdf(pol_fig_title,width=5,height=4)
   print(ggplot()+ggtitle(pol_title)+geom_smooth(data=data_temp,aes(inv_hcl_mean,pol),method='lm') +
                    coord_cartesian(xlim = xlim, ylim = ylim) +
                    geom_point(data=v,aes( mean_inv_hcl,mean,ymin=q25,ymax=q75,group=1))+theme_bw()+ylab(paste0(pol_title))+
                    xlab("Inverse of boundary layer"))
   dev.off()

}

plot_one_corr_uncond()
plot_one_corr_uncond('NO2_mean_city',"NO2","no2_uncond_invhcl.pdf")
plot_one_corr_uncond('O3_mean_city',"O3","o3_uncond_invhcl.pdf")
plot_one_corr_uncond('SO2_mean_city',"SO2","so2_uncond_invhcl.pdf")
plot_one_corr_uncond('CO_mean_city',"CO","co_uncond_invhcl.pdf")
plot_one_corr_uncond('PM10_mean_city',"PM10","pm10_uncond_invhcl.pdf")

# ==================================== Figure: Thermal Inversion Occurence and Pollutant Concentrations ============================================ #

plot_one_corr_2 <-function(pol = "PM2_5_mean_city", pol_title = "PM2.5", pol_fig_title = "pm25_invth.pdf"){
  data_temp <- data[!is.na(data[,pol]),]
  
  r<-lm(as.formula(paste0(pol,"~ ",m,"+",EFT)),data=data_temp)
  t<-lm(as.formula(paste0("invth_simple_mean~ ",m,"+",EFT)),data=data_temp)
  data_temp$res_pol<-residuals(r)
  data_temp$res_invth_simple<-residuals(t)
  
  data_temp$inv_th_mean_bins<-cut(data_temp$res_invth_simple,quantile(data_temp$res_invth_simple,seq(0,1,0.05)),include.lowest = T,right=T)
  
  v<-data_temp %>% group_by(inv_th_mean_bins) %>% summarise(mean=mean(res_pol),
                                                             mean_inv_th = mean(res_invth_simple),
                                                             q25=quantile(res_pol,0.25),
                                                             q75=quantile(res_pol,0.75))
  xlim <- quantile(data_temp$res_invth_simple, c(0.01,0.99))
  ylim <- quantile(data_temp$res_pol, c(0.2,0.8))
  
  pdf(pol_fig_title,width=5,height=4)
  
  print(ggplot()+ggtitle(pol_title)+geom_smooth(data=data_temp,aes(res_invth_simple,res_pol),method='lm') +
          coord_cartesian(xlim = xlim, ylim = ylim) +
          geom_point(data=v,aes( mean_inv_th,mean,ymin=q25,ymax=q75,group=1))+theme_bw()+ylab(paste0(pol_title," deviation from prediction"))+
          xlab("Hours with thermal inversion (deviation from prediction)"))
  
  dev.off()
  
}


plot_one_corr_2()
plot_one_corr_2('NO2_mean_city',"NO2","no2_invth.pdf")
plot_one_corr_2('O3_mean_city',"O3","o3_invth.pdf")
plot_one_corr_2('SO2_mean_city',"SO2","so2_invth.pdf")
plot_one_corr_2('CO_mean_city',"CO","co_invth.pdf")
plot_one_corr_2('PM10_mean_city',"PM10","pm10_invth.pdf")

