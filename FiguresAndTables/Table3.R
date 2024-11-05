
m<-"neige_sol+precipitations+I(precipitations^2)+temperature+I(temperature^2)+ff+I(ff^2)+humidite+I(humidite^2)+ensoleil_glo"

data$day_ville<-factor(paste0(data$day,data$ville))
data$month_year_ville<-factor(paste0(substr(data$date,1,7),data$ville))

data$inv_hcl_mean_1000<-data$inv_hcl_mean*1000

f1<-felm(as.formula(paste0("`pollution_index_pca_imputed_F1`~jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
fs1<-felm(as.formula(paste0("`pollution_index_pca_imputed_F1`~jf+vac+wep+`mean_to`+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
f2<-felm(as.formula(paste0("rate_emergencies_tot~jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
fs2<-felm(as.formula(paste0("rate_emergencies_tot~jf+vac+wep+`mean_to`+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
f3<-felm(as.formula(paste0("`inv_hcl_mean_1000`~jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
fs3<-felm(as.formula(paste0("`inv_hcl_mean_1000`~`mean_to`+jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
f4<-felm(as.formula(paste0("`invth_50p_mean`~jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
fs4<-felm(as.formula(paste0("`invth_50p_mean`~`mean_to`+jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)

stargazer(f1,fs1,f2,fs2,f3, fs3,f4, fs4, type = "text",
          omit.stat=c("ser","F","adj.rsq"),dep.var.labels.include = F,
          column.labels = c("Air Pollution","All Hospiral Emergencies","Inverse of PBL Height","Thermal Inversions"),
          column.separate = c(2,2,2,2),
          omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))
