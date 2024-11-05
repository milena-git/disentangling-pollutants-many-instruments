rm(list=ls())

## =============== Other outcomes ============= ##
# Congestion
load("Input Data/data.Rdata")
# Nuitée hotelière
load("Input Data/TO_AUdate.Rdata")
u$LIBAU2010<-sapply(strsplit(u$LIBAU2010," "),FUN=function(x) toupper(x[[1]]))
v$ville<-toupper(v$ville)
data<-merge(v,u, by.x = c('ville','date'), by.y = c('LIBAU2010',"date"))
            
# Accidents 
load("Input Data/reg.Rdata")
v$ville<-toupper(v$ville)
data<-merge(data, v[,c('ville','date','nb_acc_par_mh')],by=c("ville","date"))

# Nombre total d'admissions en urgence par ville et date
load("Input Data/nombres_durgences.Rdata")
A$aire_urbaine<-sapply(strsplit(A$aire_urbaine,"-"),FUN=function(x) toupper(x[[2]]))
A$aire_urbaine<-sapply(strsplit(A$aire_urbaine," "),FUN=function(x) toupper(x[[1]]))
data<-merge(data,A, by.x = c('ville','date'), by.y = c('aire_urbaine',"date_entree"))
data$log_n_entree<-log(data$n_entree)

oo<-data

load("Data/base_IV_purged_FE.Rdata")
data<-merge(data,oo,by=c("ville","date"), all.x = T)


# EF temporels
EFT<-"ville+day+month_year+day:ville+month_year:ville"

# Weather
m<-paste0(c("neige_sol","verglas_gel_sol","neige_couvrant_tout_sol","precipitations","precipitations_2",
        "temperature","temperature_2","ff","ff_2","humidite","ensoleil_glo"), collapse = "+")

data$day_ville<-factor(paste0(data$day,data$ville))
data$month_year_ville<-factor(paste0(substr(data$date,1,7),data$ville))

data$inv_hcl_mean_1000<-data$inv_hcl_mean*1000


a1<-felm(as.formula(paste0("`pollution_index_pca_imputed_F1`~ jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
a2<-felm(as.formula(paste0("`pollution_index_pca_imputed_F1`~`mean_to`+jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
b1<-felm(as.formula(paste0("`log_n_entree`~ jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
b2<-felm(as.formula(paste0("`log_n_entree`~`mean_to`+jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
c1<-felm(as.formula(paste0("`inv_hcl_mean_1000`~ jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
c2<-felm(as.formula(paste0("`inv_hcl_mean_1000`~`mean_to`+jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
d1<-felm(as.formula(paste0("`invth_50p_mean`~ jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)
d2<-felm(as.formula(paste0("`invth_50p_mean`~`mean_to`+jf+vac+wep+",m,"|day_ville+month_year_ville|0|month_year_ville")),data=data)


stargazer(a1,a2,b1,b2,c1,c2,d1,d2,
                         omit.stat=c("ser","F","adj.rsq"),dep.var.labels.include = F,column.labels = c("Pollution","Hospital entries","IBLH","InvTH"),
                         column.separate = c(2,2,2,2),title="Proxy population present",
                         omit=c("month_year","Constant","neige","verglas","temp","precip","ff","dd","humidi","pstat","tsv","ensoleil","ville","^day","week","^an"))

