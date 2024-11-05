############### ========================= Pollutants and weather ================================== ########################

p<-read_dta("Input Data/bases_aires_urbaines_pollution_horaire_toutes_villes.dta")

setDT(p)

pollutant_list <-  c("PM2_5_mean_city","PM10_mean_city","NO2_mean_city","O3_mean_city","CO_mean_city","SO2_mean_city")

p <- p[a_supprimer_passage_h_hiver == 0][,.(local_time = ymd_hms(local_time),ville = toupper(city),PM2_5_mean_city,PM10_mean_city,NO2_mean_city,O3_mean_city,CO_mean_city,SO2_mean_city)]

p[, paste0(pollutant_list,'_scaled') := lapply(.SD, FUN = scale), .SDcols = pollutant_list]

scaled_pollutant_list <- paste0(pollutant_list,'_scaled')

## PCA ## 

pollutants<-apply(p[, ..scaled_pollutant_list],MARGIN=2,FUN=function(x) as.numeric(unlist(x)))

u1<-imputePCA(pollutants,ncp=2,threshold=10^-6,maxiter=1000)
pca<-prcomp(u1$completeObs)

p$PM2_5_pca_imputed<-u1$completeObs[,1]
p$PM10_pca_imputed<-u1$completeObs[,2]
p$NO2_pca_imputed<-u1$completeObs[,3]
p$O3_pca_imputed<-u1$completeObs[,4]
p$CO_pca_imputed<-u1$completeObs[,5]
p$SO2_pca_imputed<-u1$completeObs[,6]


p$pollution_index_pca_imputed_F1<-  (u1$completeObs[,1]*pca$rotation[1,1]+u1$completeObs[,2]*pca$rotation[2,1]+u1$completeObs[,3]*pca$rotation[2,1]+u1$completeObs[,4]*pca$rotation[4,1]+u1$completeObs[,5]*pca$rotation[5,1]+u1$completeObs[,6]*pca$rotation[6,1])
p$pollution_index_pca_imputed_F2<-  (u1$completeObs[,1]*pca$rotation[1,2]+u1$completeObs[,2]*pca$rotation[2,2]+u1$completeObs[,3]*pca$rotation[2,2]+u1$completeObs[,4]*pca$rotation[4,2]+u1$completeObs[,5]*pca$rotation[5,2]+u1$completeObs[,6]*pca$rotation[6,2])

### Table de corrélations des polluants et de l'indice construit ###

cor(p[,c("pollution_index_pca_imputed_F1","PM10_mean_city","PM2_5_mean_city","NO2_mean_city","CO_mean_city","SO2_mean_city","O3_mean_city")],use="pairwise.complete.obs")

# Base Météo France

m<-read_dta("Input Data/weather_surface_all.dta")

setDT(m)
m <- m[a_supprimer_passage_h_hiver == 0][,
                                         .(local_time = ymd_hms(local_time),
                                           ville = toupper(city),
                                           precipitations,temperature,ff,dd,humidite,pstat,tsv,ensoleil_glo,
                                           ensoleil_uv,nebulosite,visibilite,hauteur_neige_totale,duree_gel_horaire_minutes,etat_ciel_present,
                                           etat_ciel_passe1,etat_ciel_passe2,etat_sol_sans_neige,etat_sol_avec_neige,neige_sol,
                                           neige_couvrant_tout_sol,verglas_gel_sol)]

########################
### Base journalière ###
########################

b <- merge(m, p, by = c('ville','local_time'))

setDT(b)
b[,date:=as.character(as.Date(ymd_hms(local_time)))]

b0 <- b

b1 <- b[,lapply(.SD,mean),by=c("ville","date"),.SD=names(b)[!grepl("(^neige_)|(etat)|(ville)|(date)|(verglas)",names(b))]]
b2 <- b[,lapply(.SD,any),by=c("ville","date"),.SD=names(b)[(grepl("neige",names(b)) | grepl("etat",names(b)) | grepl("verglas",names(b))) & !grepl("hauteur",names(b))]]
b3 <- b[,lapply(.SD,max),by=c("ville","date"),.SD=pollutant_list]
b4 <- b[,lapply(.SD,min),by=c("ville","date"),.SD=pollutant_list]
colnames(b3) <- c('ville','date', paste0(pollutant_list,'_max24h' ) )
colnames(b4) <- c('ville','date', paste0(pollutant_list,'_min24h' ) )

setkeyv(b1,c("date","ville"))
setkeyv(b2,c("date","ville"))
setkeyv(b3,c("date","ville"))
setkeyv(b4,c("date","ville"))


b <- b1[b2][b3][b4]

##
## Quelques variables météos en plus: les min et max, des bins pour les tests de robustesse. 
##

bMAX<-b0[, {cols=lapply(.SD,max) ;names(cols) = paste0(names(cols),'_max'); cols},by=c("ville","date"), .SDcols = c("temperature","precipitations","humidite","ff")]
bMIN<-b0[, {cols=lapply(.SD,min) ;names(cols) = paste0(names(cols),'_min'); cols},by=c("ville","date"), .SDcols = c("temperature","precipitations","humidite","ff")]
setkeyv(b,c("date","ville"))
setkeyv(bMAX,c("date","ville"))
setkeyv(bMIN,c("date","ville"))


b <- b[bMAX][bMIN]


cut_bins<-function(df, n_bins = 10, name_var, name_bin , breaks = NULL){
  #' How many hours of the day in each bin of the continuous variable
  #' 
  df_r <- copy(df)
  
  if(is.null(breaks)){ breaks <- unique(quantile( df_r[,  get(name_var)] , na.rm = T, probs =  seq(0, 1,1/n_bins)))} 
  
  nbins <- length(breaks) - 1 
  
  df_r$var <- cut(df_r[,  get(name_var)] ,breaks = breaks, include.lowest =  TRUE, right = TRUE)
  df_r <- dcast(df_r[, .(.N) , by = c('ville','date', 'var')], ville + date ~ var, value.var = 'N', sep = '_')
  df_r$`NA` <- NULL
  colnames(df_r) <- c('ville', 'date', paste0(name_var,"_",name_bin,"_",seq(1,nbins)))
  
  df_r <- df_r[, lapply(.SD, FUN = function(x) ifelse(is.na(x),0,x))] 
  setDT(df_r,key = c('ville','date'))
  
  
}

setDT(b, key = c('ville','date'))

## Alternative transformation of temperature: decile or 5-degrees bins.
bt<- cut_bins(b0, n_bins = 10, 'temperature', 'bin_decile')
bt1<- cut_bins(b0, n_bins = 11, 'temperature', 'bin')
bt2<- cut_bins(b0, n_bins = 4, 'temperature', 'bin_set3')
bt3<- cut_bins(b0, n_bins = 7, 'temperature', 'bin_breaks', breaks = c(-Inf, 0,5,10,15,20,25, Inf))

b <- b[bt][bt1][bt2][bt3]

### Alternative transformation of humidity, wind strength
bh <- cut_bins(b0, n_bins = 4, 'humidite', 'bin_set3')
bff <- cut_bins(b0, n_bins = 4, 'ff', 'bin_set3')

b <- b[bh][bff]

### Add wind-direction dummies
bdd <- cut_bins(b0, n_bins = 6, 'dd','binwd',breaks = c(0,60,120,180,240,300,360))
b <- b[bdd]

## Precip  ##
quantile(b$precipitations,seq(0.5,1,0.1),na.rm=T)

bp<- cut_bins(b0, n_bins = 3, 'precipitations','bin', breaks = c(0,0.025,0.258,7))

b <- b[bp]

## Ensoleil ##
be<- cut_bins(b0, n_bins = 3, 'ensoleil_glo', 'bin')

b <- b[be]

save(b,file="Data/base_meteo_france_polluants.Rdata")


########################================== Transform LMDZ output ===========================############################


# Model Ouputs (hourly). 
load("Input Data/ncdf_var_temporelles_toutes_villes_complete_final.Rdata") # 41 variables
df0 <- df
load("Input Data/ncdf_var_time_pression_toutes_villes_complete_final.Rdata") # 9 variables * 79 niveaux de pression

setDT(df)
setDT(df0)

# Axe vertical des niveaux de pression.
pres <- df[, .(layer = .N) , by = presniv][, layer:=.I]
save(pres,file="Input Data/pres_level.Rdata")

df <- merge(df, pres, by="presniv")

df$time <- df$presniv <- NULL

df<-dcast(df,local_time+ville~layer,value.var=c("ovap","temp","pres","paprs","theta","tke","vitu","vitv","zfull"))

df0 <- merge(df0, df, by = c('local_time','ville'))

rm(df)
gc()

# Add: Thermal inversions. #
# - From potential temperature.
temp_profile <- colnames(df0)[ grepl("(^theta_[0-9]$)|(^theta_1[0-8]$)", colnames(df0), fixed = FALSE)]
invtheta <- df0[,..temp_profile]

## At least 15/17, 10/17 or 5/17 mesures of theta between p2 et p18 are greater than theta at p1
temp<-apply(invtheta ,MARGIN=1, FUN=function(x) sum(x[1]<x[2:length(x)]))
invt0<-temp>14
invt1<-temp>9
invt2<-temp>5
df0$invtheta_30p<-invt2
df0$invtheta_50p<-invt1
df0$invtheta_80p<-invt0

# - From temperature.
temp_profile <- colnames(df0)[ grepl("(^temp_[0-9]$)|(^temp_1[0-8]$)", colnames(df0), fixed = FALSE)]
invth <- df0[,..temp_profile]
temp<-apply(invth, MARGIN=1, FUN=function(x) sum(x[1]<x[2:length(x)]))
invt0<-temp>14
invt1<-temp>9
invt2<-temp>5
df0$invth_30p<-invt2
df0$invth_50p<-invt1
df0$invth_80p<-invt0

# - From Jans et al. 
df0$invth_simple<-(df0$temp_9-df0$temp_1)>0
df0$invth_strength<-(df0$temp_9-df0$temp_1)

## ============= Build daily data =============== ##

df0[,date:=date(local_time)]

#Simple average, over 24h or over daytime.
c<-df0[,lapply(.SD,mean),by=c("ville","date")]
c2<-df0[hour(local_time)<21 & hour(local_time)>6]
c2<-c2[,lapply(.SD,mean),by=c("ville","date")]
colnames(c)[3:length(colnames(c))]<-paste0(colnames(c)[3:length(colnames(c))],"_mean")
colnames(c2)[3:length(colnames(c2))]<-paste0(colnames(c2)[3:length(colnames(c2))],"_daymean")

#For important instruments, keep within day variations
# - Thermal inversion, 0-4h; 4-8h, 8-12h, 12-16h, 16-20, 20-24h
e<-df0[,c("date","tranche"):=.(date(local_time),hour(local_time)  %/% 4)]
e<-e[,lapply(.SD,mean),by=c("ville","date","tranche"),.SDcols=names(e)[grepl("^invth",names(e))]]
e<-dcast(e,ville+date~tranche,value.var=names(e)[grepl("^invth",names(e))])

# - IPBLH, 0-4h; 4-8h, 8-12h, 12-16h, 16-20, 20-24h
f <- df0[,c("date","tranche"):=.(date(local_time),hour(local_time)  %/% 4)]
f<-f[,.(inv_hcl=mean(1/s_pblh)),by=c("ville","date","tranche")]
f<-dcast(f,ville+date~tranche,value.var=c("inv_hcl"))
names(f)[3:8]<-paste0("inv_hcl_",names(f)[3:8])
f<-f[,paste0("l_inv_hcl",seq(0,5,1)):=.(dplyr::lag(inv_hcl_0),dplyr::lag(inv_hcl_1),dplyr::lag(inv_hcl_2),
                                        dplyr::lag(inv_hcl_3),dplyr::lag(inv_hcl_4),dplyr::lag(inv_hcl_5)),by=c("ville")]

# - PBLH itself 0-4h; 4-8h, 8-12h, 12-16h, 16-20, 20-24h
g<-df0[,c("date","tranche"):=.(date(local_time),hour(local_time)  %/% 4)]
g<-g[,lapply(.SD,mean),by=c("ville","date","tranche"),.SDcols=names(g)[grepl("^s_pblh",names(g))]]
g<-dcast(g,ville+date~tranche,value.var="s_pblh")
colnames(g)[3:length(colnames(g))]<-paste0('s_pblh4hours_',colnames(g)[3:length(colnames(g))])


# Plus, hourly data
h <- df0[,heure := hour(local_time)]
h[, ipblh := 1/s_pblh]
h <- dcast(h, ville+date ~ heure, fun.aggregate = mean, value.var=c("ipblh","s_pblh","s_pblt","invth_50p")) # fun.aggregate because for some hour, there may be 2 observations (winter time changes)
colnames(h)[grepl('invth_50p',colnames(h))] <- paste0('hourly_',colnames(h)[grepl('invth_50p',colnames(h))])
#Merge

setkeyv(c,c("date","ville"))
setkeyv(c2,c("date","ville"))
setkeyv(e,c("date","ville"))
setkeyv(f,c("date","ville"))
setkeyv(h,c("date","ville"))
setkeyv(g,c("date","ville"))


d<-c[c2][e][f][h][g]

# ====== Some additional instruments ====== #

# Wind norm
alti<-paste0(seq(10,79),'_mean')
for(j in alti){  eval(parse(text=paste0("d$norm_v_",j,"<-sqrt(d$vitu_",j,"^2+d$vitv_",j,"^2)")))}


# IPBLH day averages
d[,c("inv_hcl_mean","inv_hcl_daymean"):=.(1/s_pblh_mean,1/s_pblh_daymean)]
d[,lag_inv_hcl_mean:=dplyr::lag(inv_hcl_mean),by=c("ville")]
d[,lag_inv_hcl_daymean:=dplyr::lag(inv_hcl_daymean),by=c("ville")]
d$precip_group<-factor(findInterval(d$precip_mean,quantile(d$precip_mean[d$precip_mean>0],seq(0,0.8,0.2))))

alti<-paste0(seq(10,79,10),'_mean') # Sample altitude levels by suffix.

# Average wind direction x norm.
for(alt in alti){
  eval(parse(text=paste0("d$angle_",alt,"<- get_angle_cat(d$vitu_",alt,",d$vitv_",alt,")")))
}

for(alt in alti){
  for(cat in c("SO","SE","NO","NE","O","S","E","N")){
    eval(parse(text=paste0("d$norm_v_angle_",cat,"_alt_",alt," <- d$norm_v_",alt,"*(d$angle_",alt,"==\"",cat,"\")")))
  }
}

# Main instruments x city
lapply(c("inv_hcl_daymean","inv_hcl_mean",
         "invth_simple_mean","invth_simple_daymean", 
         paste0("inv_hcl_",seq(0,5)), 
         paste0("l_inv_hcl",seq(0,5)), 
         paste0("s_pblh4hours_",seq(0,5)), 
         paste0("invth_50p_",seq(0,5)), 
         paste0("invtheta_50p_",seq(0,5)), 
         names(d)[grepl('norm_v_angle',names(d))]), FUN = function(var) interact_city_continuous(d, var))


d$local_time_daymean <- d$time_mean <- d$time_daymean <- d$local_time_mean <- NULL
#names(d) <- c('ville','date',paste0('INST_',names(d)[which(! names(d) %in% c('ville','date'))]))


d <- d[year(date) %in% seq(2010,2015)]
save(d, file = 'data/daily_lmd_base.Rdata')

##### Other outcomes #### 

# Vacances (& Congestion)
load("Input Data/data.Rdata")
# Nuitée hotelière
load("Input Data/TO_AUdate.Rdata")
u$ville<-sapply(strsplit(u$LIBAU2010," "),FUN=function(x) toupper(x[[1]]))
v$ville<-toupper(v$ville)
data<-merge(v,u,by=c("ville","date"))
# Accidents 
load("Input Data/reg.Rdata")
v$ville<-toupper(v$ville)
setDT(v)
oo<- merge(data , v[,.(ville,date,nb_acc_par_mh)],by=c("ville","date"))
save(oo,file="Data/other_outcomes.Rdata")


##### =========================================================================================================== ##### 

## === Merge with health outcomes, surface weather and pollutant outcomes. === # 
load('data/daily_lmd_base.Rdata')
load("data/base_sante.Rdata")
load("data/base_meteo_france_polluants.Rdata")

sante$ville <- as.character(sante$ville)
sante$date <- as.Date(sante$date)
d$ville <- as.character(d$ville)
d$date <- as.Date(d$date)
b$ville <- as.character(b$ville)
b$date <- as.Date(b$date)

setDT(sante,key = c('ville','date'))
setDT(d,key = c('ville','date'))
setDT(b,key = c('ville','date'))

data <- d[sante][b]
#data$date <- as.Date(data$date)

data <- data[date != "2010-01-01"]
data[,c('day','month_year','ville') := .( factor(wday(date)), factor(substr(as.character(date),1,7)), factor(ville))]


data[,c('precipitations_2','temperature_2','ff_2','temperature_3','precipitations_3','ff_3',"humidite_2"):=.(precipitations^2, temperature^2, ff^2,precipitations^3, temperature^3, ff^3, humidite^2)]
data0 <- copy(data) 

### =========== Variables to orthogonalize ======= ###
ovaps<-names(data)[grepl("^ovap",names(data))]
vitu<-names(data)[grepl("^vitu_",names(data))]
vitv<-names(data)[grepl("^vitv_",names(data))]
normv<-names(data)[grepl("^norm_v_[0-9]+_.+ean",names(data))]
zfull<-names(data)[grepl("^zfull_",names(data))]
other <- c("LWdnSFC_mean", "LWupSFC_mean", "SWdnSFC_mean", "SWdnSFCclr_mean","SWdnTOAclr_mean" ,   
           "SWupSFC_mean"   ,     "bils_mean"       ,    "cldh_mean"    ,       "cldl_mean"   ,        "cldm_mean"      ,     "cldt_mean"     ,            
           "flat_mean"      ,    "phis_mean"    ,       "pluc_mean"   ,        "plul_mean"      ,     "precip_mean"   ,      "prw_mean" ,          
           "psol_mean"       ,    "q2m_mean"         ,   "rh2m_mean"     ,      "s_lcl_mean"   ,       "sens_mean" ,         
           "sens_ter_mean"    ,   "soll_mean"         ,  "sols_mean"      ,     "sols0_mean"    ,      "t2m_mean"         ,   "topl_mean"       ,    "tsol_mean"  ,        
           "tsol_ter_mean"    ,   "u10m_mean"         ,  "u10m_ter_mean"  ,     "ustar_mean"    ,      "ustar_ter_mean"   ,   "v10m_mean"       ,    "v10m_ter_mean",      
           "zmax_th_mean"  )


instruments<-c(names(data)[grepl('(^invtheta_)|(^invth_simple)|(^invth_strength)|(^inv_hcl_[0-5]$)|(^l_inv_hcl[0-5]$)', names(data))],
               names(data)[grep("intcity_invth_50p_[0-5]_[A-Z]",names(data))],
               names(data)[grep("intcity_inv_hcl_[0-5]_[A-Z]",names(data))],
               names(data)[grep("intcity_l_inv_hcl[0-5]_[A-Z]",names(data))],
               names(data)[grep("intcity_s_pblh4hours_[0-5]_[A-Z]",names(data))],
               names(data)[grepl('(^ipblh_)|(^s_pblh)|(^s_pblt)|(^hourly_invth)',names(data))],
               "inv_hcl_daymean","inv_hcl_mean","lag_inv_hcl_daymean","lag_inv_hcl_mean",
               "evap_mean","lat_ter_mean",
               ovaps,
               vitu,
               vitv,
               zfull,
               other,
               normv)


weather <- c("neige_sol","verglas_gel_sol","neige_couvrant_tout_sol","precipitations","precipitations_2","precipitations_3","precipitations_min","precipitations_max",
             "temperature","temperature_2","temperature_3","temperature_min","temperature_max","ff","ff_2","ff_3","ff_min","ff_max","humidite","ensoleil_glo", "pstat","tsv",
             names(data)[grepl('bin',names(data))])

health <- c(names(data)[grepl('rate_adm_C_10_Full',names(data))],
            names(data)[grepl('rate_adm_C_09_Full',names(data))],
            names(data)[grepl('rate_dec_',names(data))],
            'rate_adm_C_11_Full_all',
            names(data)[grepl('^n_[(I|J|K)]',names(data))],
            "rate_emergencies_tot",
            "rate_emergencies_9_10",
            "rate_emergencies_not_9_10"
            )

pollutants <- names(data)[grepl('mean_city$',names(data)) |grepl('_max24h$',names(data)) |grepl('_min24h$',names(data))]

lapply(weather, FUN = function(x) purge_FE(data, x))

lapply(health, FUN = function(x) purge_FE(data, x))

lapply(pollutants, FUN = function(x) purge_FE(data, x))

lapply(instruments, FUN = function(x) purge_FE(data, x))


var_list <- c('date', 'ville',instruments, health, pollutants, weather)
data <- data[, ..var_list]
names(data)[!grepl('(^date$)|(^ville$)', names(data))] <- paste0('FEP_',names(data)[!grepl('(^date$)|(^ville$)', names(data))])

data <- merge(data, data0, by = c('ville','date'))

## Other Outcomes
load("Data/other_outcomes.Rdata")
oo$day <- oo$LIBAU2010 <- NULL
data <- merge(data,oo, by = c('ville','date'))

save(data,file="data/base_IV_purged_FE.Rdata")
