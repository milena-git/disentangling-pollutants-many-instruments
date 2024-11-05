rm(list=ls())

path <- "//abra/current/CRT_PRJ_AIR_POLLUTION/"
setwd(path)

## Normalisation
load("population/population_cl_age_AU.Rdata")
r0$ville<-r0$LIBAU2010
r0$ville[r0$ville=="Strasbourg (partie française)"]<-"Strasbourg"
r0$ville[r0$ville=="Lille (partie française)"]<-"Lille"
r0$ville[r0$ville=="Marseille - Aix-en-Provence"]<-"Marseille"
r0$cl_age<-gsub("pop","",r0$cl_age)
setDT(r0)
r0[as.numeric(substr(cl_age,1,2))<60 & as.numeric(substr(cl_age,1,2))>=15,cl_age:="15-59"]
r0$cl_age[r0$cl_age=="0004"]<-"0-4"
r0$cl_age[r0$cl_age=="0509"]<-"5-14"
r0$cl_age[r0$cl_age=="1014"]<-"5-14"
r0$cl_age[r0$cl_age=="6064"]<-"60-79"
r0$cl_age[r0$cl_age=="6569"]<-"60-79"
r0$cl_age[r0$cl_age=="7074"]<-"60-79"
r0$cl_age[r0$cl_age=="7579"]<-"60-79"
r0$cl_age[r0$cl_age=="80P"]<-"80-PP"


r00 <- r0[,.(pop =sum(pop)), by = .(ville)]
r0<-r0[,.(pop =sum(pop)), by = .(ville,cl_age)]

young <- r0[cl_age %in% c('0-4','5-14')][,.(pop = sum(pop), cl_age = "0-14"),by = .(ville)]
r0<-rbind(r0, young )


######################### -------------------- Emergencies by age & diseases -  PMSI ------------------------- ######################################
B<-read_sas("PMSI urgence/transmission ATIH/demande 2018/base_pmsi_2010_2017/base_pmsi_2010_2017.sas7bdat")
setDT(B)
B <- B[entree_urgence == 1][,.(effectif = sum(effectif)), by = .(cl_age, pathologie, date_entree, aire_urbaine)]
B[,sum(effectif)]

z<-data.frame(pathologie=as.character(unique(B$pathologie)),chap=substr(unique(B$pathologie),1,2))
detail<-sapply(strsplit(as.character(z[,1])," - "),FUN=function(x) gsub("[(-)]","",x[2]))
detail[is.na(detail)]<-"Full"
detail<-sapply(strsplit(detail," "),FUN=function(x) x[1])
detail[detail=="J21"]<-"J21-h0"
z<-data.frame(z,detail=detail)

B <- merge(B, z,by = "pathologie")
B <- B[, .(effectif = sum(effectif)) , by =.(cl_age, chap, date_entree, aire_urbaine)]

# Classe d'Ã¢ge
B$date<-as.Date(dmy(B$date_entree))

B$ville<-""
B$ville[B$aire_urbaine=="001-Paris"]<-"Paris"
B$ville[B$aire_urbaine=="006-Bordeaux"]<-"Bordeaux"
B$ville[B$aire_urbaine=="011-Rennes"]<-"Rennes"
B$ville[B$aire_urbaine=="008-Nantes"]<-"Nantes"
B$ville[B$aire_urbaine=="009-Strasbourg (partie française)"]<-"Strasbourg"
B$ville[B$aire_urbaine=="004-Toulouse"]<-"Toulouse"
B$ville[B$aire_urbaine=="003-Marseille - Aix-en-Provence"]<-"Marseille"
B$ville[B$aire_urbaine=="007-Nice"]<-"Nice"
B$ville[B$aire_urbaine=="002-Lyon"]<-"Lyon"
B$ville[B$aire_urbaine=="005-Lille (partie française)"]<-"Lille"

table(B$ville)


## Décompte des urgences par grand chapitres
A<-B[,.(effectif = sum(effectif)) , by =.(chap,date, ville)]
save(A,file="PMSI urgence/data/nombres_durgences_par_chapitre.Rdata")


B[cl_age %in% c("1 an","3 ans","4 ans", "nouveau-nés","29 jours à m","2 ans"), cl_age_new:="0-4"]
B[cl_age %in% c( "5 à 9 ans" , "10 à 14 ans"), cl_age_new:="5-14"]
B[cl_age %in% c("15 à 19 ans","20 à 59 ans"), cl_age_new:="15-59"]
B[cl_age %in% c("60 à 64 ans",
                "65 à 69 ans",
                "70 à 74 ans",
                "75 à 79 ans"), cl_age_new:="60-79"]
B[cl_age %in% c( "80 à 84 ans", "85 à 89 ans","90 à 94 ans","95 à 99 ans", "100 ans et +"), cl_age_new:="80-PP"]

table(B$cl_age, B$cl_age_new)

## Chapter Selections
B <- B[chap %in% c("09","10","11")]


# /!\ Erreur dans les agrégations dans les codes précédents :
# - (100 ans et +) dans le code précédent
# x 2 (nn & 29j-1 ans) 
# X 2 les chapitres déjà présents dans la base: 11
# Comparaisons avant/aprés des effectifs.
#   09 763334        761 722 
#   10 751029        868 428 
#   11 681733	      1 397 528 

B[,.(sum(effectif*(cl_age=="100 ans et +"))/sum(effectif))] # 0.003% des admissions du champ.
B[,.(sum(effectif*(cl_age=="100 ans et +"))), by = chap]
B[,.(sum(effectif*(cl_age=="nouveaux-nés" | cl_age == "29 jours à m"))), by = chap]
B[,.(sum(effectif)), by = .(cl_age_new, chap)] 
B[,.(sum(effectif)), by = .(chap)] 


B$cl_age <- B$cl_age_new
B$cl_age_new <- NULL

B<- B[,.(effectif = sum(effectif)), by = .(cl_age, chap, date, ville)]
Byoung <- B[cl_age %in% c('0-4','5-14'),.(effectif=sum(effectif), cl_age = "0-14")  ,by =.(ville,date, chap)]
B <- rbind(B, Byoung)
B[,sum(effectif) ,by=cl_age] # Attention une classe d'âge redondante ici.

## Panel ##
B<- B[,.(effectif = sum(effectif)), by = .(cl_age, chap, date, ville)]

sk<-expand.grid(chap=unique(B$chap),date=unique(B$date),ville=unique(B$ville),cl_age=unique(B$cl_age))

df <- merge(sk, B, by= c("date","ville","chap","cl_age"),all.x = TRUE, all.y = FALSE) 
setDT(df)

df[is.na(effectif), effectif:=0]
df[, detail:="Full"]


df <- merge(df,r0, by = c('cl_age','ville'))
# Attention à enlever la classe d'âge redondante dans le calcul du total.
tot <- merge(df[cl_age != "0-14",.(cl_age= "all", effectif = sum(effectif),detail="Full"),by=c('date','ville','chap')],r00, by="ville")

B[cl_age!="0-14",sum(effectif)] == tot[,sum(effectif)]

df <- rbind(df, tot)

df$var<-paste("C",df$chap,"Full",df$cl_age,sep="_")

df$rate_adm<- df$effectif/df$pop*100000

df<-df[,c("cl_age","chap","detail"):=NULL]


df<-dcast(df,date+ville~var,value.var=c("rate_adm","pop","effectif"))

save(df,file="PMSI urgence/data/pmsi_panel_chap.Rdata")

######################### -------------------- Death - Sources: Insee (by Age), Inserm (By Cause) ------------------------- ######################################

### By Age, Source : Insee
load("Z:/Deces/deces20102015_AU/base_deces.Rdata")
setDT(ag)

ag[order(GROUPAGE),.(sum(n_dec)),by=GROUPAGE]

ag0 <- ag[,.(n_deces=sum(n_dec)) , by = c('LIBAU2010','date')]


ag[,AGE:=ifelse(as.numeric(GROUPAGE)<14,"0-14",
                  ifelse(as.numeric(GROUPAGE)>79,"80-PP",
                         ifelse(as.numeric(GROUPAGE)>14 & as.numeric(GROUPAGE)<60, "15-59",
                                ifelse(as.numeric(GROUPAGE)>59 & as.numeric(GROUPAGE)<80, "60-79","80-PP"))))]
table(ag$AGE, ag$GROUPAGE)


ag <- ag[,.(n_deces=sum(n_dec)) , by = c('LIBAU2010','date','AGE')][,.(n_deces, LIBAU2010, date, AGE = AGE)]

sum(ag$n_deces) == sum(ag0$n_deces)  #970 439

ag0$AGE<-"all"
ag<-rbind(ag0,ag)
ag <- ag[, .(ville = LIBAU2010, date, cl_age = AGE, n_deces)]

#### ----- Death By Causes ---- ####
path <-'Z:/Deces_par_cause/Fichiers sources CSV/'
file1 <- fread(paste0(path, 'airpolluants_1015.csv'))
file2 <- fread(paste0(path, 'airpolluants_2015_cause_mult.csv'))

dim(file1) # 940 024 décès, 3.5 causes en moyennes .  # 30 000 décès de plus dans la source Insee (~ 3%) (enregistrement dans la commune du dÃ©cÃ¨s). Est-ce que le lieu d'enregistrement (lieu du dÃ©cÃ¨s vs rÃ©sidence?)
file2[,.(n_cause=.N),by=id][,.(mean(n_cause))]
look <- merge(file2,file1[,.(id, CauseInitialeCode)],by="id")

# Outcomes:
# (1) Cause initiale I ou J ou K (cause initiale = tableau de mortalité statistique officielle)
# (2) I, J ou K apparaÃ®t au moins une fois dans une cause figurant sur le certificat de décès
file2[,chapitre:=substr(CodeCIM,1,1)]
file2 <- file2[chapitre %in% c("I","J","K")]
file2 <- file2[,.(I_in_CertD=any(chapitre == "I"), J_in_CertD=any(chapitre == "J"), K_in_CertD=any(chapitre == "K")),by=id]

file1 <- merge(file1,file2, by = 'id', all.x = T, all.y = F)

for (i in names(file1)[grepl('CertD$',names(file1))]){ file1[is.na(get(i)), (i):= FALSE] }

file1[, date:= as.Date(DateDecesRetenu)]
file1[,chapitre:=substr(CauseInitialeCode,1,1)]

table(file1$chapitre, file1$I_in_CertD)
# I     52 219850
# J  38315  22105
# K  24203  13182
table(file1$chapitre, file1$J_in_CertD)
# I 180875  39027
# J      2  60418
# K  30768   6617
table(file1$chapitre, file1$K_in_CertD)
# I 212007   7895
# J  58957   1463
# K     32  37353
table(file1$J_in_CertD, file1$K_in_CertD)
table(file1$I_in_CertD, file1$K_in_CertD)

deathCauses <- file1[,.(n_deces_inserm =.N, 
                        n_I = sum(chapitre == 'I'), n_J = sum(chapitre == 'J'), n_K = sum(chapitre == 'K'), 
                        n_I_in_CertD = sum(I_in_CertD),n_J_in_CertD = sum(J_in_CertD), n_K_in_CertD = sum(K_in_CertD),
                        n_IJ_in_CertD = sum(I_in_CertD & J_in_CertD),
                        n_K_not_IJ_in_CertD = sum(K_in_CertD & !I_in_CertD & !J_in_CertD)),by=.(date, LIBAU2010)]
deathCauses[, lapply(.SD,sum), .SDcols = names(deathCauses)[grepl('^n_',names(deathCauses))]]
#    n_deces_inserm    n_I   n_J   n_K n_I_in_CertD n_J_in_CertD n_K_in_CertD n_IJ_in_CertD n_K_not_IJ_in_CertD
# 1:         940024 219902 60420 37385       401858       216753        90365         97968               48917
deathCauses$LIBAU2010<-toupper(deathCauses$LIBAU2010)
deathCauses$LIBAU2010[deathCauses$LIBAU2010=="LILLE (PARTIE FRANÇAISE)"]<-"LILLE"
deathCauses$LIBAU2010[deathCauses$LIBAU2010=="MARSEILLE - AIX-EN-PROVENCE"]<-"MARSEILLE"
deathCauses$LIBAU2010[deathCauses$LIBAU2010=="STRASBOURG (PARTIE FRANÇAISE)"]<-"STRASBOURG"
deathCauses[, ville:=LIBAU2010]
deathCauses$LIBAU2010<-NULL

# 
r0[, ville:= toupper(ville)]
r00[, ville:= toupper(ville)]
ag<- merge(ag, rbind(r0, r00[,.(ville, cl_age = "all",pop)]),by=c("ville","cl_age"), all.x = TRUE)
ag[,n_deces_sur_pop_annuelle:=n_deces/pop*100000]
deathCauses <-merge(deathCauses, r00[,.(ville,pop)], by = "ville") 

for(i in names(deathCauses)[grepl('^n_', names(deathCauses))]){ deathCauses[,(i):= get(i)/pop*100000]}


## Panel  ville * date
totinsee <- ag[cl_age=='all',.(ville,date, n_deces)]
ag <- dcast(data = ag, ville + date ~ cl_age, value.var = 'n_deces_sur_pop_annuelle')
names(ag)[-which(names(ag) %in% c("ville","date"))]<-paste0("rate_dec_",names(ag)[-which(names(ag) %in% c("ville","date"))])
for (i in names(ag)[grepl('^rate_dec',names(ag))]){ ag[is.na(get(i)), (i):= 0] }

death <- merge(ag, deathCauses, by = c("ville","date"))
death <- merge(death, totinsee, by = c("ville","date") )
cor(death$rate_dec_all,death$n_deces_inserm) # 0.94

### Total Emergencies ### 
load("PMSI urgence/data/nombres_durgences_par_chapitre.Rdata")
setDT(A)
r0[, ville:= toupper(ville)]
r00[, ville:= toupper(ville)]
A[, ville:=toupper(ville)]
A<- A[,.(n_entree_tot=sum(effectif), n_entree_9_10=sum(effectif*(as.numeric(chap) %in% c("9","10"))),n_entree_not_9_10=sum(effectif*(!as.numeric(chap) %in% c("9","10")))), by=.(ville,date)]
A<-merge(A[,.(ville,date,n_entree_tot, n_entree_9_10, n_entree_not_9_10)],r00, by ="ville")
A[, rate_emergencies_tot:= n_entree_tot/pop*100000]
A[, rate_emergencies_9_10:= n_entree_9_10/pop*100000]
A[, rate_emergencies_not_9_10:= n_entree_not_9_10/pop*100000]

A$ville <- toupper(A$ville)


### Merge All Health Outcomes ### 
load('PMSI urgence/data/pmsi_panel_chap.Rdata')
df$ville <- toupper(df$ville) 

sante <- merge(df, death, by=c("ville","date")) 
sante <- merge(sante, A[,.(ville,date, rate_emergencies_tot,rate_emergencies_9_10,rate_emergencies_not_9_10)], by=c("ville","date")) 

save(sante,file="Disentangling/data/base_sante.Rdata")

