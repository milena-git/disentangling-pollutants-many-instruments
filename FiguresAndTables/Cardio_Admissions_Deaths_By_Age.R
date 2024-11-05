rm(list=ls())

load("data/base_IV_purged_FE.Rdata")

names(data)[grepl('rate_dec',names(data))]

tab1 <- data[,lapply(.SD, FUN = function(x) sum(x)/6), .SDcols = c("rate_adm_C_09_Full_0-14","rate_adm_C_09_Full_15-59",
                                                                  "rate_adm_C_09_Full_60-79","rate_adm_C_09_Full_80-PP")]
tab2 <- data[,lapply(.SD, FUN = function(x) sum(x)/6), .SDcols = c("rate_dec_0-14","rate_dec_15-59",
                                                                  "rate_dec_60-79","rate_dec_80-PP")]

year <- data.frame(age= factor(rep(c("0-14","15-59","60-79",">=80"),2), levels = c("0-14","15-59","60-79",">=80")), events = unlist(c(tab1,tab2)), outcome = rep(c("Cardiovascular Emergency","Deaths (All Causes)"),each = 4))

ggplot(year,aes(age, events, col = outcome))+geom_point()+theme_bw() + ylab("Yearly Events")


tab3 <- data[,lapply(.SD, FUN = function(x) sum(x)/6), .SDcols = c("rate_adm_C_10_Full_0-14","rate_adm_C_10_Full_15-59",
                                                                   "rate_adm_C_10_Full_60-79","rate_adm_C_10_Full_80-PP")]
