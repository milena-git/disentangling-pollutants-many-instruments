randomization_estimate <- function(data, fixedcols = c("ville","date",outcomes, pollutants, m), shufflecols = c("ville","date", selected_set), grouping = "date", random_set = unique(data$ville), random_var = "ville"){
  allcols <- unique(c(fixedcols, shufflecols, grouping))
  setDT(data)
  dataRI <- data[,..allcols]
  dataRI[, random:= sample(random_set, min(.N,length(random_set))),by = grouping]

  shuffle <- unique(dataRI[,.(random,date,ville)])
  shuffle[,.(mean(random==get(random_var)))]
  newvars <- dataRI[,.(random = get(random_var), grp = get(grouping), .SD), .SDcols = selected_set]
  names(newvars) <- c("random", grouping, selected_set)
  
  shuffle <- merge(shuffle,newvars, by = c("random", grouping))
  
  shcols <- c(shufflecols, "random")
  dataRI <- merge(shuffle[,..shcols], dataRI[,..fixedcols], by = c('date','ville'))
  dataRI[, month_year_ville:= paste0(substr(date,1,7),ville)]
  
  a <- lapply(outcomes, FUN = function(vd) twoSLS_on_selected(dataRI, vd, m, pollutants, selected_set))
  r <- rbindlist(lapply(a, FUN = function(i) {res <- data.frame(summary(i)$coefficients[grepl('mean_city', rownames(summary(i)$coefficients)),]); 
  res$pol<- pollutants; names(res)<- c('est','se','tstat','pval',"pol"); data.frame(res)}))
  r$outcome<- rep(outcomes, each = length(pollutants))
  

  return(r)
}

randomization_estimates_ville <- function(N_draw,data, a_ref, grouping = "date"){
  r <- rbindlist(lapply(1:N_draw, FUN = function(x) randomization_estimate(data, grouping = grouping )))
  s <- rbindlist(lapply(a_ref, FUN = function(i) {res <- data.frame(summary(i)$coefficients[grepl('mean_city', rownames(summary(i)$coefficients)),]); 
  res$pol<- pollutants; names(res)<- c('est_ref','se_ref','tstat_ref','pval_ref',"pol"); data.frame(res)}))
  s$outcome<- rep(outcomes, each = length(pollutants))
  r<- merge(r,s, by=c('pol','outcome'))
  return(r)
}