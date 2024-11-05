bootstrap_mono_multi_equality_test<-function(N_b,data,outcomes, m, pollutants,selected_set, cluster = "month_year_ville"){
  allcols <- c(outcomes, m, pollutants,selected_set)
  b <- data[,..allcols]
  b$random <- b[,get(cluster)]
  clusters <- unique(b[,get(cluster)])
  
  ## What sampling should we use?
  samples <- rbindlist(lapply(1:N_b, FUN = function(x) data.frame(bootstrap = x, random = sample(clusters, size = length(clusters), replace = TRUE))))
  
  datasamples <- merge(samples, b, by = "random", allow.cartesian = TRUE)
  res <- rbindlist(lapply(1:N_b, FUN = function(x) get_coefs(datasamples[bootstrap == x],outcomes, m, pollutants,selected_set)))
}


get_coefs <- function(sample,outcomes, m, pollutants,selected_set){
  
  multi <- lapply(outcomes, FUN = function(vd) { u<- coefficients(twoSLS_on_selected(sample, vd, m, pollutants,selected_set)); u[grepl('mean_city',names(u))]})
  
  mono <-lapply(outcomes, FUN = function(vd) sapply(pollutants, FUN = function(x) { u<- coefficients(twoSLS_on_selected(sample, vd, m, x,selected_set)); u[grepl('mean_city',names(u))]}))
  
  data.frame(outcomes = rep(outcomes, each = length(pollutants)), pollutant = rep(pollutants, length(outcomes)*2), estimate = c(unlist(multi),unlist(mono)), case = rep(c('multi','mono'),each = length(outcomes)*length(pollutants)) )           
}
