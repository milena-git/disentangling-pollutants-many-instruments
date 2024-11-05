#' Residualized from location-specific seasonal fixed effects.
purge_FE <- function(data, var){
  vector <- data[,get(var)]
  r<-felm(vector ~ day+day:ville|month_year*ville,data=data)
  data[!is.na(vector),c(var):= .(r$residuals)]
  print(var)
}
