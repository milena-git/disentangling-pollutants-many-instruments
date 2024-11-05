MHT_corrections<-function(alpha, grouping = "", outc = c("FEP_rate_adm_C_10_Full_all","FEP_rate_adm_C_09_Full_all","FEP_rate_dec_all")){
  
  a <- lapply(outc, FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants, selected_set))
  
  if(grouping != ""){
    r <- rbindlist(lapply(a, FUN = function(i) {res <- data.frame(summary(i)$coefficients[grepl('mean_city', rownames(summary(i)$coefficients)),]); res$pol<- pollutants; names(res)<- c('est','se','tstat','pval',"pol"); data.frame(res)}))
    r$outcome<- rep(c(outc), each = length(pollutants))
    r<- r[order(get(grouping),pval)]
    r[,FWERthreshold:=alpha/(.N:1), by = .(get(grouping))][order(pval),FWEReject:= pval < FWERthreshold]
    r[,FDRthreshold:=(1:.N)*alpha/(.N), by = .(get(grouping))][order(pval),FDRReject:= pval < FDRthreshold]
  }else{
    r <- rbindlist(lapply(a, FUN = function(i) {res <- data.frame(summary(i)$coefficients[grepl('mean_city', rownames(summary(i)$coefficients)),]); res$pol<- pollutants; names(res)<- c('est','se','tstat','pval',"pol"); data.frame(res)}))
    r$outcome<- rep(c(outc), each = length(pollutants))
    r<- r[order(pval)]
    r[,FWERthreshold:=alpha/(.N:1)][order(pval),FWEReject:= pval < FWERthreshold]
    r[,FDRthreshold:=(1:.N)*alpha/(.N)][order(pval),FDRReject:= pval < FDRthreshold]  
  }
  return(r)
}

FWER_one_model <- function(res){
  coef <- summary(res)$coefficients
  p_vals <-coef[grepl('mean_city',rownames(coef)),4]
  ordered <- p_vals[order(p_vals)]
  fwer_c<- lapply(c(0.1,0.05,0.01), FUN =function(alpha) ordered[ordered < alpha/(length(ordered):1)])
}

FWER_models <- function(list_res, names_res = pollutants){
  p_vals<- unlist(lapply(list_res, FUN = function(res) summary(res)$coefficients[grepl('mean_city',rownames(summary(res)$coefficients)),4] ))
  names(p_vals) <- pollutants
  ordered <- p_vals[order(p_vals)]
  fwer_c<- lapply(c(0.1,0.05,0.01), FUN =function(alpha) ordered[ordered < alpha/(length(ordered):1)])
}