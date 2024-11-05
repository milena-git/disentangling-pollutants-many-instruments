#' Run robust lasso selection for the given polluant under 4 estimations
#' (i) selected_robust. Robust lasso selection of both instruments and controls, and reintroduction of controls post selection
#' (ii) selected_robustC. As selected_robust, but implement clustered-lasso
#' (iii) selected_robustunp. Robust lasso selection of instruments, controls being unpenalized in the lasso problem 
#' (iv) selected_robustunpC. As selection_robustunp, but implement clustered-lasso
#' 
#' Prefered estimation: (iv)
selected_instruments <- function(data, meteo, instrument, polluant, clustervar = 'month_year_ville', baseline_only = FALSE){
  
  var_list <- c(polluant , meteo, instrument, clustervar)
  var_x <- c(meteo, instrument)
  df <- data[,..var_list]
  df <- df[complete.cases(df)]
  df$temp <- 1
  
  fs0_c <- rlasso_unpenalized_controls(y = as.matrix(df[,get(polluant)]), 
                                       x = as.matrix(df[, ..var_x]), 
                                       index_x_no_selection = 1:length(m), 
                                       clusters = as.matrix(df[,get(clustervar)]))
  selected_robustunpC <- unique(var_x[fs0_c$coefficients[-1]!=0])
  
  if(baseline_only == FALSE){
    fs <- rlasso_reintroduce_controls(y = as.matrix(df[,get(polluant)]), 
                                 x = as.matrix(df[, ..var_x]), 
                                 index_x_no_selection = 1:length(m), 
                                 clusters = NULL)
    
    fs_c <- rlasso_reintroduce_controls(y = as.matrix(df[,get(polluant)]), 
                                   x = as.matrix(df[, ..var_x]), 
                                   index_x_no_selection = 1:length(m), 
                                   clusters = as.matrix(df[,get(clustervar)]))
    
    
    fs0 <- rlasso_unpenalized_controls(y = as.matrix(df[,get(polluant)]), 
                                       x = as.matrix(df[, ..var_x]), 
                                       index_x_no_selection = 1:length(m), 
                                       clusters = NULL)
    
    
    selected_robust <- unique(c(meteo, var_x[fs$coefficients[-1]!=0]))
    selected_robustC <- unique(c(meteo,var_x[fs_c$coefficients[-1]!=0]))
    selected_robustunp <- unique(var_x[fs0$coefficients[-1]!=0])
    
    return(list( selected_robust = selected_robust, 
                 selected_robustC = selected_robustC,
                 selected_robustunp = selected_robustunp,
                 selected_robustunpC = selected_robustunpC)) 
  }else{
    return(list(selected_robustunpC = setdiff(selected_robustunpC,meteo)))
  }
  
}
