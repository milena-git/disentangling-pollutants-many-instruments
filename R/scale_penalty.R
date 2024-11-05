#' From the penalty of a robust-clustered lasso, scale it to select less/more instruments.

scale_penalty <- function(data, meteo, instrument, polluant, clustervar = 'month_year_ville', grid = seq(2, 10, by = 4)){
  
  var_list <- c(polluant , meteo, instrument, clustervar)
  var_x <- c(meteo, instrument)
  df <- data[,..var_list]
  df <- df[complete.cases(df)]
  df$temp <- 1
  
  fs_c <- rlasso_unpenalized_controls(y = as.matrix(df[,get(polluant)]), 
                                      x = as.matrix(df[, ..var_x]), 
                                      index_x_no_selection = 1:length(m), 
                                      clusters = as.matrix(df[,get(clustervar)]))
  
  ## Around the robust penalty ##
  RPenC <- lapply(grid, FUN = function(scale) LassoShooting.fit(x = as.matrix(df[, ..var_x]), y = as.matrix(df[,get(polluant)]), lambda = fs_c$lambda*scale))
  
  selected_around_robustC <- lapply(RPenC, FUN = function(res) setdiff(var_x[res$coefficients!=0],meteo))
  
  return(list(selected_around_robustC = selected_around_robustC))
  
}
