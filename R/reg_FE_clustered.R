#' Run a regression with the chosen set of fixed effect and clustering
reg_FE_clustered <- function(data, x2, y, x, FE = "day*ville+month_year_ville", cluster = "month_year_ville"){
  felm(as.formula(paste0("`",y,"`~",paste(x, collapse = "+" ), "+", 
                         paste(x2,collapse = "+"),"|",FE,"|0|",cluster)),data=data)
}

