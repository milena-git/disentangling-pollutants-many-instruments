#' Run 2SLS for a set of instruments
twoSLS_on_selected <- function(data, vd, m, pollutants, selected, FE = '0',clust = 'month_year_ville'){
  felm(as.formula(paste0("`",vd,"`~", paste(m,collapse = "+"),"|",FE,"|
                             (",paste0(pollutants,collapse="|"),"~",paste(selected,collapse="+"),")|", clust)),data=data)
}