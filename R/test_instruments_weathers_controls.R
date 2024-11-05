test_instruments_weathers_controls<-function(m, instruments){
  
  selected_set<-lapply(pollutants, FUN = function(x) selected_instruments(data5, m, instruments, x, clustervar = 'month_year_ville', baseline_only = TRUE))
  selected_set <-unique(unlist(selected_set))
  out <- c( "FEP_rate_adm_C_10_Full_all", "FEP_rate_adm_C_09_Full_all", "FEP_rate_emergencies_tot",  "FEP_rate_adm_C_11_Full_all", "FEP_rate_dec_all", "FEP_n_I_in_CertD", "FEP_n_J_in_CertD"   ,"FEP_n_K_in_CertD"    )
  a <- lapply(out, FUN = function(vd) twoSLS_on_selected(data5, vd, m, pollutants,selected_set))
  return(a)
}