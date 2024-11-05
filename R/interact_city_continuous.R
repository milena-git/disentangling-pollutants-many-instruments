interact_city_continuous<-function(data,var_text){
  for(v in unique(data$ville)){
    eval(parse(text=paste0("data[,intcity_",var_text,"_",v,":=",var_text,"*(ville==v)]")))
  }
}