get_angle_cat<-function(x,y){
  theta<-acos(x/sqrt(x^2+y^2))*(y>0)+(pi+pi-acos(x/sqrt(x^2+y^2)))*(y<=0)
  cat<-ifelse(theta<pi/8 |  theta>=15/8*pi,"E",
              ifelse(theta>=pi/8 & theta<3*pi/8,"NE",
                     ifelse(theta>=3*pi/8 & theta<5*pi/8,"N",  
                            ifelse(theta>=5*pi/8 & theta<7*pi/8,"NO",  
                                   ifelse(theta>=7*pi/8 & theta<9*pi/8,"O",  
                                          ifelse(theta>=9*pi/8 & theta<11*pi/8,"SO",   
                                                 ifelse(theta>=11*pi/8 & theta<13*pi/8,"S","SE")))))))
  return(as.factor(cat))
}
