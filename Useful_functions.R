#extract kde from distributional data
Extr_dens<-function(distr, poi=200){
  vec<-HistDAWass:::COMP_MQ(distr,c(0:poi)/poi)
  mm<-density(vec)
}