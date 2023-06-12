#create tibble Customer air
load("Customer_air.RData")
n<-max(CLASSIFIED$group_id)
VARS<-list()

  cc<-0
  for (j in 8:21){
    cc<-cc+1
    name_v<-names(CLASSIFIED)[j]
    VARS[[name_v]]<-list()
    for(i in 1:n){
        tmp<-CLASSIFIED %>% filter(group_id==1) %>% group_by_at(j) %>% tally()
    names(tmp)[1]<-"x"
    names(tmp)[2]<-"abs"
    tmp$freq<-tmp$abs/sum(tmp$abs)
    VARS[[name_v]][[i]]<-tmp
  }
  }
  
  MYTIB<-as_tibble(VARS)
  