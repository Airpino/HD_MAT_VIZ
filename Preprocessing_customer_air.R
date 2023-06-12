# customer satisfaction air
library(tidyverse)
load("Customer_air.RData")
data<-train
train<-train %>% mutate(ageclass=cut(Age,breaks = c(0,18,34,50,65,100)))
train<-train %>% mutate(Distance=cut(`Flight Distance`,breaks = c(0,1000,10000)))

CLASSIFIED<-train %>% group_by(Gender,ageclass,`Customer Type`,`Type of Travel`,Class,Distance) %>%
  mutate(group_id = cur_group_id()) 

CLASSIFIED_G<-CLASSIFIED %>%  group_by(group_id) %>% tally() %>% filter(n>100) %>% ungroup()
CLASSIFIED<-CLASSIFIED %>% right_join(CLASSIFIED_G,by = c("group_id"))
CLASSIFIED<-CLASSIFIED %>% select(-c(group_id,n))
CLASSIFIED<-CLASSIFIED %>% group_by(Gender,ageclass,`Customer Type`,`Type of Travel`,Class,Distance) %>%
  mutate(group_id = cur_group_id())
CLASSIFIED_G<-CLASSIFIED %>%  group_by(group_id) %>% tally() %>% filter(n>50) %>% ungroup()
CLASSIFIED<-CLASSIFIED %>% right_join(CLASSIFIED_G,by = c("group_id"))
save(CLASSIFIED,train,file="Customer_air.RData")
#%>% 
#  mutate(count=n(group_id))
#create 