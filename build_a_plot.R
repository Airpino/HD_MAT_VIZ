# arrange plot
#preparo la matrice
# df4<-prep_for_plot_freq(B,v=sels,
#                         n = 32,#kernel density parameter
#                         order_of_rows = order_of_rows,
#                         order_of_cols = order_of_cols)


# df4 str
# 'data.frame':	19968 obs. of  7 variables:
#   $ x      : num  0.789 3.114 5.438 7.763 10.088 ...
# $ y      : num  0 0 0 0.000849 0.012701 ...
# $ ID     : Factor w/ 78 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ ID_VAR : Factor w/ 4 levels "Temperature.C",..: 3 3 3 3 3 3 3 3 3 3 ...
# $ ID_name: Factor w/ 78 levels "I28","I40","I82",..: 21 21 21 21 21 21 21 21 21 21 ...
# $ NX     : num  0 0.0159 0.0317 0.0476 0.0635 ...
# $ VAR    : Factor w/ 5 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
library(tidyverse)



nrows<-length(levels(df4$ID_name))
vars<-length(levels(df4$ID_VAR))
### set layout
per_left<-20 #perc for left panel header
per_right<-80 #perc for right panel header
per_up <-10  #perc for top panel header
per_bott<-90  #perc for bottom panel header


x_ini<-x_fin<-y_ini<-y_fin<-val<-numeric()
for (i in 1:nrows){
  for (j in 1:vars){
    tmp<-df4 %>% 
      filter(ID_name==levels(ID_name)[i],
             ID_VAR==levels(ID_VAR)[j])
    if (nrow(tmp)==0) browser()
    tmp_v<-unlist(tmp %>% 
                    select(y))
    x_ini<-c(x_ini,unlist(tmp %>% 
                            select(NX))+j)
   # x_fin<-c(x_fin,unlist(tmp %>% 
    #                       select(NX))+j)
    y_ini<-c(y_ini,rep(i,length(tmp_v)))
    #y_fin<-c(y_fin,c(1:length(tmp_v))/length(tmp_v)+i)
    val<-c(val,tmp_v)
    
  }
}
new_df<-data.frame(x_ini=x_ini,y_ini=y_ini,value=val)
new_df<-new_df %>% 
  mutate(val=if_else(val==0,NA_real_,val),
         x_ini=(x_ini-min(x_ini))/diff(range(x_ini))*per_right+per_left,
         y_ini=(y_ini-min(y_ini))/diff(range(y_ini))*per_bott)


ppp<-ggplot(new_df,aes(x=x_ini,y=y_ini,fill=val))+
  geom_raster()+scale_fill_gradient(low="yellow",high = "red",na.value = "white")
# montiamo il dendrogramma sopra
segs_of_vars<-ddata_v$segments
rx<-range(c(unlist(segs_of_vars$x),unlist(segs_of_vars$xend)))
rx[1]<-rx[1]-0.5
rx[2]<-rx[2]+0.5
ry<-range(c(unlist(segs_of_vars$y),unlist(segs_of_vars$yend)))
segs_of_vars<-segs_of_vars%>% mutate(x=(x-rx[1])/diff(rx)*per_right+per_left,
                                          xend=(xend-rx[1])/diff(rx)*per_right+per_left,
                                          y=(y-min(ry[1]))/diff(ry)*per_up+per_bott+6,
                                     yend=(yend-min(ry[1]))/diff(ry)*per_up+per_bott+6)

ppp<-ppp +geom_segment(data=segs_of_vars,inherit.aes = F,aes(x=x,xend=xend,y=y,yend=yend),color="black")
#montiamo il dendrogramma a sinistra
segs_of_inds<-ddata$segments
rx<-range(c(unlist(segs_of_inds$x),unlist(segs_of_inds$xend)))
rx[1]<-rx[1]-0.5
rx[2]<-rx[2]+0.5
### aggiusta qui!!!!
segs_of_inds2<-segs_of_inds%>% mutate(x=(x-min(rx[1]))/diff(rx)*per_left,
                                     xend=(x-min(rx[1]))/diff(rx)*per_left,
                                     y=(y-0.5)/(max(y)+1)*per_bott,
                                     yend=(yend-0.5)/(max(yend)+1)*per_bott)

ppp<-ppp +geom_segment(data=segs_of_inds,inherit.aes = F,
                  aes(y=(x-rx[1])/(diff(rx))*per_bott,
                      yend=(xend-rx[1])/(diff(rx))*per_bott,
                      x=(1-(y+0.5)/(max(y)+1))*per_left,
                      xend=(1-(yend+0.5)/(max(yend)+1))*per_left))+
  theme_void()
## attacchiamo le etichette righe e colonne

ID_n<-data.frame(names=levels(df4$ID_name))
ID_v<-data.frame(names=levels(df4$ID_VAR))

# ggplot(ID_n,aes(x=0,y=c(1:nrow(ID_n))))+
#   geom_text_repel(aes(label=names),size=2)

ppp<-ppp+geom_text(data=ID_n,inherit.aes = F,
               aes(x=102,y=(c(1:(nrow(ID_n)))-1)/(nrow(ID_n))*per_bott,
              label=names,vjust=0),size=2)+
  geom_text(data=ID_v,inherit.aes = F,
            aes(x=(c(1:nrow(ID_v))-0.5)/nrow(ID_v)*per_right+per_left,
                y=per_bott+2,
                label=names,vjust=0),size=2)+theme(legend.position = "none")
show(ppp)