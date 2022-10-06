#BLOOD to tibble
library(HistDAWass)
library(tidyverse)
library(tibble)
library(ggthemes)
source("functions.R")
runfollow=F #prova gli Horizon plots
if (runfollow){
  library(ggHoriPlot)
  
  df0<-one_density(B,v = 8,n=256)
  cutpoints <- df0  %>% 
    mutate(
      outlier = between(
        y, 
        quantile(y, 0.25, na.rm=T)-
          1.5*IQR(y, na.rm=T),
        quantile(y, 0.75, na.rm=T)+
          1.5*IQR(y, na.rm=T))) %>% 
    filter(outlier)
  
  ori <- sum(range(cutpoints$y))/2
  sca <- seq(range(cutpoints$y)[1], 
             range(cutpoints$y)[2], 
             length.out = 7)[-4]
  show(df0 %>% ggplot() +
         geom_horizon(aes(x, 
                          y,
                          fill = ..Cutpoints..), 
                      origin = ori, horizonscale = sca,show.legend = F) +
         scale_fill_hcl(palette = 'Reds 3', reverse = T)+
         facet_grid(ID~.)+
         theme_few() +
         theme(
           panel.spacing.y=unit(0, "lines"),
           strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
           axis.text.y = element_blank(),
           axis.title.y = element_blank(),
           axis.ticks.y = element_blank(),
           panel.border = element_blank(),
           panel.background = element_rect(fill = "lightblue")#,
           #                                           colour = "lightblue",
           #                                           size = 0.5, linetype = "solid")
         )
  )
}

# una sola variabile
# qui, invece calcolo i rettangoli da colorare (funziona con le densità)
create_new_df_to_plot<-function(df2,col=1){
  # $ x      : num  -0.373 2.663 5.699 8.736 11.772 ...
  # $ y      : num  0.00727 0.0173 0.03051 0.04325 0.05869 ...
  # $ ID     : Factor w/ 78 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
  # $ ID_VAR : Factor w/ 4 levels "Ozone.Conc.ppb",..: 1 1 1 1 1 1 1 1 1 1 ...
  # $ ID_name: Factor w/ 78 levels "I1","I2","I3",..: 1 1 1 1 1 1 1 1 1 1 ...
  # $ NX     : num  0 0.0204 0.0408 0.0612 0.0816 ...
  # $ VAR    : Factor w/ 5 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
  df<-df2 %>% filter(ID_VAR==levels(df2$ID_VAR)[col])
  # quanti valori per ogni osservazione
  N<-nrow(df %>% filter(ID_name==levels(df2$ID_name)[1]))
  range_of_var<-range(df$x)
  x<-seq(from=range_of_var[1],to=range_of_var[2],length.out=N+1)
  x_ini<-rep(x[1:N],length(levels(df$ID_name)))
  x_fin<-rep(x[2:(N+1)],length(levels(df$ID_name)))
  xx<-seq(0,1,length.out=N)
  myc<-rep(scales::seq_gradient_pal(low="red",high="blue")(xx),length(levels(df$ID_name)))
  y_ini<-y_fin<-alp<-numeric()
  ID_name<-character()
  c<-0
  for(i in levels(df2$ID_name)){
    c<-c+1
    y_ini<-c(y_ini,rep(c-1,N))
    y_fin<-c(y_fin,rep(c,N))
    alp<-c(alp,round(unlist(df %>% filter(ID_name==i)%>% select(y)),3))
    ID_name<-c(ID_name,rep(i,N))
  }
  #ID_var<-
  #browser()
  df<-data.frame(x_ini=x_ini,x_fin=x_fin,
                 y_ini=y_ini,y_fin=y_fin,
                 cols=myc,alp=alp,
                 ID_name=factor(ID_name,levels=levels(df2$ID_name)),
                 VAR_name=rep(levels(df2$ID_VAR)[col],length(x_ini)))
  df<-df %>% mutate(alp=if_else(alp==0,NA_real_,alp),cols=if_else(alp==0,"white",cols))
  return(df)
}

# tutte le variabili
# qui, invece calcolo i rettangoli da colorare (funziona con le densità)
create_new_df_to_plot_M<-function(df2){
  
  DDF<-create_new_df_to_plot(df2,col=1)
  c=1
  for(j in levels(df2$ID_VAR)[2:length(levels(df2$ID_VAR))]){
    c=c+1
    DDFtmp<-create_new_df_to_plot(df2,col=c)
    DDF<-rbind(DDF,DDFtmp)
  }
  DDF$VAR_name<-factor(DDF$VAR_name,levels=levels(df2$ID_VAR))
  return(DDF)
}

# tutte le variabili
# qui, invece calcolo i rettangoli da colorare (funziona con le frequenze)
create_new_df_to_plot_M_freq<-function(df2){
  
  DDF<-create_new_df_to_plot(df2,col=1)
  c=1
  for(j in levels(df2$ID_VAR)[2:length(levels(df2$ID_VAR))]){
    c=c+1
    DDFtmp<-create_new_df_to_plot(df2,col=c)
    DDF<-rbind(DDF,DDFtmp)
  }
  DDF$VAR_name<-factor(DDF$VAR_name,levels=levels(df2$ID_VAR))
  return(DDF)
}
plots_strips_distr<-function(B,v=c(1,2),n=64,poi=500,pp=200,NORM=TRUE){
  
  plots<-list()
  require(RColorBrewer)
  
  col_strip <- brewer.pal(9, "Reds")
  df2<-mul_density(B,v = v,n=n,poi=poi,pp = pp,NORM=NORM)
  
  for(j in 1:length(v)){
    
    plots[[j]]<-ggplot(df2 %>% filter(VAR==as.character(v[j])),
                       aes(x = x, y = 1, fill = y))+
      geom_tile(show.legend = F)+#show.legend = F)+
      scale_fill_gradientn(colors = (col_strip))+
      facet_grid(rows = vars(ID),
                 #            cols =vars(VAR),
                 #            scales = "free_y"
      )+
      theme_few() +
      theme(
        panel.spacing.y=unit(0, "lines"),
        strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank())
  }
  return(plots)
}




####  APPLY -----
RUNME<-F
if (RUNME){
  library(RColorBrewer)
  theme_strip <- theme_minimal()+
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_blank(),
          axis.text.x = element_text(vjust = 3),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 14, face = "bold"),
          panel.spacing.y=unit(0, "lines"),
          strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_blank()
    )
  
  
  col_strip <- brewer.pal(9, "Reds")
  
  brewer.pal.info
  #  df1 <- df0%>%filter(ID==1) 
  df2<-one_density(B,v = 1,n=64)
  df2<-mul_density(B,v = c(1:4),n=64)
  ggplot(df2,
         aes(x = x, y = 1, fill = y))+
    geom_tile()+#show.legend = F)+
    scale_fill_gradientn(colors = (col_strip))+
    facet_grid(rows = vars(ID)#,
               #            cols =vars(VAR),
               #            scales = "free_y"
    )+
    theme_few() +
    theme(
      panel.spacing.y=unit(0, "lines"),
      strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank())
  #theme_strip
}

# library(easyGgplot2)
 # B<-MATH2tibble(China_Seas)
 # plo<-plots_strips_distr(B,v=c(5:8),poi = 50,pp = 50,n = 32)
 # show(easyGgplot2::ggplot2.multiplot(plo[[1]]+theme(strip.text.y = element_blank()) ,
 #                                     plo[[2]]+theme(strip.text.y = element_blank()),
 #                                     plo[[3]]+theme(strip.text.y = element_blank()),
 #                                     plo[[4]],cols=4))
# library(patchwork)
