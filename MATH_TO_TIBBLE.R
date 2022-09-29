#BLOOD to tibble
MATH2tibble<-function(data){
  require(HistDAWass)
  require(tidyverse)
  require(tibble)
  #data<-BLOOD
  r<-nrow(data@M)
  c<-ncol(data@M)
  A<-tibble()
  tmpl<-list()
  for (j in 1:c){
    tmpl[[j]]<-list()
    for (i in 1:r){
      tmpl[[j]][[i]]<-tibble(x=data@M[i,j][[1]]@x,y=data@M[i,j][[1]]@p)
    }
    
  }
  
  A<-tibble(V1=tmpl[[1]])
  for (j in 2:c){
    tmp<-tibble(V2=tmpl[[j]]) 
    tmp<-dplyr::rename(tmp, !!paste0("V",j) := V2)
    
    #    paste0("V",j)=tmpl[[j]])
    A<-cbind(A,tmp)
  }
  
  names(A)<-get.MatH.varnames(data)
  A<-A %>% mutate(ID_name=get.MatH.rownames(data))
  
  return(A)
}


one_density<-function(A,v=1,poi=1000,pp=200, n=1024,NORM=TRUE,NX=TRUE){
  ID_names<-unlist(A[,ncol(A)])
  A<-A %>% select(all_of(v))
  tmp<-unnest(A,cols = 1)
  ran<-range(tmp$x)
  vals <- seq(from=ran[1],to=ran[2],length.out=poi)
  nr<-nrow(A)
  my_x<-my_y<-numeric()
  my_id<-var_name<-name_ID<-character()
  
  for(i in 1:nr){
    qua<-HistDAWass:::compQ_vect(distributionH(x=A[,1][[i]]$x,
                                               p=A[,1][[i]]$y),
                                 vp=c(0:pp)/pp)
    tmp1<-density(qua,
                  from=ran[1], to=ran[2],n = n)
    my_x<-c(my_x,tmp1$x)
    my_y<-c(my_y,tmp1$y)
    my_id<-c(my_id,rep(as.character(i),length(tmp1$x)))
    var_name<-c(var_name,rep(names(A)[1],length(tmp1$x)))
    name_ID<-c(name_ID,rep(ID_names[i],length(tmp1$x)))
    # df0<-data.frame(x=c(tmp1$x,vals),
    #                 y=c(tmp1$y,rep(NA,length(vals))))
    # 
  }
  
  if(NORM) my_y<-(my_y/max(my_y))#my_y<-sqrt(my_y/max(my_y))
  if(NX) NX<-((my_x-min(my_x))/diff(range(my_x)))
  df0<-data.frame(x=my_x,y=my_y,
                  ID=factor(my_id,levels = as.character(c(1:nr))),
                  ID_VAR=var_name,ID_name=name_ID,
                  NX=NX)
  return(df0)
}
one_density_freq<-function(A,v=1,poi=1000,pp=200, n=1024){
  
  A<-A %>% select(all_of(v))
  tmp<-unnest(A,cols = 1)
  ran<-range(tmp$x)
  vals <- seq(from=ran[1],to=ran[2],length.out=poi)
  nr<-nrow(A)
  my_x<-my_y<-numeric()
  my_id<-name_var<-name_ID<-character()
  for(i in 1:nr){
    qua<-HistDAWass:::compQ_vect(distributionH(x=A[,1][[i]]$x,
                                               p=A[,1][[i]]$y),
                                 vp=c(0:pp)/pp)
    tmp1<-density(qua,
                  from=ran[1], to=ran[2],n = n)
    # sapply(vals,function(i)compP(tmp,i)) for computing probs
    #browser()
    my_x<-c(my_x,tmp1$x)
    my_y<-c(my_y,tmp1$y)
    my_id<-c(my_id,rep(as.character(i),length(tmp1$x)))
    browser()
    name_var<-c(name_var,rep(names(A)[1],length(tmp1$x)))
    # df0<-data.frame(x=c(tmp1$x,vals),
    #                 y=c(tmp1$y,rep(NA,length(vals))))
    # 
  }
  #browser()
  if(NORM) my_y<-sqrt(my_y/max(my_y))
  df0<-data.frame(x=my_x,y=my_y,
                  ID=factor(my_id,levels = as.character(c(1:nr))),
                  ID_VAR=name_var)
  return(df0)
}
mul_density<-function(A,v=c(1,2),poi=1000,pp=200, n=256,NORM=TRUE){
  v<-c(v,ncol(A))
  A<-A %>% select(all_of(v))
  NVARS<-length(v)-1
  vnames<-names(A)[1:(length(v)-1)]
  id_names<-unlist(A$ID_name)
  DF1<-one_density(A,v=1,poi=poi,pp=pp, n=n,NORM=NORM)
  #browser()
  DF1<-DF1 %>% mutate(VAR=rep(as.character(v[1]),nrow(DF1)))
  
  for (j in 2:NVARS){
    DF2<-one_density(A,v=j,poi=poi,pp=pp, n=n)
    DF2<-DF2 %>% mutate(VAR=rep(as.character(v[j]),nrow(DF2)))
    DF1<-rbind(DF1,DF2)
  }
  DF1$VAR=factor(DF1$VAR,levels=as.character(v))
  DF1$ID_VAR=factor(DF1$ID_VAR,levels=vnames)
  DF1$ID_name=factor(DF1$ID_name,levels=id_names)

  return(DF1)
}

plots_strips_distr_2<-function(B,v=c(1,2),n=64,poi=500,pp=200,NORM=TRUE,
                               order_of_rows=NA,order_of_cols=NA){
  
  require(RColorBrewer)
  
  col_strip <- brewer.pal(9, "Reds")
  df2<-prep_for_plot(B,v,n,poi,pp,NORM,
                     order_of_rows,order_of_cols)
  
  
  #browser()
  
  
    plots<-ggplot(df2 ,
                       aes(x = x, y = 1, fill = y))+
      geom_tile(show.legend = T)+#show.legend = F)+
      scale_fill_gradientn(colors = (col_strip))+
      facet_grid(rows = vars(ID_name),
                             cols =vars(ID_VAR),
                  scales = "free_x"
      )+
      theme_few() +
      theme(plot.margin=margin(r=2.5,unit = "cm"),
        panel.spacing.y=unit(0, "lines"),
        #strip.text =element_text(margin = margin(r=20)),
        strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
        strip.text.x = element_text(size = 7, angle = 0),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        legend.position='none')
  
  return(plots)
}

prep_for_plot<-function(B,v=c(1,2),n=64,poi=500,pp=200,NORM=TRUE,
                        order_of_rows=NA,order_of_cols=NA){
  df2<-mul_density(B,v = v,n=n,poi=poi,pp = pp,NORM=NORM)
  if(!is.na(order_of_rows[1])){
    new_levels<-levels(df2$ID_name)[order_of_rows]
    df2$ID_name<-factor(df2$ID_name,levels=new_levels)
  }
  if(!is.na(order_of_cols[1])){
    new_levels<-levels(df2$ID_VAR)[order_of_cols]
    df2$ID_VAR<-factor(df2$ID_VAR,levels=new_levels)
  }
  return(df2)
}


plots_strips_distr_3<-function(B,v=c(1,2),n=64,poi=500,pp=200,NORM=TRUE,
                               order_of_rows=NA,order_of_cols=NA){
  
  require(RColorBrewer)
  
  col_strip <- brewer.pal(9, "Reds")
  df2<-prep_for_plot(B,v,n,poi,pp,NORM,
                     order_of_rows,order_of_cols)

  
  #browser()
  scale_fill_gradient(low = "blue", high = "red", na.value = NA)
  plots<-ggplot(df2 ,
                aes(x = x, y = 1, fill = NX))+
    geom_raster(show.legend = T,aes(alpha=y))+#show.legend = F)+
    scale_fill_gradient(low="blue",high="red")+
    facet_grid(rows = vars(ID_name),
               cols =vars(ID_VAR),
               scales = "free_x"
    )+
    theme_few() +
    theme(plot.margin=margin(r=2.5,unit = "cm"),
          panel.spacing.y=unit(0, "lines"),
          #strip.text =element_text(margin = margin(r=20)),
          strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
          strip.text.x = element_text(size = 7, angle = 0),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 7, angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_blank(),
          legend.position='none')
  
  return(plots)
}



library(ggthemes)

runfollow=F
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
  df2<-one_density(B,v = 8,n=64)
  df2<-mul_density(B,v = c(5:8),n=64)
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
