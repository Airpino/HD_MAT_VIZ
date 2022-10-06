#BLOOD to tibble
library(HistDAWass)
library(tidyverse)
library(tibble)

# Trasforma una tibble in una math in base ad una variabile factor
tibble2MATH<-function(data, #a dataframe
                      factorv=1, # the factor variable used for the new row names
                      num_vars=c(1:3), # numerical variables to translate into distributions
                      qua=100 #a number of quantiles in wich divide the distribution
){
  
  np=c(0:qua)/qua
  
  data[[factorv]] <- factor(data[[factorv]])
  names_of_units <- as.character(levels(data[[factorv]]))
  names_of_vars <- colnames(data)[num_vars]
  
  data2<-data %>% select_at(c(factorv,num_vars))
  MyMath<-MatH(nrows=length(names_of_units),ncols=length(num_vars),
               rownames = names_of_units, varnames=names_of_vars)
  for (i in 1:length(names_of_units)){
    
    for (j in 1:length(names_of_vars)){
      #  browser()
      tmp<-data2 %>% filter(.[[1]]==names_of_units[i]) %>% select_at(j+1)
      
      tmpx<-unname(quantile(unlist(tmp),p=np))
      tmpd <-distributionH(x=tmpx,p=np)
      MyMath@M[i,j][[1]]<-tmpd
    }
  }
  return(MyMath)
}
# Trasforma una MATH in un oggetto tibble
MATH2tibble<-function(data){
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

# per una sola variabile da una tibble 
# si estrae un dataframe che descrive le density functions
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

# per un set di variabili da una tibble 
# si estrae un dataframe che descrive le density functions
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
    DF2<-one_density(A,v=j,poi=poi,pp=pp, n=n,NORM=NORM)
    DF2<-DF2 %>% mutate(VAR=rep(as.character(v[j]),nrow(DF2)))
    DF1<-rbind(DF1,DF2)
  }
  DF1$VAR=factor(DF1$VAR,levels=as.character(v))
  DF1$ID_VAR=factor(DF1$ID_VAR,levels=vnames)
  DF1$ID_name=factor(DF1$ID_name,levels=id_names)

  return(DF1)
}

# per una sola variabile da una tibble 
# si estrae un dataframe che descrive le frequenze 
# negli n bins associate ai centri
one_density_freq<-function(A,v=1,n=50,NORM=T, NX=T){
  ID_names<-unlist(A[,ncol(A)])
  A<-A %>% select(all_of(v))
  tmp<-unnest(A,cols = 1)
  ran<-range(tmp$x)
  vals <- seq(from=ran[1],to=ran[2],length.out=n+1)
  doms_c<-(vals[1:(length(vals)-1)]+vals[2:(length(vals))])*0.5
  nr<-nrow(A)
  my_x<-my_y<-numeric()
  my_id<-name_var<-name_ID<-character()
  for(i in 1:nr){
    tmp_d<-distributionH(x=A[,1][[i]]$x,
                         p=A[,1][[i]]$y)
    freqs<-diff(sapply(vals,function(i) compP(tmp_d,q=i)))
    my_x<-c(my_x,doms_c)
    my_y<-c(my_y,freqs)
    my_id<-c(my_id,rep(as.character(i),length(freqs)))
    name_ID<-c(name_ID,rep(ID_names[i],length(freqs)))
    name_var<-c(name_var,rep(names(A)[1],length(freqs)))
  }
  my_y_old<-my_y
  if(NORM) {my_y<-(my_y/max(my_y))}
  if(NX) NX<-((my_x-min(my_x))/diff(range(my_x)))
  df0<-data.frame(x=my_x,y=my_y,y_orig=my_y_old,
                  ID=factor(my_id,levels = as.character(c(1:nr))),
                  ID_VAR=name_var,ID_name=name_ID,
                  NX=NX)
  # df0<-data.frame(x=my_x,y=my_y,
  #                 ID=factor(my_id,levels = as.character(c(1:nr))),
  #                 ID_VAR=name_var)
  return(df0)
}

# per molte variabili da una tibble 
# si estrae un dataframe che descrive le frequenze 
# negli n bins associate ai centri
mul_freq<-function(A,v=c(1,2),poi=1000,pp=200, n=256,NORM=TRUE){
  v<-c(v,ncol(A))
  A<-A %>% select(all_of(v))
  NVARS<-length(v)-1
  vnames<-names(A)[1:(length(v)-1)]
  id_names<-unlist(A$ID_name)
  DF1<-one_density_freq(A,v=1, n=n,NORM=NORM)
  #browser()
  DF1<-DF1 %>% mutate(VAR=rep(as.character(v[1]),nrow(DF1)))
  
  for (j in 2:NVARS){
    DF2<-one_density_freq(A,v=j, n=n,NORM=NORM)
    DF2<-DF2 %>% mutate(VAR=rep(as.character(v[j]),nrow(DF2)))
    DF1<-rbind(DF1,DF2)
  }
  DF1$VAR=factor(DF1$VAR,levels=as.character(v))
  DF1$ID_VAR=factor(DF1$ID_VAR,levels=vnames)
  DF1$ID_name=factor(DF1$ID_name,levels=id_names)
  
  return(DF1)
}

# usa le densità. Prepara il dataframe per plottare mettendo 
# in ordine righe e colonne in base al clustering
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

# usa le frequenze. Prepara il dataframe per plottare mettendo 
# in ordine righe e colonne in base al clustering
prep_for_plot_freq<-function(B,v=c(1,2),n=64,NORM=TRUE,
                             order_of_rows=NA,order_of_cols=NA){
  df2<-mul_freq(B,v = v,n=n,NORM=NORM)
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

# crea il plot 
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


#extract kde from distributional data
Extr_dens<-function(distr, poi=200){
  vec<-HistDAWass:::COMP_MQ(distr,c(0:poi)/poi)
  mm<-density(vec)
}


