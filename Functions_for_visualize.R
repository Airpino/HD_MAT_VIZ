library(HistDAWass)
library(tidyverse)
library(ggtext)
library(ggfittext)
library(plotly)
source("functions.R")
## this function takes a MATH and discretizes the distributions
## by cutting the domain of each variable (col) into n equi-width classes
## Can be improved for avoiding outliers!(TO BE DONE!)
discretize_data<-function(MAT,variable=1,n=10,absolute=F){
  TIB<-MATH2tibble(MAT)
  TIBN<-TIB
  for(i in 1:ncol(TIB)){
    if(is(TIB[[i]])[[1]]=="list") class(TIB[[i]])="continuous distr"
  }
  
  if(absolute){
    min_ab<- Inf
    max_ab<- -Inf
    for(i in 1:ncol(TIB)){
      if(is(TIB[[i]])[[1]]=="continuous distr"){
        t_min_dom<-min(sapply(TIB[[i]], function(x) min(x[,1])))
        t_max_dom<-max(sapply(TIB[[i]], function(x) max(x[,1])))
        if(t_min_dom<min_ab) min_ab<-t_min_dom
        if(t_max_dom>max_ab) max_ab<-t_max_dom
      }
    }
    
    for(i in 1:ncol(TIB)){
      if(is(TIB[[i]])[[1]]=="continuous distr"){
        min_dom<-min_ab
        max_dom<-max_ab
        limits <- seq(min_dom,max_dom,length.ou=n+1)
        for(r in 1:nrow(TIB)){
          tmp<-MAT@M[r,i][[1]]
          cdf<-sapply(limits,function(x) compP(tmp,x))
          TIBN[[i]][[r]]=data.frame(xmin=limits[1:n],
                                    xmax=limits[2:(n+1)],
                                    cl_lab=paste0("[",round(limits[1:n],3),"-",round(limits[2:(n+1)],3),"]"),
                                    cdf_min=cdf[1:n],
                                    cdf_max=cdf[2:(n+1)],
                                    cdf=cdf[2:(n+1)],
                                    x_cat=c(1:n),
                                    freq=cdf[2:(n+1)]-cdf[1:n])%>% filter(freq>0) %>% 
            select(x_cat,freq,cdf,everything())
        }
        
      }
    }
  }else{
    for(i in 1:ncol(TIB)){
      
      if(is(TIB[[i]])[[1]]=="continuous distr"){
        min_dom<-min(sapply(TIB[[i]], function(x) min(x[,1])))
        max_dom<-max(sapply(TIB[[i]], function(x) max(x[,1])))
        limits <- seq(min_dom,max_dom,length.ou=n+1)
        for(r in 1:nrow(TIB)){
          tmp<-MAT@M[r,i][[1]]
          cdf<-sapply(limits,function(x) compP(tmp,x))
          TIBN[[i]][[r]]=data.frame(xmin=limits[1:n],
                                    xmax=limits[2:(n+1)],
                                    cl_lab=paste0("[",round(limits[1:n],3),"-",round(limits[2:(n+1)],3),"]"),
                                    cdf_min=cdf[1:n],
                                    cdf_max=cdf[2:(n+1)],
                                    cdf=cdf[2:(n+1)],
                                    x_cat=c(1:n),
                                    freq=cdf[2:(n+1)]-cdf[1:n])%>% filter(freq>0) %>% 
            select(x_cat,freq,cdf,everything())
        }
        
      }
    }
  }
  return(list(TIBN=TIBN,oriTIB=TIB)) 
}


discretize_normalized_data<-function(MAT,n=10){
  #we cut -3 3
  STD_H<-sapply(c(1:get.MatH.ncols(MAT)), 
                function(x) sqrt(HistDAWass::WH.var.covar(MAT[,x])))
  # sapply(c(1:get.MatH.ncols(MAT)), function(x) WH.vec.mean(data[,x])@m)
  MM_H<-sapply(c(1:get.MatH.ncols(MAT)), function(x) mean(sapply(c(1:get.MatH.nrows(MAT)),function(y)MAT[y,x]@M[[1]]@m)))
  
  for (i in 1:get.MatH.nrows(MAT)){
    for(j in 1:get.MatH.ncols(MAT)){
      MAT@M[i,j][[1]]@x<-(MAT@M[i,j][[1]]@x-MM_H[j])/STD_H[j]
      MAT@M[i,j][[1]]@m<-(MAT@M[i,j][[1]]@m-MM_H[j])/STD_H[j]
      MAT@M[i,j][[1]]@s<-MAT@M[i,j][[1]]@s/STD_H[j]
      
    }
  }
  
  TIB<-MATH2tibble(MAT)
  TIBN<-TIB
  for(i in 1:ncol(TIB)){
    if(is(TIB[[i]])[[1]]=="list") class(TIB[[i]])="continuous distr"
  }
  
  #if(absolute){
  min_ab<- -Inf
  max_ab<- Inf
  limits <- c(-Inf,seq(-3,3,length.out=n-1),Inf)
  for(i in 1:ncol(TIB)){
    if(is(TIB[[i]])[[1]]=="continuous distr"){
      min_dom<-min_ab
      max_dom<-max_ab
      
      for(r in 1:nrow(TIB)){
        tmp<-MAT@M[r,i][[1]]
        cdf<-sapply(limits,function(x) compP(tmp,x))
        TIBN[[i]][[r]]=data.frame(xmin=limits[1:n],
                                  xmax=limits[2:(n+1)],
                                  cl_lab=paste0("[",round(limits[1:n],3)," ; ",round(limits[2:(n+1)],3),"]"),
                                  cdf_min=cdf[1:n],
                                  cdf_max=cdf[2:(n+1)],
                                  cdf=cdf[2:(n+1)],
                                  x_cat=c(1:n),
                                  freq=cdf[2:(n+1)]-cdf[1:n])%>% filter(freq>0) %>% 
          select(x_cat,freq,cdf,everything())
      }
      
    }
  }
  #}else{
  # for(i in 1:ncol(TIB)){
  #   
  #   if(is(TIB[[i]])[[1]]=="continuous distr"){
  #     min_dom<-min(sapply(TIB[[i]], function(x) min(x[,1])))
  #     max_dom<-max(sapply(TIB[[i]], function(x) max(x[,1])))
  #     limits <- seq(min_dom,max_dom,length.ou=n+1)
  #     for(r in 1:nrow(TIB)){
  #       tmp<-MAT@M[r,i][[1]]
  #       cdf<-sapply(limits,function(x) compP(tmp,x))
  #       TIBN[[i]][[r]]=data.frame(xmin=limits[1:n],
  #                                 xmax=limits[2:(n+1)],
  #                                 cl_lab=paste0("[",round(limits[1:n],3),"-",round(limits[2:(n+1)],3),"]"),
  #                                 cdf_min=cdf[1:n],
  #                                 cdf_max=cdf[2:(n+1)],
  #                                 cdf=cdf[2:(n+1)],
  #                                 x_cat=c(1:n),
  #                                 freq=cdf[2:(n+1)]-cdf[1:n])%>% filter(freq>0) %>% 
  #         select(x_cat,freq,cdf,everything())
  #     }
  #     
  #   }
  # }
  #}
  return(list(TIBN=TIBN,oriTIB=TIB,NEWMAT=MAT)) 
}




## Squish data
squish_data<-function(Tib,labels=c(1:nrow(Tib))) {
  Tib<-as_tibble(Tib)
  n<-nrow(Tib)
  col<-ncol(Tib)
  colna<-colnames(Tib)
  IDr<-IDc<-valueD<-cdf<-freq<-numeric()
  labID<-labVar<-character()
  for(i in 1:n){
    for(j in 1:col){
      tmpx<-Tib[[j]][[i]][[1]]
      tmpcdf<-Tib[[j]][[i]][["cdf"]]
      tmpf<-Tib[[j]][[i]][["freq"]]
      tmpl<-length(tmpx)
      IDr<-c(IDr,rep(i,tmpl))
      IDc<-c(IDc,rep(j,tmpl))
      labID<-c(labID,rep(as.character(labels[i]),tmpl))
      labVar<-c(labVar,rep(colna[j],tmpl))
      valueD<-c(valueD,tmpx)
      freq<-c(freq,tmpf)
      cdf<-c(cdf,tmpcdf)
    }
  }
  DF<-data.frame(IDr=IDr,IDc=IDc,
                 labID=labID,labVar=labVar,
                 dom=valueD, freq=freq,cdf=cdf)
  return(DF)  
} 


## Create the new green eyes plots

GEI_plot<-function(Tib,selected=1,
                   labels=c(1:nrow(Tib)),TITLE=T,notick=F,
                   skewness_plo=FALSE,
                   alpha=1,
                   bg="transparent",BW=F,
                   polar=T,legend=F,pupil=0.2, pupil.color="black"){
  DF<-squish_data(Tib[selected,],labels = labels[selected])
  if(skewness_plo){
    bstat<-DF %>% 
      group_by(labVar) %>% 
      summarize(
        m=sum(dom*freq),
        s=sqrt(sum(dom^2*freq)-m^2),
        sk=(sum(((dom-m)/s)^3*freq))
      ) %>% ungroup() %>% mutate(sk=sign(sk)*abs(sk)^(1/3), coord=(sk+3)/6) 
  }
  DF$labVar<-factor(DF$labVar,levels=colnames(Tib))
  #stacked perc bars
  # cc<-c( "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B","#D9EF8B", 
  #                 "#A6D96A" ,"#66BD63", "#1A9850", "#006837")
  if (BW){
    
    cols <- c("0" = "#fefefe", "1" = "#f1f1f1", "2" = "#e1e1e1", "3" = "#d9d9d9",
              "4" = "#dddddd", "5" = "#cccccc", "6" = "#a1a1a1", "7" = "#818181",
              "8" = "#666666", "9" = "#414141", "10" = "#000000")  
  }else{
    cols <- c("0" = "#A50026", "1" = "#A50026", "2" = "#D73027", "3" = "#F46D43",
              "4" = "#FDAE61", "5" = "#FEE08B", "6" = "#D9EF8B", "7" = "#66BD63",
              "8" = "#1A9850", "9" = "#006837", "10" = "#003300")
  }
  p<-ggplot(DF,aes(x=labVar,y=freq))
  if (legend){
    p<-p+geom_bar(stat="identity",
                  aes(fill=fct_rev(as.factor(dom))),
                  width=1,alpha=alpha)
  }else{  
    p<-p+geom_bar(stat="identity",
                  aes(fill=fct_rev(as.factor(dom))),
                  width=1,
                  show.legend = F,alpha=alpha)
  }
  
  p<-p+ylim(c(-pupil,1))+
    geom_rect(aes(xmin=0.5,ymin = -pupil,xmax=max(IDc)+.5,ymax=0),fill=pupil.color)+
    scale_fill_manual(
      values = cols)+
    #    scale_fill_brewer(palette = "Spectral",direction=1)+
    #    scale_fill_brewer(palette = "RdYlGn",direction=1)+
    geom_hline(yintercept = 0.5,
               linetype="dashed",
               #color="white",
               #size=1,
               alpha=0.6)
  
  if(skewness_plo){
    p<-p+geom_point(inherit.aes = F,data=bstat,
                    aes(x=labVar,y=coord),
                    color="black",fill="pink",shape=21,
                    alpha=0.6,show.legend = F)
  }
  
  if(polar) p<-p+coord_polar()+theme_minimal()
  p<-p+theme_minimal()
  if(TITLE) p<-p+
    labs(#title = labels[selected]#,
      #subtitle = "subtitle",
      caption = labels[selected]
    )+
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5))
  if(notick){
    p<-p+
      theme(axis.text.x=element_blank())
  }
  p<-p+
    theme(axis.line=element_blank(),
          #        axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.border=element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(color="transparent",fill = "transparent"),
          panel.background = element_rect(color="transparent",fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = bg, color = NA),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0))
  return(p) 
}


GEI_plot2<-function(Tib,selected=1,
                    labels=c(1:nrow(Tib)),TITLE=T,notick=F,
                    skewness_plo=FALSE,
                    alpha=1,
                    bg="transparent",BW=F,
                    polar=T,legend=F,pupil=0.2, pupil.color="black",levels_of_colors=10,
                    palette=0){
  DF<-squish_data(Tib[selected,],labels = labels[selected])
  if(skewness_plo){
    bstat<-DF %>% 
      group_by(labVar) %>% 
      summarize(
        m=sum(dom*freq),
        s=sqrt(sum(dom^2*freq)-m^2),
        sk=(sum(((dom-m)/s)^3*freq))
      ) %>% ungroup() %>% mutate(sk=sign(sk)*abs(sk)^(1/3), coord=(sk+3)/6) 
  }
  DF$labVar<-factor(DF$labVar,levels=colnames(Tib))
  #stacked perc bars
  # cc<-c( "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B","#D9EF8B", 
  #                 "#A6D96A" ,"#66BD63", "#1A9850", "#006837")
  if (BW){
    
    cols <- c("0" = "#fefefe", "1" = "#f1f1f1", "2" = "#e1e1e1", "3" = "#d9d9d9",
              "4" = "#dddddd", "5" = "#cccccc", "6" = "#a1a1a1", "7" = "#818181",
              "8" = "#666666", "9" = "#414141", "10" = "#000000")  
  }else{
    cc<-as.vector(paletteer_d(`"MetBrewer::Paquin"`,n = levels_of_colors,type="continuous"))
    if (palette==1) cc<-as.vector(paletteer_c("grDevices::RdYlGn", levels_of_colors))
    if (palette==2) cc<-as.vector(paletteer_d("MetBrewer::Benedictus",n = levels_of_colors,type="continuous"))
    if (palette==3) cc<-as.vector(paletteer_d("MetBrewer::Hiroshige",n = levels_of_colors,type="continuous"))
    if (palette==4) cc<-as.vector(paletteer_d("MetBrewer::Hiroshige",n = levels_of_colors,type="continuous",direction=-1))
    if (palette==5) cc<-as.vector(paletteer_d("MetBrewer::Tiepolo",n = levels_of_colors,type="continuous"))
    if (palette==6) cc<-as.vector(paletteer_d("PNWColors::Bay",n = levels_of_colors,direction=-1,type="continuous"))
    if (palette==7) cc<-as.vector(paletteer_c("scico::vik",n = levels_of_colors))
    
    
    #as.vector(paletteer_c("grDevices::RdYlGn", levels_of_colors))
    #as.vector(paletteer_d(`"feathers::rose_crowned_fruit_dove"`,n = levels_of_colors,type="continuous"))
    as.vector(paletteer_d(`"MetBrewer::Paquin"`,n = levels_of_colors,type="continuous"))
    names(cc)<-c(1:levels_of_colors)
    #to try "grDevices::Fall"  grDevices::Tropic (rev)  RdYlGn
    
    #DF<-DF %>% left_join(cc,by=c("dom"="val"))
    
    cols <- c("0" = "#A50026", "1" = "#A50026", "2" = "#D73027", "3" = "#F46D43",
              "4" = "#FDAE61", "5" = "#FEE08B", "6" = "#D9EF8B", "7" = "#66BD63",
              "8" = "#1A9850", "9" = "#006837", "10" = "#003300")
  }
  p<-ggplot(DF,aes(x=labVar,y=freq))
  if (legend){
    p<-p+geom_bar(stat="identity",
                  aes(fill=fct_rev(as.factor(dom))),
                  width=1,alpha=alpha)
  }else{  
    p<-p+geom_bar(stat="identity",
                  aes(fill=fct_rev(as.factor(dom))),
                  width=1,
                  show.legend = F,alpha=alpha)
  }
  
  p<-p+ylim(c(-pupil,1))+
    geom_rect(aes(xmin=0.5,ymin = -pupil,xmax=max(IDc)+.5,ymax=0),fill=pupil.color)+
    scale_fill_manual(
      values = cc)+
    #    scale_fill_brewer(palette = "Spectral",direction=1)+
    #    scale_fill_brewer(palette = "RdYlGn",direction=1)+
    geom_hline(yintercept = 0.5,
               linetype="dashed",
               #color="white",
               #size=1,
               alpha=0.6)
  
  if(skewness_plo){
    p<-p+geom_point(inherit.aes = F,data=bstat,
                    aes(x=labVar,y=coord),
                    color="black",fill="black",shape=21,
                    alpha=0.6,show.legend = F)
  }
  
  if(polar) p<-p+coord_polar()+theme_minimal()
  p<-p+theme_minimal()
  if(TITLE) p<-p+
    labs(#title = labels[selected]#,
      #subtitle = "subtitle",
      caption = labels[selected]
    )+
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5))
  if(notick){
    p<-p+
      theme(axis.text.x=element_blank())
  }
  p<-p+
    theme(axis.line=element_blank(),
          #        axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.border=element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(color="transparent",fill = "transparent"),
          panel.background = element_rect(color="transparent",fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = bg, color = NA),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0))
  return(p) 
}

Signed_MAD<-function(x,p,type="cont"){
  if(type=="cont")  {
    tmp<-distributionH(x=x,p=p)
    med<-compQ(tmp,0.5)
    x1<-c(x[x<med],med)
    p1<-c(p[which(x<med)],0.5)*2
    tmp1<-distributionH(x=x1,p=p1)
    x2<-c(med, x[x>med])
    p2<-(c(0.5,p[which(x>med)])-0.5)*2
    tmp2<-distributionH(x=x2,p=p2)
    left_MAD<-med-compQ(tmp,0.25)#compQ(tmp1,0.5)
    right_MAD<-compQ(tmp,0.75)-med#compQ(tmp2,0.5)-med
    
    M<-HistDAWass::meanH(tmp)
    pval<-compP(tmp,M)
    x1_1<-c(x[x<M],M)
    p1_1<-c(p[which(x<M)],pval)/pval
    tmp1_1<-distributionH(x=x1_1,p=p1_1)
    x2_1<-c(M, x[x>M])
    p2_1<-(c(pval,p[which(x>M)])-pval)/(1-pval)
    tmp2_1<-distributionH(x=x2_1,p=p2_1)
    L_B<- pval*(HistDAWass::meanH(tmp1_1)-M)^2
    R_B<- (1-pval)*(HistDAWass::meanH(tmp2_1)-M)^2
    L_W<-pval*HistDAWass::stdH(tmp1_1)^2
    R_W<-(1-pval)*HistDAWass::stdH(tmp2_1)^2
  }
  if(type=="discr")  {
    #to be done!!
  }
  
  
  return(c(L_M=left_MAD,R_M=right_MAD, L_V=L_B+L_W,R_V=R_B+R_W))
}

Signed_VARS<-function(distr,med,M){
  # x1<-c(distr@x[distr@x<med],med)
  # p1<-c(distr@p[which(distr@x<med)],0.5)*2
  # tmp1<-distributionH(x=x1,p=p1)
  # x2<-c(med, distr@x[distr@x>med])
  # p2<-(c(0.5,distr@p[which(distr@x>med)])-0.5)*2
  # tmp2<-distributionH(x=x2,p=p2)
  left_MAD<-med-compQ(distr,0.25)#compQ(tmp1,0.5)
  right_MAD<-compQ(distr,0.75)-med#compQ(tmp2,0.5)-med
  
  pval<-compP(distr,M)
  x1_1<-c(distr@x[distr@x<M],M)
  p1_1<-c(distr@p[which(distr@x<M)],pval)/pval
  tmp1_1<-distributionH(x=x1_1,p=p1_1)
  x2_1<-c(M, distr@x[distr@x>M])
  p2_1<-(c(pval,distr@p[which(distr@x>M)])-pval)/(1-pval)
  tmp2_1<-tmp1_1
  tmp2_1@x<-x2_1
  tmp2_1@p<-p2_1
  #  distributionH(x=x2_1,p=p2_1)
  L_B<- pval*(HistDAWass::meanH(tmp1_1)-M)^2
  R_B<- (1-pval)*(HistDAWass::meanH(tmp2_1)-M)^2
  L_W<-pval*HistDAWass::stdH(tmp1_1)^2
  R_W<-(1-pval)*HistDAWass::stdH(tmp2_1)^2
  
  return(c(L_M=left_MAD,R_M=right_MAD, L_V=L_B+L_W,R_V=R_B+R_W))
}



#compute some informations about left and right variability wrt median and mean
compute_LR_VARS_cont<-function(MAT){#Hmat continuous 
  n<-get.MatH.nrows(MAT)
  p<-get.MatH.ncols(MAT)
  meas<-c("Mean","Median","SD","VAR", "SKEW", "L_MAD","R_MAD","L_VAR","R_VAR")
  res<-array(0,dim = c(n,p,length(meas)))
  for(i in 1:n){
    for(j in 1:p){
      tmp<-MAT@M[i,j][[1]]
      res[i,j,1]<-meanH(tmp)
      res[i,j,2]<-compQ(tmp,0.5)
      res[i,j,3]<-stdH(tmp)
      res[i,j,4]<-(res[i,j,3])^2
      res[i,j,5]<-skewH(tmp)
      tmp2<-Signed_VARS(tmp,res[i,j,2],res[i,j,1])
      res[i,j,6]<-tmp2[1]
      res[i,j,7]<-tmp2[2]
      res[i,j,8]<-tmp2[3]
      res[i,j,9]<-tmp2[4]
      
    }
  }
  
  dimnames(res) <- list(get.MatH.rownames(MAT), 
                        get.MatH.varnames(MAT),
                        meas)
  
  return(res)
}

draw_VARS<-function(res2,SDS=sapply(c(1:get.MatH.ncols(data)),function(x) sqrt(WH.var.covar((data[,x]))))){
  #meas<-c("Mean","Median","SD","VAR", "SKEW", "L_MAD","R_MAD","L_VAR","R_VAR")
  p<-dim(res2)[2]
  MAXvars<-sapply(as_tibble(res2[,,"VAR",drop=T]),max)
  MAXsd<-sapply(as_tibble(res2[,,"SD",drop=T]),max)
  Vars<-res2[,,"VAR",drop=T]
  Sds_s<-res2[,,"SD",drop=T]
  Perc_L_Vars<-res2[,,"L_VAR",drop=T]/Vars
  Perc_R_Vars<-res2[,,"R_VAR",drop=T]/Vars
  V1<-as.vector(Vars)
  L1<-as.vector(Perc_L_Vars)
  R1<-as.vector(Perc_R_Vars)
  SDS_s<-as.vector(Sds_s)
  IDr<-rep(c(1:(dim(res2)[1])),dim(res2)[2])
  Names_r<-rep(dimnames(res2)[[1]],dim(res2)[2])
  IDc<-rep(c(1:(dim(res2)[2])),each=dim(res2)[1])
  Names_c<-rep(dimnames(res2)[[2]],each=dim(res2)[1])
  MAX_V<-rep(MAXvars,each=dim(res2)[1])
  # browser()
  MAX_S<-rep(MAXsd,each=dim(res2)[1])
  DF<-data.frame(IDr,Names_r,IDc,Names_c,VAR=V1,SD=SDS_s,MAX_V,L1,R1,MAX_S)
  
  DF<-DF %>% mutate(xmin=IDc,xmax=IDc,
                    y1min=0,y1max=SD/MAX_S*0.25,
                    y2min=y1max, y2max=SD/MAX_S*0.25,PROP=SD/MAX_S*0.5)
  
  # pl<-ggplot(DF, aes(x=IDc))+
  #   geom_rect(aes(xmin=xmin,ymin = y1min,xmax=xmax,ymax=y1max),fill="white")+
  #   geom_rect(aes(xmin=xmin,ymin = y2min,xmax=xmax,ymax=y2max),fill="black")
  # 
  # ggplot(DF %>% filter(IDr==3), aes(x=IDc))+
  #   geom_segment(aes(x=(xmin+xmax)*0.5,xend=(xmin+xmax)*0.5,y = y1min+(1-PROP)*0.5,yend=y1max+(1-PROP)*0.5),
  #             color="grey80",size=3)+
  #   geom_segment(aes(x=(xmin+xmax)*0.5,xend=(xmin+xmax)*0.5,y = y2min+(1-PROP)*0.5,yend=y2max+(1-PROP)*0.5),
  #             color="black",size=3)+
  #   geom_segment(aes(x=0.5,xend=p+0.5,y=0,yend=0))+
  #   ylim(c(0,1))+theme_void()+
  #   coord_polar()
  # 
  return(DF)
}

GEI_plot3<-function(Tib,selected=1,
                    labels=c(1:nrow(Tib)),
                    TITLE=T,
                    notick=F,
                    skewness_plo=FALSE,
                    alpha=1,
                    bg="transparent",
                    BW=F,
                    polar=T,
                    legend=F,
                    pupil=0.2, pupil.color="black",
                    levels_of_colors=10,
                    palette=1, 
                    var_bar=T,
                    var_DF=NULL){
  DF<-squish_data(Tib[selected,],labels = labels[selected])
  if(skewness_plo){
    bstat<-DF %>% 
      group_by(labVar) %>% 
      summarize(
        m=sum(dom*freq),
        s=sqrt(sum(dom^2*freq)-m^2),
        sk=(sum(((dom-m)/s)^3*freq))
      ) %>% ungroup() %>% mutate(sk=sign(sk)*abs(sk)^(1/3), coord=(sk+3)/6) 
  }
  DF$labVar<-factor(DF$labVar,levels=colnames(Tib))
  #stacked perc bars
  # cc<-c( "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B","#D9EF8B", 
  #                 "#A6D96A" ,"#66BD63", "#1A9850", "#006837")
  if (BW){
    cc<-as.vector(paletteer_c("ggthemes::Classic Gray", levels_of_colors))
    # cols <- c("0" = "#fefefe", "1" = "#f1f1f1", "2" = "#e1e1e1", "3" = "#d9d9d9",
    #           "4" = "#dddddd", "5" = "#cccccc", "6" = "#a1a1a1", "7" = "#818181",
    #           "8" = "#666666", "9" = "#414141", "10" = "#000000")  
  }else{
    cc<-as.vector(paletteer_d(`"MetBrewer::Paquin"`,n = levels_of_colors,type="continuous"))
    if (palette==1) cc<-as.vector(paletteer_c("grDevices::RdYlGn", levels_of_colors))
    if (palette==2) cc<-as.vector(paletteer_d("MetBrewer::Benedictus",n = levels_of_colors,type="continuous"))
    if (palette==3) cc<-as.vector(paletteer_d("MetBrewer::Hiroshige",n = levels_of_colors,type="continuous"))
    if (palette==4) cc<-as.vector(paletteer_d("MetBrewer::Hiroshige",n = levels_of_colors,type="continuous",direction=-1))
    if (palette==5) cc<-as.vector(paletteer_d("MetBrewer::Tiepolo",n = levels_of_colors,type="continuous"))
    if (palette==6) cc<-as.vector(paletteer_d("PNWColors::Bay",n = levels_of_colors,direction=-1,type="continuous"))
    if (palette==7) cc<-as.vector(paletteer_c("scico::vik",n = levels_of_colors))
    if (palette==8) cc<-as.vector(paletteer_c("ggthemes::Temperature Diverging", n = levels_of_colors,direction=-1))
    
    
    names(cc)<-c(1:levels_of_colors)
    #to try "grDevices::Fall"  grDevices::Tropic (rev)  RdYlGn
    
    #DF<-DF %>% left_join(cc,by=c("dom"="val"))
    
    # cols <- c("0" = "#A50026", "1" = "#A50026", "2" = "#D73027", "3" = "#F46D43",
    #           "4" = "#FDAE61", "5" = "#FEE08B", "6" = "#D9EF8B", "7" = "#66BD63",
    #           "8" = "#1A9850", "9" = "#006837", "10" = "#003300")
  }
  p<-ggplot(DF,aes(x=labVar,y=freq))
  if (legend){
    p<-p+geom_bar(stat="identity",
                  aes(fill=fct_rev(as.factor(dom))),
                  width=1,alpha=alpha)
  }else{  
    p<-p+geom_bar(stat="identity",
                  aes(fill=fct_rev(as.factor(dom))),
                  width=1,
                  show.legend = F,alpha=alpha)
  }
  
  p<-p+ylim(c(-pupil,1))+
    geom_rect(aes(xmin=0.5,ymin = -pupil,xmax=max(IDc)+.5,ymax=0),fill=pupil.color)+
    scale_fill_manual(
      values = cc)+
    #    scale_fill_brewer(palette = "Spectral",direction=1)+
    #    scale_fill_brewer(palette = "RdYlGn",direction=1)+
    geom_hline(yintercept = 0.5,
               linetype="dashed",
               #color="white",
               #size=1,
               alpha=0.6)
  
  if(var_bar){
    tmp<-var_DF
    p<-p+geom_segment(data=tmp,aes(x=(xmin+xmax)*0.5,xend=(xmin+xmax)*0.5,
                                   y = 0.5-PROP*0.5,yend=0.5+PROP*0.5),
                      color="grey20",linewidth=1.5)
    # +
    #      geom_segment(data=tmp,aes(x=(xmin+xmax)*0.5,xend=(xmin+xmax)*0.5,
    #                                y = y2min+(1-PROP)*0.5,yend=y2max+(1-PROP)*0.5),
    #                color="grey30",linewidth=1.5)#+
    #   #   geom_segment(aes(x=0.5,xend=p+0.5,y=0,yend=0))+
    #   ylim(c(0,1))+theme_void()+
    #   coord_polar()
  }
  
  if(skewness_plo){
    bstat$x<-c(1:nrow(bstat))
    
    p<-p+geom_segment(inherit.aes = F,data=bstat,
                      aes(x=x+0.2,xend=x-0.2,y=coord,yend=coord),
                      color="black",#fill="orange",shape=8,
                      alpha=0.7,show.legend = F,
                      #stroke=1,
                      linewidth=1)
  }
  
  if(polar) p<-p+coord_polar()+theme_minimal()
  p<-p+theme_minimal()
  if(TITLE) p<-p+
    labs(#title = labels[selected]#,
      #subtitle = "subtitle",
      caption = labels[selected]
    )+
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5))
  if(notick){
    p<-p+
      theme(axis.text.x=element_blank())
  }
  p<-p+
    theme(axis.line=element_blank(),
          #        axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.border=element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(color="transparent",fill = "transparent"),
          panel.background = element_rect(color="transparent",fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = bg, color = NA),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0))
  return(p) 
}

Dist_HMAP_c<-function(Dmat, # original distributional data
                      lev.col=20,
                      hclu_R=F,
                      hclu_C=F,
                      exp_y_bott=0.2,
                      exp_x_right=0,
                      cut_labs_x=0,
                      cut_labs_y=0,
                      size_lab_x=3,
                      size_lab_y=3,
                      space_x=0.02,
                      space_y=0.02,
                      angle_x_lab=45,
                      angle_y_lab=0,
                      gg_fit_text=F,
                      interactive=F){
  n<-get.MatH.nrows(Dmat)
  p<-get.MatH.ncols(Dmat)
  row_labs<-get.MatH.rownames(Dmat)
  col_labs<-get.MatH.varnames(Dmat)
  
  
  res<-discretize_normalized_data(Dmat,n = lev.col)
  Dmat<-res$NEWMAT
  
  row_ord<-c(1:n)
  col_ord<-c(1:p)
  if(hclu_R){ 
    HC_r<-WH_hclust(Dmat)
    row_ord<-HC_r$order
    
  }
  if(hclu_C){
    
    HC_c<-hclust(as.dist(sqrt(1-round(WH.correlation(Dmat),8)^2)))
    col_ord<-HC_c$order
  }
  
  res<-discretize_data(Dmat,n = lev.col)# uses the minimum and the maximum of each variable
  DF<-squish_data(res$TIBN[,1:p],labels=res$TIBN[,(p+1)])
  
  #create rectangles
  
  DF<- DF %>% mutate(IDc_fin=as.numeric(factor(IDc,levels=col_ord)),
                     IDr_fin=as.numeric(factor(IDr,levels=row_ord)))
  DF<-DF %>% mutate(min_cdf=if_else(freq==cdf,0,cdf-freq))
  DF<-DF %>% mutate(Xmin=IDc_fin-0.5+min_cdf*(1-space_x),
                    Xmax=IDc_fin-0.5+cdf*(1-space_x),
                    Ymin=IDr_fin-(0.5-space_y),
                    Ymax=IDr_fin+(0.5-space_y))
  
  
  cc<-as.vector(paletteer_d(`"MetBrewer::Paquin"`,n = lev.col,type="continuous"))
  p0<-ggplot(DF)+geom_rect(aes(xmin=Xmin,xmax=Xmax,ymin=Ymin, ymax=Ymax, 
                               fill=as.factor(dom),
                               text=paste0("Id: ",labID,
                                           "</br></br>Var: ",labVar)),
                           show.legend = F)+scale_fill_manual(values=cc)
  
  labr_df<-data.frame(y=c(1:n),lab=row_labs[row_ord])
  
  
  if(gg_fit_text){  
    p0<-p0+geom_fit_text(data=labr_df,aes(xmin=rep(p+0.5,n),
                                          xmax=rep(p+0.5+min(1,0.1*p),n),
                                          ymin=y-0.5,
                                          ymax=y+0.5,label=lab))
  }else{
    p0<-p0+scale_y_continuous(breaks=seq(1,n,1),labels= labr_df$lab,position = "right",
                              expand = expansion(mult = c(0.01,0)))
      # geom_text(data=labr_df, aes(x=rep(p+0.5,n),y=y,label=lab),
      #           hjust=0,
      #           vjust=0.5,#nudge_x = 0.1,
      #           angle=angle_y_lab,
      #           size=size_lab_y)
    
  }
  #p0<-p0+scale_x_continuous(expand = expansion(mult = c(0.01, .2+(exp_x_right))))
 
  labc_df<-data.frame(x=c(1:p),lab=col_labs[col_ord])
  xv_j<-0.5
  xh_j<-1
  if(abs(angle_x_lab<10)){
    xv_j<-1
    xh_j<-0.5
  }
  
  if(gg_fit_text){ 
    p0<-p0+
       geom_fit_text(data=labc_df,aes(xmin=x-0.5,xmax=x+0.5,ymin=rep(-0.5-(p*0.1),p),ymax=0.5,
                                      x=x,y=rep(0,p),
                                      label=lab))
  }else{
    p0<-p0+scale_x_continuous(breaks=seq(1,p,1),labels= labc_df$lab,position = "bottom",
                              expand = expansion(mult = c(0.01,0)))
       # geom_text(data=labc_df, aes(y=rep(0,p),
       #                             x=x,label=lab),
       #                vjust=xv_j, hjust=xh_j,angle=angle_x_lab,size=size_lab_x)
    }
  
   # p0<-p0+ scale_y_continuous(expand = expansion(mult = c(exp_y_bott,0)))
  
  if(hclu_R){
    dendro_data<-ggdendro::dendro_data(HC_r)
    ss<-dendro_data$segments
    RY<-diff(range(c(ss$y,ss$yend)))
    minY<-min(c(ss$y,ss$yend))
    
    ss$y<-(ss$y-minY)/RY*0.25*(p+1)
    ss$yend<-(ss$yend-minY)/RY*0.25*(p+1)
    
    p0<-p0+
      geom_segment(data = ss, 
                   aes(y = x, x = -y+0.45, yend = xend, xend = -yend+0.45)
      )
  }
  
  if(hclu_C){
    dendro_data<-ggdendro::dendro_data(HC_c)
    ss<-dendro_data$segments
    RY<-diff(range(c(ss$y,ss$yend)))
    minY<-min(c(ss$y,ss$yend))
    
    ss$y<-(ss$y-minY)/RY*0.25*(n+1)
    ss$yend<-(ss$yend-minY)/RY*0.25*(n+1)
    
    p0<-p0+
      geom_segment(data = ss, 
                   aes(x = x, y = y+n+0.6, xend = xend, yend = yend+n+0.6)
      )
  }
  
  
  p0<-p0+theme(panel.border = element_blank(), 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               axis.line = element_blank(),
               panel.background = element_blank(),
               axis.ticks.x=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title=element_blank(),
               axis.text.x = element_text(size = size_lab_x, angle = angle_x_lab,
                                          #vjust=xv_j, 
                                          hjust=xh_j),
               axis.text.y = element_text(size = size_lab_y, angle = angle_y_lab,))
  
  if(interactive){
    p0<-ggplotly(p0,tooltip = c("text")) %>% plotly::layout(yaxis=list(side="right"),showlegend = FALSE)
  }
  return(p0)
}

ggplotly(Dist_HMAP_c(BLOOD,hclu_R = T,hclu_C = T,angle_x_lab = 0)+theme_bw(),
         tooltip = c("text")) %>% 
  plotly::layout(showlegend = FALSE)

# EXAMPLE  -----


run_that<-F
if (run_that){
  data<-China_Month[,13:24]
  data<-BLOOD
  
  # browser()
  nl<-50
  #res<-discretize_data(BLOOD,n = nl)
  res<-discretize_data(data,n = nl)
  res<-discretize_normalized_data(data,n = nl)
  data<-res$NEWMAT
  res2<-compute_LR_VARS_cont(data)
  tmp_DF<-draw_VARS(res2)
  #res<-discretize_data(China_Month[,13:24],n = nl,absolute = T)
  
  p<-GEI_plot2(res$TIBN[,1:get.MatH.ncols(data)],
               selected = 1,
               pupil = 0.3,
               pupil.color="white", labels = res$TIBN$ID_name,levels_of_colors = nl,
               palette=4)
  # 
  pp<-list()
  for(i in 1:nrow(res$TIBN)){
    # pp[[i]]<-GEI_plot3(res$TIBN[,1:12],
    #                    selected = i,
    #                    pupil = 0.4,
    #                    pupil.color="black",notick = T, 
    #                    labels = res$TIBN$ID_name,levels_of_colors = nl,palette=7)
    pp[[i]]<-GEI_plot3(res$TIBN[,1:get.MatH.ncols(data)],
                       selected = i,
                       pupil = 0.4,
                       pupil.color="grey30", notick = T,
                       labels = res$TIBN$ID_name,levels_of_colors = nl,
                       palette=7,
                       var_bar = T,var_DF = tmp_DF %>% filter(IDr==i),
                       skewness_plo = T)
    
  }
  library(patchwork)
  pp[[1]]+pp[[2]]+pp[[3]]+pp[[4]]+pp[[5]]+pp[[6]]+pp[[7]]+pp[[8]]+pp[[9]]+pp[[10]]+pp[[11]]+pp[[12]]+pp[[13]]+pp[[14]]
}
