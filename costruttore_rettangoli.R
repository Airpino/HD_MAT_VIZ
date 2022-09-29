# n<-32
# x<-seq(from=0,to=1,length.out=n+1)
# x_ini<-x[1:n]
# x_fin<-x[2:(n+1)]
# xx<-seq(0,1,length.out=32)
# myc<-scales::seq_gradient_pal(low="red",high="blue")(xx)
# y_ini<-rep(0,n)
# y_fin<-rep(1,n)
# set.seed(1234)
# alp<-runif(32)
# 
# df<-data.frame(x_ini=x_ini,x_fin=x_fin,y_ini=y_ini,y_fin=y_fin,cols=myc,alp=alp)
# 
# ppp<-ggplot(df)+geom_rect(aes(xmin=x_ini,xmax=x_fin,ymin=y_ini,ymax=y_fin,fill=cols,alpha=alp),show.legend = F)+
#   scale_fill_identity()+theme_bw()


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


runme<-F
if(runme){
df3<-prep_for_plot(B,v=sels,poi = 50,
                                       pp = 50,
                                       n = 32,#kernel density parameter
                                       order_of_rows = order_of_rows,
                                       order_of_cols = order_of_cols)
DDF_to<-create_new_df_to_plot_M(df3)
# ggplot(DDF)+geom_rect(aes(xmin=x_ini,xmax=x_fin,ymin=y_ini,ymax=y_fin,fill=cols,alpha=alp),show.legend = F)+
#   scale_fill_identity()+theme_few()+ylim(c(0,78))

ggplotly(ggplot(DDF_to)+geom_rect(aes(xmin=x_ini,xmax=x_fin,ymin=y_ini,ymax=y_fin,fill=cols,alpha=alp),show.legend = F)+
   scale_fill_identity()+facet_grid(col=vars(VAR_name),scales = "free_x")+theme_few()+
     theme(legend.position = "none"))

ggplotly(ggplot(DDF_to)+geom_point(aes(x=x_ini,y=y_ini,color=cols,alpha=alp),show.legend = F)+
           scale_fill_identity()+facet_grid(col=vars(VAR_name),scales = "free_x")+
           theme_few()
           )
#
# # library(plotly)
# # fig <- plot_ly(DDF_to, x = ~x_ini, y = ~y_ini, name = "myplot")
# # 
# # # add shapes to the layout
# # fig <- layout(fig, title = 'rectangles',
# #               shapes = list(
# #                 list(type = "rect",
# #                      fillcolor = "blue", line = list(color = "blue"), opacity = 0.3,
# #                      x0 = "1980-01-01", x1 = "1985-01-01", xref = "x",
# #                      y0 = 4, y1 = 12.5, yref = "y"),
# #                 list(type = "rect",
# #                      fillcolor = "blue", line = list(color = "blue"), opacity = 0.2,
# #                      x0 = "2000-01-01", x1 = "2005-01-01", xref = "x",
# #                      y0 = 4, y1 = 12.5, yref = "y")))
# # 
# # fig
# 
# 
ppp<-ggplot(DDF_to)+
  geom_rect(aes(xmin=x_ini,xmax=x_fin,
                ymin=y_ini,ymax=y_fin,
                fill=cols,alpha=alp),show.legend = F)+
  scale_fill_identity()+
  facet_grid(col=vars(VAR_name),scales = "free_x")+
  theme_few()+
  theme(#plot.margin=margin(r=2.5,unit = "cm"),
        panel.spacing.y=unit(0, "lines"),
        #strip.text =element_text(margin = margin(r=20)),
        #strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
        strip.text.x = element_text(size = 7, angle = 0),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        legend.position='none')
tmp1<-DDF_to %>% group_by(ID_name) %>% summarise(y_ini=min(y_ini),y_fin=min(y_fin))
myth<-theme(#plot.margin=margin(r=2.5,unit = "cm"),
  #panel.spacing.y=unit(0, "lines"),
  #strip.text =element_text(margin = margin(r=20)),
  #strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
  #strip.text.x = element_text(size = 7, angle = 0),
  axis.text.y = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.border = element_blank(),
  legend.position='none')
pplabs<-ggplot(tmp1) +
  geom_text(aes(x=1,y=(y_fin+y_ini)*0.5,label=ID_name),size=2)+
  theme_few()+myth
gridExtra::grid.arrange(ppp,
                        pplabs,#+theme(strip.text.x = element_blank()),
                        ncol=2, widths = c(5,1))
show(plotly::subplot(
  nrows =2,
  plotly::plotly_empty(type = "scatter",mode   = 'markers'),
  plotly::ggplotly(pv),
  plotly::plotly_empty(type = "scatter",mode   = 'markers'),
  plotly::ggplotly(p),
  plotly::ggplotly(ppp+scale_y_continuous(limits=c(0,78),expand=c(0,0))), 
  plotly::ggplotly(pplabs+scale_y_continuous(limits=c(0,78),expand=c(0,0))),
  heights = c(0.2,0.8),
  widths = c(0.15,0.75, 0.1)))
}

library(plotly)
DDF_to%>%
  group_by(VAR_name) %>% 
  do(p=plot_ly(., x = ~x_ini, y = ~y_ini, 
               color = ~cols, 
               colors= ~tolower(cols), 
               group=x,
               opacity= ~alp,
               type = "scatter")) %>%
  subplot(nrows = 1, shareX = FALSE, shareY = TRUE)
