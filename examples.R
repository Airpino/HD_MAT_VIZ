### Questo file contiene esempi di rappresentazioni
### 
### carico le funzioni che servono

source("MATH_TO_TIBBLE.R")


######## Complete application-------
library(plotly)
library(ggdendro)
library(HistDAWass)
source("MATH_TO_TIBBLE.R")
data<-OzoneFull
B<-MATH2tibble(data)
sels=c(1:4)#,17:20,29:32)
#plo<-plots_strips_distr_2(B,v=sels,poi = 50,pp = 50,n = 32)
#ggplotly(plo)#funziona benissimo!!

## Row Hclustering

## Hclustering of rows
HCLU<-WH_hclust(data[,sels],standardize = T) #do Hclust for dd

order_of_rows<-HCLU$order # the order of rows
dhc <- as.dendrogram(HCLU)
ddata <- dendro_data(dhc, type = "rectangle")

## Columns Hclust
c_mat<-WH.correlation(data[,sels])
vhclu<-hclust(as.dist(round(1-c_mat,5)))#using correlation distance
order_of_cols<-vhclu$order
dh_v <- as.dendrogram(vhclu)
ddata_v <- dendro_data(dh_v, type = "rectangle")


## generate the dendrogram of rows
nels<-nrow(ddata$labels)
p_r <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0, 0))+
  scale_x_reverse(expand = c(0, 0),limits = c((nels+0.5),0.5))+
  theme_dendro()

## generate the dendrogram of vars
nva<-nrow(ddata_v$labels)
p_c <- ggplot(segment(ddata_v)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  scale_x_continuous(expand = c(0, 0),limits = c(0.5,(nva+0.5)))+
  theme_dendro()





runme<-F
if(runme){
df3<-prep_for_plot(B,v=sels,poi = 50,
                                       pp = 50,
                                       n = 32,#kernel density parameter
                                       order_of_rows = order_of_rows,
                                       order_of_cols = order_of_cols)
df4<-prep_for_plot_freq(B,v=sels,
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
# show(plotly::subplot(
#   nrows =2,
#   plotly::plotly_empty(type = "scatter",mode   = 'markers'),
#   plotly::ggplotly(pv),
#   plotly::plotly_empty(type = "scatter",mode   = 'markers'),
#   plotly::ggplotly(p),
#   plotly::ggplotly(ppp+scale_y_continuous(limits=c(0,78),expand=c(0,0))), 
#   plotly::ggplotly(pplabs+scale_y_continuous(limits=c(0,78),expand=c(0,0))),
#   heights = c(0.2,0.8),
#   widths = c(0.15,0.75, 0.1)))
# 
}

# library(plotly)
# DDF_to%>%
#   group_by(VAR_name) %>% 
#   do(p=plot_ly(., x = ~x_ini, y = ~y_ini, 
#                color = ~cols, 
#                colors= ~tolower(cols), 
#                group=x,
#                opacity= ~alp,
#                type = "scatter")) %>%
#   subplot(nrows = 1, shareX = FALSE, shareY = TRUE)
