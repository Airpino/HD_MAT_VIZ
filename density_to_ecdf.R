#### density to ECDF ----

dd<-BLOOD@M[1,1][[1]]
quants<-HistDAWass:::compQ_vect(dd,vp=c(0:100)/100)
den<-density(quants,from = min(dd@x),to=max(dd@x))

area<-sum((den$y[1:(length(den$y)-1)]+den$y[2:length(den$y)])/2*diff(den$x))
cfreq<-cumsum(c(0,(den$y[1:(length(den$y)-1)]+den$y[2:length(den$y)])/2*diff(den$x)))
cfreq<-cfreq/max(cfreq)
plot(cfreq)

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



## generate plot for reordered rows
plo<-plots_strips_distr_2(B,v=sels,poi = 50,
                          pp = 50,
                          n = 32,#kernel density parameter
                          order_of_rows = order_of_rows,
                          order_of_cols = order_of_cols)

## generate the dendrogram of rows
nels<-nrow(ddata$labels)
p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0, 0))+
  scale_x_reverse(expand = c(0, 0),limits = c((nels+0.5),0.5))+
  theme_dendro()

## generate the dendrogram of rows
nva<-nrow(ddata_v$labels)
pv <- ggplot(segment(ddata_v)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  scale_x_continuous(expand = c(0, 0),limits = c(0.5,(nva+0.5)))+
  theme_dendro()
#pv


show(plotly::subplot(
  nrows =2,
  plotly::plotly_empty(type = "scatter",mode   = 'markers'),
  plotly::ggplotly(pv),
  plotly::ggplotly(p), 
  plotly::ggplotly(plo),
  heights = c(0.25,0.75),
  widths = c(0.25, 0.75)))

runme<-F
if(runme){
## df2<-mul_density(B,v = v,n=n,poi=poi,pp = pp,NORM=NORM)
df2<-mul_density(B,v = sels,n=50,poi=50,pp = 32,NORM=T)
ggplot(df2 %>% filter(VAR=="5") ,
       aes(x = x, y = 1))+
  geom_tile(colour="#FFFFFF00",aes(alpha=y*0.001),fill = "red",show.legend = T)+
  scale_alpha(range = c(0,1))

ggplot(df2 ,
       aes(x = x, y = 1))+
  geom_tile(aes(fill = y), show.legend = T)+
  scale_fill_gradient(low="yellow",high="red")+
  facet_grid(rows = vars(ID_name),cols = vars(ID_VAR),scales = "free_x")+
  theme_few()+
  theme(panel.spacing.y=unit(0, "lines"))

col_strip <- brewer.pal(9, "Reds")
ggplot(df2 ,
              aes(x = x, y = 1, fill = y))+
  geom_raster(show.legend = T)+#show.legend = F)+
  scale_fill_gradient(low="white",high="red")+
  facet_grid(rows = vars(ID_name),
             cols =vars(ID_VAR),
             scales = "free_x"
  )
+
  theme_few() +
  theme(plot.margin=margin(r=2.5,unit = "cm"),
        panel.spacing.y=unit(0, "lines"),
        #strip.text =element_text(margin = margin(r=20)),
        strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
        strip.text.x = element_text(size = 7, angle = 20),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        legend.position='none')
}