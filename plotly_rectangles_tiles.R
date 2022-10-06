library(plotly)

library(ggdendro)
library(HistDAWass)
source("MATH_TO_TIBBLE.R")
data<-HistDAWass:::c_Prepare(China_Month,simplify=T,qua=50,standardize=F)$x#China_Seas
B<-MATH2tibble(data)
sels=c(13:24)#c(5:8,17:20,29:32)
### set layout
per_left<-20 #perc for left panel header
per_right<-80 #perc for right panel header
per_up <-10  #perc for top panel header
per_bott<-90  #perc for bottom panel header

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

### prepare data for plot
df4<-prep_for_plot_freq(B,v=sels,
                        n = 64,#kernel density parameter
                        order_of_rows = order_of_rows,
                        order_of_cols = order_of_cols)
df5<-df4 %>% mutate(x_n=(NX+as.numeric(ID_VAR)-min(NX+as.numeric(ID_VAR)))/diff(range(NX+as.numeric(ID_VAR)))*per_right+per_left,
                    y_n=(as.numeric(ID_name)-min(as.numeric(ID_name)))/diff(range(as.numeric(ID_name)))*per_bott,
                    fc=y)
df6<-df5 %>% mutate(fc=if_else(round(fc,4)==0,NA_real_,round(fc,4)))

segs_of_vars<-ddata_v$segments
rx<-range(c(unlist(segs_of_vars$x),unlist(segs_of_vars$xend)))
rx[1]<-rx[1]-0.5
rx[2]<-rx[2]+0.5
ry<-range(c(unlist(segs_of_vars$y),unlist(segs_of_vars$yend)))
segs_of_vars<-segs_of_vars%>% mutate(x=(x-rx[1])/diff(rx)*per_right+per_left,
                                     xend=(xend-rx[1])/diff(rx)*per_right+per_left,
                                     y=(y-min(ry[1]))/diff(ry)*per_up+per_bott+6,
                                     yend=(yend-min(ry[1]))/diff(ry)*per_up+per_bott+6)
#montiamo il dendrogramma a sinistra
segs_of_inds<-ddata$segments
rx<-range(c(unlist(segs_of_inds$x),unlist(segs_of_inds$xend)))
rx2<-range(c(unlist(segs_of_inds$x),unlist(segs_of_inds$xend)))
rx[1]<-rx[1]-0.5
rx[2]<-rx[2]+0.5
### aggiusta qui!!!!
segs_of_inds2<-segs_of_inds%>% mutate(x=(x-rx[1])/(diff(rx))*(per_bott),
                                      xend=(xend-rx[1])/(diff(rx))*per_bott,
                                      y=(1-(y+0.5)/(max(y)+1))*per_left,
                                      yend=(1-(yend+0.5)/(max(yend)+1))*per_left)
segs_of_inds3<-segs_of_inds%>% mutate(x=(x-rx2[1])/(diff(rx2))*(per_bott),
                                      xend=(xend-rx2[1])/(diff(rx2))*per_bott,
                                      y=(1-(y+0.5)/(max(y)+1))*per_left,
                                      yend=(1-(yend+0.5)/(max(yend)+1))*per_left)
## attacchiamo le etichette righe e colonne

ID_n<-data.frame(names=levels(df4$ID_name))
ID_v<-data.frame(names=levels(df4$ID_VAR))
ID_n <-ID_n %>% mutate(x=rep(102,nrow(.)),y=(c(1:(nrow(.)))-1)/(nrow(.)-1)*per_bott)
ID_v <-ID_v %>% mutate(x=(c(1:nrow(.))-0.5)/nrow(.)*per_right+per_left,y=per_bott+2)

nrows<-length(levels(df4$ID_name))
vars<-length(levels(df4$ID_VAR))



ray_x<-min(round(diff(sort(unique(df6$x_n))),3))/2
ray_y<-min(round(diff(sort(unique(df6$y_n))),3))/2

df6 <- df6 %>% mutate(x_min=x_n-ray_x,x_max=x_n+ray_x,
                      y_min=y_n-ray_y,y_max=y_n+ray_y)

var_limits<-df4 %>% group_by(ID_VAR) %>% summarize(min=min(x),maxs=max(x)) %>%
  mutate(ord=c(1:nrow(.)),
         new_m=((ord-1)+0.1)/(nrow(.))*per_right+per_left,
         new_max=(ord-0.1)/(nrow(.))*per_right+per_left, 
         new_mid=((ord-1)+0.5)/(nrow(.))*per_right+per_left,
         mid_f=prettyNum((maxs+min)*0.5,digits=4),
         mins_f=prettyNum(min, digits=4),
         maxs_f=prettyNum(maxs, digits=4))

move<-per_bott/((nrow(ID_n)-1)*2)

fig <- plot_ly(
  x = df6$x_n, y =df6$y_n,
  z = df6$fc,
  type = "heatmap",
  hoverinfo='none',
  colors = colorRamp(c("yellow", "red")),
) %>%   add_markers(
  inherit = F,
  x = df6$x_n, y = df6$y_n,
  #data = df.risk,
  showlegend = F,
  text = paste0("V: ",as.character(df6$ID_VAR),"<br>",
                "O: ",as.character(df6$ID_name),"<br>",
                "val: ",prettyNum(df6$x,digits=4),"<br>",
                "r_f: ",prettyNum(df6$y_orig,digits=3),"<br>"),
  color = I("transparent"),# or whatever color you prefer
  hovertemplate = paste0("%{text}<br>")
) %>%  add_segments(
  #fig,
                    x = segs_of_vars$x,
                    y = segs_of_vars$y+move,
                    z =NULL,
                    xend = segs_of_vars$xend,
                    yend = segs_of_vars$yend+move,
                    data = NULL,
                    inherit = TRUE
) %>% add_segments(
  #fig,
                   x = segs_of_inds3$y,
                   y = segs_of_inds3$x,
                   z =NULL,
                   xend = segs_of_inds3$yend,
                   yend = segs_of_inds3$xend,
                   data = NULL,
                   inherit = TRUE
)%>% layout(annotations=list(
  #add_annotations(
  x = ID_n$x,
  y = ID_n$y,
  z = NULL,
  text = ID_n$names,
  font=list(size=10),
  xref = "x",
  yref = "y",
  xanchor='left',
  showarrow = FALSE)
)%>% add_annotations(x = ID_v$x,
                    y = ID_v$y+move,
                    z=NULL,
                    text = ID_v$names,
                    font=list(size=8),
                    xref = "x",
                    yref = "y",
                    showarrow = FALSE
) %>% add_segments(
  #fig,
             x = var_limits$new_m,
             y = -(1+move),
             z =NULL,
             xend = var_limits$new_max,
             yend = -(1+move),
             data = NULL,
             inherit = TRUE
)%>% add_annotations(x = var_limits$new_m,
                     y = -(5+move),
                     z=NULL,
                     text = var_limits$mins_f,
                     xref = "x",
                     yref = "y",
                     showarrow = FALSE,
                     align="left",
                     font=list(size=8),
                     textangle=90,
                     valign='bottom'
)%>% add_annotations(x = var_limits$new_max,
                     y = -(5+move),
                     z=NULL,
                     text = var_limits$maxs_f,
                     xref = "x",
                     yref = "y",
                     showarrow = FALSE,
                     align='left',
                     font=list(size=8),
                     textangle=90,
                     valign='bottom')

show(fig %>%  
   layout(yaxis = list(showgrid = FALSE,
                              showticklabels = FALSE,
                              visible=FALSE),
                 xaxis = list(showgrid = FALSE,
                              showticklabels = FALSE,
                              visible=FALSE),
                 showlegend = FALSE) %>% hide_colorbar()
)
# plot_ly(
#   x = df6$x_n, y =df6$y_n,
#   type = "scatter",mode="markers",
#   #  z = df6$fc, type = "heatmap",
#   color=df6$x,
#   colors = colorRamp(c("yellow", "red")),#,
#   opacity=df6$fc
# )
##################DOTS WITH TRASPARENCY#################################################
tmp1<-df6 %>% select(x_n,y_n,NX,fc) %>% 
  mutate(NX_2=ceiling(NX*31)+1,
         COL=colorRampPalette(c("yellow", "red"),space="rgb")(32)[NX_2]) %>% mutate(
           opac=if_else(is.na(fc),0,round(sqrt(fc),2)),
           red=col2rgb(COL)[1,],
           green=col2rgb(COL)[2,],
           blue=col2rgb(COL)[3,],
           fincol=paste0("rgba(",red,", ",green,", ",blue,", ",opac,")"))
f2<-plot_ly(type = 'scatter', mode = 'markers',x=tmp1$x_n,y=tmp1$y_n,
            marker=list(
              symbol="square",
              color=tmp1$fincol,
              #opacity=0.5,
              line = list(width=0))
)%>% layout(plot_bgcolor='rgb(0, 0, 0)')

fig2<-plot_ly(type = 'scatter', mode = 'markers',x=c(1,2),y=c(4,3)) %>% 
  layout(
  annotations=list(x = c(1,2),
                  y = c(4,3),
                  text = c("Ci+ao", "come-stai"),
           showarrow=FALSE,
           xanchor='left'))
fig2


