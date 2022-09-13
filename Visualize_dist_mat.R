library(HistDAWass)
data<-China_Month
#get min and max for each column
mins<-get.MatH.stats(data,stat="min")$mat
cmins<-sapply(data.frame(mins), min) 
maxs<-get.MatH.stats(data,stat="max")$mat
cmaxs<-sapply(data.frame(maxs), max)
data_st<-data
for(i in 1:nrow(data@M)){
  for (j in 1:ncol(data@M)){
    data_st@M[i,j][[1]]@x<-(data@M[i,j][[1]]@x-cmins[j])/(cmaxs[j]-cmins[j])
    data_st@M[i,j][[1]]@m<-(data@M[i,j][[1]]@m-cmins[j])/(cmaxs[j]-cmins[j])
    data_st@M[i,j][[1]]@s<-(data@M[i,j][[1]]@s)/(cmaxs[j]-cmins[j])
  }
}

#divide into 10 bins
# nb<- 10
# bins<-c(0:nb)/nb
# 
# # calculate empirical frequency
# for(i in 1:nrow(data@M)){
#   for (j in 1:ncol(data@M)){
#     freq<-diff(compP(data_st@M[i,j][[1]]))
#   }
# }

#only temp
CMAT<-WH.correlation(data[,13:24])
dm<-HistDAWass::WH_MAT_DIST(data[,13:24])#only temp