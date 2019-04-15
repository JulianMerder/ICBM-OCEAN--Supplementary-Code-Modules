### precise join:: a sample junction algorithm
#required package data.table & dplyr
#'Data' represents output from MDL
#column mz = mass to charge ratio
#column I = Intensity
#column index = sample id
tolp<-0.5 #your tolerance
library(data.table)
library(dplyr)
Data$m1<-Data$mz #copy of mz column

mdlname<-names(Data)[3] # change column name
names(Data)[3]<-"MDL"

#sort data by I (descending)
Data<-dplyr::arrange(Data,desc(I))
#create row id
Data<-dplyr::mutate(Data,obsIdx = as.integer(1:n()))

#create matrix for corss table intensities
dg<-matrix(ncol=length(unique(Data$index)),nrow=nrow(Data))
colnames(dg)<-unique(Data$index)
dxg<-matrix(ncol=length(unique(Data$index)),nrow=nrow(Data))
colnames(dxg)<-unique(Data$index)
mmd<-Data$MDL
mmdres<-Data$ResPow
#index for already merged massesin upcoming for loop
indexx<-0


setDT(Data)
glo<-length(Data[,mz])
setkey(Data,obsIdx)

  for (i in 1:glo) {
    if(!any(indexx==i,na.rm=T)){
      
      sl <- data.table(dist=abs(((Data[,mz] - Data[i,mz])/Data[i,mz])*10^6),sample=Data[,index], int=Data[,I], indx=Data[,obsIdx],mzz=Data[,mz],m1=Data[,m1],MDL=Data[,MDL])
      sl<-sl[sl[ , .I[which.min(dist)], by = sample]$V1][dist<tolp,]
      
      dg[i,as.vector(sl[,sample])]<-as.vector(sl[,int])
      dxg[i,as.vector(sl[,sample])]<-as.vector(sl[,mzz])
      
      
      indexx<-c(indexx, as.integer(sl[,indx]))
      Data<-Data[!J(sl[,indx])]
      
    }
    
  }
  
  
  Data<-as.data.frame(Data)
  
  
  mzmean<-as.vector(sapply(1:nrow(dxg),function(x){weighted.mean(as.vector(dxg[x,]),as.vector((dg[x,])^1),na.rm=T) })) # possibility to scale weights by e.g. ^0.5
  sdmean<-apply(dxg,1,function(x) {sd(x,na.rm=T) } )
 
  colnames(dg) <- paste0("Sample", colnames(dg))
  colnames(dxg) <- paste0("Sample", colnames(dxg))
  
  Data<-data.frame(mz=mzmean,MDL=mmd,ResPow=mmdres)
  dd<-dplyr::bind_cols(Data,as.data.frame(dg))
  dd<-as.data.frame(dd)
  Data_precisejoin<-dd[rowSums(is.na(dd[,-c(1,2,3),drop=F])) != (ncol(dd)-3),]
  row.names(Data_precisejoin)<-NULL
  
  # if you like you can create a matrix for the original masses instead of sample intensities
  #dxg<-dplyr::bind_cols(Data,as.data.frame(dxg))
  #dxg<-as.data.frame(dd)
  #Data_precisejoin<-dxg[rowSums(is.na(dxg[,-c(1,2),drop=F])) != (ncol(dxg)-2),]
  #row.names(dxg)<-NULL
  
  
  names(Data_precisejoin)[2]<-mdlname
  