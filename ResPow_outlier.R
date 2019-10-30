# dataset: your spectrum
# column_mass: which is your mass column? (number)
# column_respow: which is your Resolution Power column (number)
# method: outlier test based on kernel density estimation or mirroring of negative residuals? numbers: 1= kernel density, 2 = min. residuals


ResPow_outlier<-function(dataset,column_mass=1,column_respow=4,method=1){
  package <- c("quantreg","ggplot2")
  install <- package[!(package %in% installed.packages()[,"Package"])]
  if(length(install)) install.packages(install)
  
  column_mass<-as.integer(column_mass)
  column_respow<-as.integer(column_respow)
  library(quantreg)
  library(ggplot2)
  if(!is.integer(column_mass) | !is.integer(column_respow)) stop("no columns selected")
  
  name_mass<-names(dataset)[column_mass]
  name_respow<-names(dataset)[column_respow]
  names(dataset)[column_mass]<-"m.z"
  names(dataset)[column_respow]<-"ResPow"
  
  
  d<-rq(log(ResPow)~log(m.z),data=dataset,method="fn")
  
  
  
  if(method==1){
  j<-density(d$residuals)
  e<-j$x[c(-1)][j$x[c(-1)]>0][which(diff(j$y)[j$x[c(-1)]>0]>0)[1]]  #first derivative of kernel sensity estimation, where it becomes positive again is set as the outlier treshold
  }
  else{
  e<-min(d$residuals)  
  }
  
  dataset$residual_color<-NA
  dataset$residual_color[d$residuals>abs(e)]<-"removed"
  dataset$residual_color[d$residuals<=abs(e)]<-"retained"
  PLOT<-ggplot(dataset,aes_(x=dataset$m.z,y=dataset$ResPow,color=dataset$residual_color))+geom_point()+scale_color_manual(values=c("removed"="red","retained"="lightblue"))+labs(color="peak",x=name_mass,y=name_respow)
  
  names(dataset)[column_mass]<-name_mass
  names(dataset)[column_respow]<-name_respow
  
  dataset<-dataset[d$residuals>abs(e),]
  dataset$residual_color<-NULL
  return(list(dataset=dataset,plot=PLOT))
  
}

# execute function, where D is your dataset
D<-read.csv("mydata.csv")
KK<-ResPow_outlier(D)

# function returns a list with a plot and the new dataset without outliers

#plot can take some time to pop up because of the huge number of peaks
KK$plot

# the new dataset without outliers  
my_new_data<-KK$dataset



