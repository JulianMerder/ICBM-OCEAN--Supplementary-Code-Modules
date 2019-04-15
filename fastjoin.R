### fast join:: a sample junction algorithm
#required package dplyr & tidyr
#'Data' represents output from MDL
#column mz = mass to charge ratio
#column I = Intensity
#column index = sample id
tolp<-0.5 #your tolerance
mdlname<-names(Data)[3]
names(Data)[3]<-"MDL"
library(dplyr)
library(tidyr)
# Sort data by mass
Data_fastjoin<-dplyr::arrange(Data,mz) %>%
  # Start with the smallest mass and check if the next mass is the tolerance range defined by 'tolp'. If Yes group togetherinto a mass cluster, if not start a new group (mass cluster). 
  dplyr::mutate(mz1=mz,group = cumsum( abs(((mz - lag(mz, default = mz[1]) )/lag(mz, default = mz[1]))*10^6) > tolp)) %>%
  # Group by found mass cluster and check for duplicate samples within a mass cluster
  dplyr::group_by(group) %>%
  dplyr::mutate(max=(mz[which.max(I)]-mz)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group, index) %>% dplyr::arrange(max)%>% dplyr::mutate(gjk=as.character(ifelse(dplyr::row_number()==1,paste0(group,"x",dplyr::row_number()),paste0(group,"x",index,"x",dplyr::row_number()))))%>%            #dplyr::summarize(mz=(mz[which.min(max)]),I=(I[which.min(max)]),MDL=(MDL[which.min(max)]),ResPow=(ResPow[which.min(max)]),mz1=mz1[which.min(max)],refe=refe[which.min(max)],m1=m1[which.min(max)])%>%
  dplyr::ungroup() %>% 
  dplyr::select(c(-max,-group))%>%
  dplyr::group_by(gjk) %>%
  dplyr::mutate(
  # calculate weighted average mass within each mass group 
  mz = weighted.mean(mz, sqrt(I), na.rm=T),MDL=mean(MDL), ResPow=mean(ResPow),SE=sd(mz1)/mz*1e6,mz1=mean(mz1)) %>%  
  
  # Take out the  group assignment and columns not needed
  dplyr::ungroup() %>%
  dplyr::select(c(-mz1,-gjk)) %>%
  # Add the "B" leader to the index
  dplyr::mutate(index = paste0("Sample", index)) %>%
  # Put intensity values of each sample into separate columns
  spread(index, I) %>% as.data.frame()

names(Data_fastjoin)[2]<-mdlname


# in case order of samples is important to preserve:
#library(readr)
#i2 <- grep("ample", names(Data_fastjoin))
#i3 <- names(Data_fastjoin)[i2]
#i4 <- readr::parse_number(i3)
#names(Data_fastjoin)[i2] <- i3[order(i4)]
#Data_fastjoin[i2] <- Data_fastjoin[i2][order(i4)]


