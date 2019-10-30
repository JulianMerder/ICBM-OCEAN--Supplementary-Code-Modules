### function to calculate isotope ratio deviances

isotope_deviance<-function(intensity_parent,intensity_child,number_of_atoms,q=0.0108,p=0.9892,number_of_isotopes=1,group=as.character(" "),color_by=as.character("black"),add_points=T,ptitle="",xlabel="",ylabel=expression(paste(delta^{13}, "C (\u2030)"))){ # ylabel: for e.g. O18 change 13 to 18 and "C (\u2030)" to "O (\u2030)"
  #intensity_parent: intensity of the isotope free peak (e.g. pure C-12); can be a vector
  #intensity_child: intensity of the isotope peak (e.g. formula with one C-13); can be a vector
  #number_of_atoms: the number of total atoms of the element of interest (e.g. for Delta C-13 in  C6H12O6 use 6); can be a vector
  #number_of_isotopes: the number of the rare isotope inside the formula; for one C-13 use 1, for two C-13 use 2
  #q natural abundance or standard of the rare isotope e.g. C-13 . For Iron and elements with more than two stable isotopes use dmultinom() instead
  #p natural abundance or standard of the common isotope e.g. C-12; this is not neccesarily 1-q e.g. for Iron. For Iron use dmultinom() instead
  
  
###calculate deviance###
  measured_intensity_ratio<-intensity_child/intensity_parent
  expected_intensity_ratio<-dbinom(number_of_isotopes,number_of_atoms,q)/dbinom(number_of_atoms,number_of_atoms,p)
  
  rdeviance<-((measured_intensity_ratio/expected_intensity_ratio)-1)*1000

  
  
###create plots:###
if (!require("ggplot2")) {install.packages("ggplot2"); library(ggplot2)}


 if(sum(!is.na(rdeviance))>2){
   
  if(isTRUE(add_points)){
    if(is.numeric(color_by)){
    P<-ggplot(data.frame(x=group,y=rdeviance),aes(x=x,y=y))+geom_violin(scale="area",trim=FALSE,adjust=1.2,position=position_dodge(width = 0.7))+geom_boxplot(aes(group=x),outlier.shape=NA,position=position_dodge(width = 0.7), fill="white",width=0.1)+labs(x=xlabel,y=ylabel)+theme_classic()+theme(text = element_text(size=16))+ggtitle(ptitle)+geom_jitter(data=data.frame(x=group,y=rdeviance,color_=color_by),aes(x=x,y=y,color=color_),height=0)+theme(legend.title=element_blank())+scale_color_gradient(low="grey",high="red") 
    } else{
    
      P<-ggplot(data.frame(x=group,y=rdeviance),aes(x=x,y=y))+geom_violin(scale="area",trim=FALSE,adjust=1.2,position=position_dodge(width = 0.7))+geom_boxplot(aes(group=x),outlier.shape=NA,position=position_dodge(width = 0.7), fill="white",width=0.1)+labs(x=xlabel,y=ylabel)+theme_classic()+theme(text = element_text(size=16))+ggtitle(ptitle)+geom_jitter(data=data.frame(x=group,y=rdeviance,color_=color_by),aes(x=x,y=y,color=color_),height=0)+theme(legend.title=element_blank()) 
      
    }
  } else{
  
    P<-ggplot(data.frame(x=group,y=rdeviance),aes(x=x,y=y))+geom_violin(scale="area",trim=FALSE,adjust=1.2,position=position_dodge(width = 0.7))+geom_boxplot(aes(group=x),outlier.shape=NA,position=position_dodge(width = 0.7), fill="white",width=0.1)+labs(x=xlabel,y=ylabel)+theme_classic()+theme(text = element_text(size=16))+ggtitle(ptitle)+geom_jitter(data=data.frame(x=group,y=rdeviance,color_=color_by),aes(x=x,y=y,color=color_),height=0)+theme(legend.title=element_blank()) 
    
  }
 }else{
  P<-ggplot(data.frame(x=1:10,y=1:10),aes(x,y))+annotate("text",x=1,y=2,label="Dear Geochemist, you have not enough ratios to plot a density",size=5)+labs(title="",x ="", y = "")+ theme(panel.background = element_blank())

 }
  
  return(list(deviance=rdeviance,plot=P))
 }

#example: measure Delta C-13 of two isotope pairs with a natural abundance of q=1.08% and for this case p=1-q. We consider the case for a single C13 in the molecule
 #q=1.08/100
 #p=1-q
 #A<-data.frame(Isotopefree_formula=c("C_6 H_12 O_6","C_10 H_22 O_1","C_5 H_12 O_4"),Isotope_containing_formula=c("C_5 H_12 O_6 C13_1","C_9 H_22 O_1 C13_1","C_4 H_12 O_4 C13_1"),intensity_parent=c(1000,2500,1500),intensity_child=c(63,270,80),C=c(6,10,5),H=c(12,22,12),O=c(6,1,4),color_by=as.character(c("Group 1","Group 2","Group 3")))
 #A$Delta_C13<-isotope_deviance(A$intensity_parent,A$intensity_child,number_of_atoms=A$C,q,p,number_of_isotopes=1,color_by=A$color_by)$deviance
 #isotope_deviance(A$intensity_parent,A$intensity_child,number_of_atoms=A$C,q,p,number_of_isotopes=1,color_by=A$color_by)$plot
 #isotope_deviance(A$intensity_parent,A$intensity_child,number_of_atoms=A$C,q,p,number_of_isotopes=1,color_by=c(1,2,3))$plot
 #C<-A[1:2,]
 #isotope_deviance(C$intensity_parent,C$intensity_child,number_of_atoms=C$C,q,p,number_of_isotopes=1,color_by=c(1,2))$plot
 #B<-rbind(A,A)
 #B$group<-c(rep("A",nrow(A)),rep("B",nrow(A)))
 #isotope_deviance(B$intensity_parent,B$intensity_child,number_of_atoms=B$C,q,p,number_of_isotopes=1,group=B$group,color_by=B$color_by,add_points=T)$plot
 #D<-isotope_deviance(A$intensity_parent,A$intensity_child,number_of_atoms=A$C,q,p,number_of_isotopes=1,color_by=A$color_by)
 #D$deviance
 #D$plot
 #A<-data.frame(Isotopefree_formula=c("C_6 H_12 O_6","C_10 H_22 O_1","C_5 H_12 O_4"),Isotope_containing_formula=c("C_6 H_12 O_5 O18_1","C_10 H_22 O18_1","C_5 H_12 O_3 O18_1"),intensity_parent=c(1000,2500,1500),intensity_child=c(12.07,5,12),C=c(6,10,5),H=c(12,22,12),O=c(6,1,4),color_by=as.character(c("Group 1","Group 2","Group 3")))
 #q=0.2/100
 #p=99.76/100
 #isotope_deviance(A$intensity_parent,A$intensity_child,number_of_atoms=A$O,q,p,number_of_isotopes=1,color_by=c(1.4,2.8,3),ylabel=expression(paste(delta^{18}, "O (\u2030)")))$plot
 #isotope_deviance(A$intensity_parent,A$intensity_child,number_of_atoms=A$O,q,p,number_of_isotopes=1,color_by=c(1.4,2.8,3))$deviance
 
