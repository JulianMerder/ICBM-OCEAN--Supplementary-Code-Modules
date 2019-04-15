# create homologous series network : here, not from Kendricks mass defekt, but from molecular formulas directly
# needed packages: igraph and vis network for latter visualization of the edgelist

get_network<-function(i,r,B,indexf,homcheck){
  
  
  
  h<-c(NA,NA,NA,NA,NA,NA,i)
  if(i<r){
    ee<-sum(B[i,1]==B[,1])+sum(B[i,1]==(B[,1]-1))
    zu<-apply(B[i,-c(1,2,12)],1,function(x) paste(x,collapse=""))
    zu1<-apply(B[(i+1):(i+min(ee,r-i)),-c(1,2,12)],1,function(x) paste(x,collapse="")) 
    #CH2
    if (as.character(1) %in% homcheck){
      if (any(B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]-1) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]-2) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12]) & zu==zu1)) 
      {
        h[1]<-as.character(indexf[(i+1):(i+min(ee,r-i))][B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]-1) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]-2) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12]) & zu==zu1])[1]
      }
    }
    #CO2
    if (as.character(2) %in% homcheck){
      if (any(B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]-1) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12]-2)& zu==zu1)) 
      {
        h[2]<-as.character(indexf[(i+1):(i+min(ee,r-i))][B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]-1) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12]-2)& zu==zu1])[1]
      }
    }
    #H2
    if (as.character(3) %in% homcheck){
      if (any(B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]-2) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12]) & zu==zu1)) 
      {
        h[3]<-as.character(indexf[(i+1):(i+min(ee,r-i))][B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]-2) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12])& zu==zu1])[1]
      }
    }
    #H2O
    if (as.character(4) %in% homcheck){
      if (any(B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]-2) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12]-1)& zu==zu1)) 
      {
        h[4]<-as.character(indexf[(i+1):(i+min(ee,r-i))][B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]-2) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12]-1)& zu==zu1])[1]
      }
    }
    #O
    if (as.character(5) %in% homcheck){
      if (any(B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12]-1)& zu==zu1)) 
      {
        h[5]<-as.character(indexf[(i+1):(i+min(ee,r-i))][B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12]-1)& zu==zu1])[1]
      }
    }
    #isotope
    if (any(B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12])& zu==zu1)) 
    {
      h[6]<-as.character(indexf[(i+1):(i+min(ee,r-i))][B[i,1]==(B[(i+1):(i+min(ee,r-i)),1]) & B[i,2]==(B[(i+1):(i+min(ee,r-i)),2]) & B[i,12]==(B[(i+1):(i+min(ee,r-i)),12])& zu==zu1])[1]
    }
  }
  return(h)
}



#D is your dataset after formula attribution
#create B, a dataset that represents the summed abundance of all isotopes of the respective element
B<-data.frame(C=D$C+D$C13,H=D$H,Br=D$Br79+D$Br81,Cl=D$Cl+D$Cl37,Co=D$Co,Cu=D$Cu.i+D$Cu.ii+D$Cu65.i+D$Cu65.ii,Fe=D$Fe.ii+D$Fe.iii+D$Fe54.ii+D$Fe54.iii,I=D$I,N=D$N+D$N15,Na=D$Na,Ni=D$Ni.ii+D$Ni60.ii,O=D$O+D$O18,P=D$P,S=D$S+D$S34,Zn=D$Zn.ii+D$Zn66.ii)
D<-D[order(B[,1],B[,2],B[,12]),]
B<-B[order(B[,1],B[,2],B[,12]),]
#run function above
jj<-sapply(1:(nrow(D)),FUN=get_network,r=nrow(D),B=B,indexf=(1:nrow(B)),homcheck) #homcheck represents a vector that indicates which links in the network should be accounted for
#create edgelist
edges1 = matrix(data=NA, nrow=0, ncol=2)
for(i in 1:ncol(jj)) {
  for(j in 1:7) {
    if(!is.na(jj[j,i])) { edges1 = rbind(edges1, c(i,as.numeric(jj[j,i]))) }
  }
}
#reformat edgelist
Gv = graph_from_edgelist(as.matrix(edges1), directed=FALSE)
Ak<-components(Gv)

#add network membership and size to Data
D$homnetworkmember<-Ak$membership

D$homseries<-sapply(Ak$membership,function(x) {sum(Ak$membership==x)})

  