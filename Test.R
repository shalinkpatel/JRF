# --- Generate data sets
source("ptmJRF.R")
nclasses=2
n1<-n2<-20
p<-20
genes.name<-paste("G",seq(1,p),sep="")   # genes name
ptm.name<-c("G1","G2","G3","G3","G4","G5","G1")   # ptm name
p.ptm<-length(ptm.name)
data1<-matrix(rnorm(p.ptm*n2),p.ptm,n1)       # generate PTM data
data2<-matrix(rnorm(p*n1),p,n1)       # generate global proteomics data
# --- Run JRF and obtain importance score of interactions for each class
out<-ptmJRF(X=list(data1,data2),genes.name=genes.name,ptm.name=ptm.name)


Plot_Modules(out, genes.name = genes.name) 
