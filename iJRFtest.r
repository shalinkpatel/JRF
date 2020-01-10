library(iJRFNet)

nclasses=2
n1<-n2<-20
p<-400
genes.name<-paste("G",seq(1,p),sep="")

data1<-matrix(rnorm(p*n1),p,n1)
data2<-matrix(rnorm(p*n2),p,n1)

out.iJRFNet<-iJRFNet(X=list(data1,data2),genes.name=genes.name, model = "iJRF")

for (j in 1:8){
  final.net<-Plot_Modules(out.iJRFNet[sample(dim(out.iJRFNet)[1],j*50),c(1,2)], genes.name=genes.name)
}