
replicates <- function(data,geneID){ 
  
  # -- When gene/protein replicates exist choose the one with highest interquantile range
  # -- data : (p x n) matrix with rowmanes indicating geneID/proteinID
  
  data<-data[geneID!="?",]; geneID<-geneID[geneID!="?"];
  dup<-duplicated(geneID)  
  geneDup<-unique(geneID[dup])
  dataNew<-data[is.na(match(geneID,geneDup)),]
  geneNew<-geneID[is.na(match(geneID,geneDup))]
  rownames(dataNew)<-geneNew
  
  
  for (j in 1:length(geneDup)){
    data.dup<-data[geneID==geneDup[j],]
    iqr<-apply(data.dup,1,function(x) {quantile(x,.75)-quantile(x,.25)})
    dataNew<-rbind(dataNew,data.dup[iqr==max(iqr),])
  }
  rownames(dataNew)<-c(geneNew,geneDup)
  
  return(dataNew)
}