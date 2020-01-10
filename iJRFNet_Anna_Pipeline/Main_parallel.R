run<-function(parallel,n,perm=FALSE,perm.seed=NULL){ 

  # -- parallel: number parallel job 
  # -- n: sample size 
  # -- perm: TRUE for permutation, FALSE for standard iJRF
  # -- perm.seed: permutation seed
    
### --- load library -------------------------------------------------------------------------------------- ###
 DIR<-"~/Dropbox/iJRFNet CPTAC3 Pipeline"
 setwd(DIR)
 install.packages("iJRFNet_1.1-4.tar.gz")
 library(iJRFNet)
 source("Remove_replicates.R")
 
### --- Load Data ----------------------------------------------------------------------------------------- ###
  # --- Load RNAseq Data (p x n matrix)
  DIR<-"~/Dropbox/CPTAC3/Pipeline/Hui Pipeline/data/Breast"
  setwd(DIR)
  RNA<-read.table("RNAseq.txt")
  gene.ID<-rownames(RNA); 
  gene.ID<-unlist(strsplit(gene.ID, split='|', fixed=TRUE))[seq(1,length(gene.ID)*2,2)]
  gene.ID<-gene.ID[rowSums(is.na(RNA))==0]; RNA<-RNA[rowSums(is.na(RNA))==0,]; 
  RNA<-replicates(RNA,gene.ID); 
  gene.ID<-rownames(RNA)

 # --- Load Proteomics (output from imputation step)
  map<-read.table("Annotation_Protein_Gene.txt",sep="\t",header=TRUE) # -- mapping genes to proteins
  Global.pro<-read.table('Proteome_Protein_remove_batch_impute.txt',header=TRUE,sep="\t") 
  Pro.ID<-rownames(Global.pro);
  m.g<-match(Pro.ID,map[,3]); Pro.ID<-map[m.g[is.na(m.g)==FALSE],4];  Global.pro<-Global.pro[is.na(m.g)==FALSE,]
  Global.pro<-replicates(Global.pro,as.character(Pro.ID));  Pro.ID<-rownames(Global.pro)

# -- Match protein ID and samples for global and RNAseq data
  m.g<-match(gene.ID,Pro.ID) 
  Global.pro<-Global.pro[m.g[is.na(m.g)==FALSE],]
  RNA<-RNA[is.na(m.g)==FALSE,]; gene.ID<-gene.ID[is.na(m.g)==FALSE]

  sample.pro<-colnames(Global.pro) 
  sample.rna<-apply(sapply(strsplit(colnames(RNA), '.',fixed=TRUE), '[', c(2,3)),2,paste,collapse=".")
  m.g<-match(sample.pro,sample.rna);
  RNA<-RNA[,m.g[is.na(m.g)==FALSE]]
  Global.pro<-Global.pro[,is.na(m.g)==FALSE]
 
  iqr.pro<-apply(Global.pro,1,function(x) {quantile(x,.75)-quantile(x,.25)})
  iqr.rna<-apply(RNA,1,function(x) {quantile(x,.75)-quantile(x,.25)})
  th.rna<-sort(iqr.rna,decreasing=TRUE)[N]; th.pro<-sort(iqr.pro,decreasing=TRUE)[N]
 
  gene.sel<-unique(c(gene.ID[iqr.rna>th.rna],gene.ID[iqr.pro>th.pro]))
  p<-length(gene.sel); nclasses=2
  PRO.sel<-Global.pro[match(gene.sel,gene.ID),] 
  RNA.sel<-RNA[match(gene.sel,gene.ID),] 

# -- Run iJRFNet

    K<-floor(length(gene.sel)/10); if (length(gene.sel)>(K*10)) K=K+1
    
    if (perm == FALSE){ # run standard iJRFNet
      out.new<-iJRFNet_parallel(X=list(RNA.sel[,seq(1,n)],PRO.sel[,seq(1,n)]),genes.name=gene.sel,
                            model="iJRF",parallel=c(parallel,K))
      save(out.new$importance,file=paste("iJRF_parallel_",parallel,".rda",sep="")) # -- save on file importance scores
      
    } else { #-- run permutation iJRFNet
      out.ijrf<-iJRFNet_parallel_permutation(X=list(RNA.sel[,seq(1,n)],PRO.sel[,seq(1,n)]),genes.name=gene.sel,
                                             model="iJRF",parallel=c(parallel,K),seed=perm.seed)
      save(out.ijrf$importance,file=paste("iJRF_permutation_",perm.seed,"_parallel_",parallel,".rda",sep="")) # -- save on file
      
    }
}
 