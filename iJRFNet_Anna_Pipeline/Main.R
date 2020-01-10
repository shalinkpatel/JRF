# ---- Input   -------------------------------------------------------------------------------------------- ###
Parallel=FALSE; # -- Parallelize: TRUE for Yes; FALSE for No
N=50;  # --- construct network for N proteins with highest interquantile range
n=20   # --- sample size to consider
M=20   # -- number of permutations for iJRFNet

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

##################################################################################################################
### --- Construct Global Proteomics and RNASeq Network using iJRFNet ----------------------------------------- ###  
##################################################################################################################

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

  if (Parallel==FALSE) { # -- If not in parallel
    out<-iJRFNet(list(RNA.sel[,seq(1,n)],PRO.sel[,seq(1,n)]),model="iJRF",genes.name=gene.sel)
    out.perm<-iJRFNet_permutation(list(RNA.sel[,seq(1,n)],PRO.sel[,seq(1,n)]),model="iJRF",M=M)
  } else { # -- if parallel

    out<-array(0,c(p,p,nclasses))
    n.var=0; 
    K<-floor(length(gene.sel)/10); if (length(gene.sel)>(K*10)) K=K+1
    for (k in 1:K){ # --- Parallelize this for loop
      out.new<-iJRFNet_parallel(X=list(RNA.sel[,seq(1,n)],PRO.sel[,seq(1,n)]),genes.name=gene.sel,
                            model="iJRF",parallel=c(k,K))
      n.target<-dim(out.new$importance)[2]
      for (c in 1:nclasses) out[,seq(n.var+1,n.var+n.target),c]<-out.new$importance[,,c];
      n.var=n.var+n.target
    }
    out<-FinalScore_parallel(importance=out, model="iJRF", genes.name=gene.sel)
    
    # --- Parallelize this for loop
    out.new<-array(0,c(p,p,nclasses)); out.perm<-array(0,c((p-1)*p/2,M,nclasses))
    for (seed in 1:M){  # -- for loop over permutation
      n.var=0;
      for (k in 1:K){ 
      out.ijrf<-iJRFNet_parallel_permutation(X=list(RNA.sel[,seq(1,n)],PRO.sel[,seq(1,n)]),genes.name=gene.sel,
                                        model="iJRF",parallel=c(k,K),seed=seed)
      n.target<-dim(out.ijrf$importance)[2]
      for (c in 1:nclasses) out.new[,seq(n.var+1,n.var+n.target),c]<-out.ijrf$importance[,,c];
      n.var=n.var+n.target
    }
     out.perm[,seed,]<-FinalScore_parallel_permutation(importance=out.new, model="iJRF", genes.name=gene.sel)
    }
  } # -- End if Parallel

  # -- Derive networks and save networks in CSV file  ----------------------------------------------------------------------------- #
  Net<-Derive_network(out,out.perm,.001)
  net1<-paste(Net[[1]][,1],Net[[1]][,2],sep="-") 
  net2<-paste(Net[[2]][,1],Net[[2]][,2],sep="-") 
  vec<-unique(c(net1,net2)); 
  vec<-cbind(vec,as.data.frame(matrix(0,length(vec),2)),stringsAsFactors=FALSE); colnames(vec)<-c("Edge","RNA-Net","PRO-Net")
  vec[is.na(match(vec[,1],net1))==FALSE,2]<-1
  vec[is.na(match(vec[,1],net2))==FALSE,3]<-1
 
  write.table(vec, file=paste("iJRFNet_RNAseq_and_Global_Proteomics_Networks.csv",sep=""), 
             append = FALSE, quote = FALSE, sep = "\t",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names =FALSE, qmethod = c("escape", "double"),
             fileEncoding = "")
 
  # -- Plot Modules and save modules in CSV file  --------------------------------------------------------------------------------- #
  pdf("Hubgenes_RNAseq_Hubs.pdf")
  Hubgenes_barplot(Net,net1=1,net2=2,name.net1="RNA-Network",name.net2="PRO-Network",genes.name=gene.sel) 
  dev.off()

  pdf("Hubgenes_Global_Proteomics_Hubs.pdf")
  Hubgenes_barplot(Net,net1=2,net2=1,name.net1="RNA-Network",name.net2="PRO-Network",genes.name=gene.sel) 
  dev.off()
 
  # -- Derive cluster for each network and save in csv file ----------------------------------------------------------------------- #
  pdf("RNASeq_modules.pdf")
  cluster1<-Plot_Modules(Net[[1]],gene.sel)
  dev.off()

  pdf("Global_proteomics_modules.pdf")
  cluster2<-Plot_Modules(Net[[2]],gene.sel)
  dev.off()  

  write.table(cluster1, file=paste("iJRFNet_RNAseq_Clusters.csv",sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names =FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

  write.table(cluster2, file=paste("iJRFNet_Global_Proteomics_Clusters.csv",sep=""), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names =FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

##################################################################################################################
### --- Construct Global Proteomics and Phospho Network using iJRFNet ----------------------------------------- ###  
##################################################################################################################
