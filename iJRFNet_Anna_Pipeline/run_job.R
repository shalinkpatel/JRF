c=as.numeric(Sys.getenv("LSB_JOBINDEX"))   
DIR<-" ..." 
setwd(DIR)

source('Main_parallel.R')
test.I(parallel=1,n=50,perm=FALSE)