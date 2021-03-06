\name{iJRFNet}
\alias{iJRFNet}
\title{Derive importance scores for function of class iJRFNet.}
\usage{
iJRFNet(X, W=NULL, ntree=NULL, mtry=NULL, model=NULL, genes.name,
        ptm.name=NULL)}
\arguments{
\item{X}{List object containing expression data for each class, \code{X=list(x_1,x_2, ... )} where \code{x_j} is a \code{(p x n_j)} matrix with rows corresponding to genes and columns to samples. Rows need to be the same across objects, while samples can vary. Missing values are not allowed. If \code{model="ptmJRF"}, the first object of the list must contain the expression of post translational modification variables. Only in this case, the number of variables in the first object might differ from that of other objects. Rows of \code{X[[1]]} does not need to be ordered in a specific way.
}

\item{W}{\code{(p x p)} Optional symmetric matrix containing sampling scores. When omitted, the standard JRF algorithm without weighted sampling scheme will be implemented. Element \code{(i,j)} contains interaction score \code{(i - j)}. Scores must be non-negative. Larger value of sampling score corresponds to higher likelihood of gene \code{i} interacting with gene \code{j}. Columns and rows of \code{W} must be in the same order as the columns of \code{X}. Sampling scores \code{W} are computed considering one prior data such as protein-protein interactions.}

\item{ntree}{Numeric value: number of trees. If omitted, 1000 trees are considered.}

\item{mtry}{Numeric value: number of predictors to be sampled at each node. If omitted, \code{mtry} is set to the square root of the number of variables.}


\item{model}{ Variable indicating which iJRFNet model needs to be run. Takes values in \code{c("iJRF", "iRafNet","ptmJRF")} }

\item{genes.name}{Vector containing genes name. The order needs to match the rows of \code{x_j}.}

\item{ptm.name}{List of post translational modification variables in protein domain. This list must be ordered as rows of \code{X[[1]]}. }

}
\value{
A three dimensional matrix (\code{I} x \code{M} x \code{C}) with \code{I} being the number of total interactions, \code{M} the number of permutations and \code{C} the number of classes. Element \code{(i,j,k)} corresponds to the importance score for interaction \code{i}, permuted data \code{j} and class \code{k}. 
}
\description{
This function computes importance score for \code{M} permuted data sets. Sample labels of target genes are randomly permuted and JRF is implemented. Resulting importance scores can be used to derive an estimate of FDR.  
}
\examples{

 # --- Generate data sets
 nclasses=2               # number of data sets / classes
 n1<-n2<-20               # sample size for each data sets
 p<-5                   # number of variables (genes)
 genes.name<-paste("G",seq(1,p),sep="")   # genes name
 
 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n1)       # generate data2
 
 ##---------------------------------------------------------------------------##
 ## ---  Run iJRFNet 
 
 # --- Obtain importance scores of gene-gene (protein-protein) interactions
  out<-iJRFNet(X=list(data1,data2),genes.name=genes.name,model="iJRF")

 ##---------------------------------------------------------------------------##
 ## ---  Run iRafNet 
 
 W<-abs(matrix(rnorm(p*p),p,p))    # generate weights for interactions

 # --- Obtain importance scores of gene-gene (protein-protein) interactions
 out<-iJRFNet(X=list(data1),W=W,genes.name=genes.name,model="iRafNet")


}
\references{
Petralia, F., Song, W.M., Tu, Z. and Wang, P. (2016). New method for joint network analysis reveals common and different coexpression patterns among genes and proteins in breast cancer. \emph{Journal of proteome research}, \bold{15}(3), pp.743-754.

A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. \emph{R News} \bold{2}, 18--22.
}
