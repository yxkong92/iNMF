% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixed_point_algorithm.R
\name{fixedpoint}
\alias{fixedpoint}
\alias{fixedpoint_alpha}
\alias{Bayes_fixedpoint}
\alias{Bayes_fixedpoint_alpha}
\title{Zero truncated/inflated Nonnegative Matrix Factorization}
\usage{
fixedpoint(X, K, W2 = NULL, initial = NULL)

fixedpoint_alpha(X, K, W2 = NULL, updateA = 0, initial = NULL, alpha)

Bayes_fixedpoint(X, K, W2 = NULL, initial = NULL)

Bayes_fixedpoint_alpha(X, K, W2 = NULL, updateA = 0, initial = NULL, alpha)
}
\arguments{
\item{X}{A gene count matrix with gene names as row names and sample(bulk) names as column names.
If single cell counts,S, are provided, combine single cell count matrix with mixed sample counts. X and S should have the same row names.
And X should be in the first few columns of the combined matrix like (X;S).}

\item{K}{An integer. This is the number of cell types.}

\item{W2}{An indicator matrix. If in matrix X, single cell count matrix is provided, W2 must be provided. The
row names of W2 are the cell types and column names are the cell names. If a cell belongs
to a cell type, then the column for the cell should have 1 for that cell type and 0 for the rest of the cell types.
If no single cell matrix provided, use NULL.}

\item{initial}{A list of two matrices or a matrix. The initial matrices for deconvoluting signiture matrix A and proportion matrix W.
If \code{updateA=1}, initial must be a signature matrix. If no initial information provided, \code{initial=NULL}.}

\item{updateA}{An integer 0 or 1. If \code{updateA=0}, then A(signature matrix) is updated in every iteration. If \code{updateA=1},
then A is fixed. A signature matrix must be provided in \code{initial} option.}

\item{alpha}{A positive vector. Must have a length of \code{K}.This is the concentration parameter of proportion matrix W.
If no information is provided for W,\code{alpha=NULL}.}
}
\value{
\code{fixedpoint} returns a list of deconvultion results containing a signiture matrix A,
a proportion matrix W and the value of the likelihood for the last iteration.

\code{fixedpoint_alpha} returns a list of deconvultion results containing a signiture matrix A,
a proportion matrix W and the value of the likelihood for the last iteration.

\code{Bayes_fixedpoint} returns a list of deconvultion results containing a signiture matrix A,
a proportion matrix W and the value of the likelihood for the last iteration.

\code{Bayes_fixedpoint_alpha} returns a list of deconvultion results containing a signiture matrix A,
a proportion matrix W and the value of the likelihood for the last iteration.
}
\description{
Deconvolute the gene count matrix X into a signiture matrix A and a proportion matrix W.
}
\examples{
#generate a count matrix X(bulk counts) without single cell counts.
list1<-samplecounts(N=10,M1=5,Ks=0,Kb=2,Nm=2,Na=0,Nh=5,alpha=c(1,1),missing="logit")
X<-list1$X
fit1<-fixedpoint(X,2)
fit1$A
fit1$W
#use the output from fit1 as initial matrices for fit1_Bayes
fit1_Bayes<-Bayes_fixedpoint(X, 2,initial=list(A=fit1$A,W=fit1$W))

#generate a count matrix X(bulk counts) with single cell counts.
list2<-samplecounts(N=15,M1=5,Ks=3,P=c(10,15,20),Kb=3,Nm=2,Na=0,Nh=5,alpha=c(1,2,3),missing="logit")
X<-list2$X
#Check here: the first 5 columns of W are for bulks and the rest of the columns are for single cells.
trueW<-list2$W
fit2<-fixedpoint_alpha(X, K=3,W2=trueW[,-c(1:5)], updateA=0,  initial=NULL, alpha=c(1,2,3))
#use the output from fit2 to estimate signature matrix A and not change A in fit2_Bayes
fit2_Bayes<-Bayes_fixedpoint_alpha(X, 3,updateA=1,initial=fit2$A, alpha=c(1,2,3))

}
