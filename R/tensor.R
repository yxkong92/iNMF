#' @title Tensor Decomposition for LDA Model
#'
#' @description Find the concentration parameter \code{alpha} of proportion matrix W.
#'
#' @param X A gene count matrix with gene names as row names and sample(bulk) names as column names. Single cell
#' data are not provided in this case.
#'
#' @param a0  A positive value. This is the sum of concentration parameter \code{alpha}.
#' @param P An integer. The number of eigenvalues of M2 matrix to use in finding \code{alpha}. When
#' \code{P} is larger than the number of positive eigenvalues of M2, the algorithm will automatically choose the
#' largest number of positive eigenvalues of M2 as \code{P}.
#' @param K An integer. The length of \code{alpha}.
#'
#' @param L A positive integer. The number of starting points to use in the algorithm.
#'
#' @param N A positive integer. The number of iterations to use in the algorithm.
#'
#' @rdname tensor
#'
#' @return \code{tensor_eigen} returns a vector of the rates of the concentration parameter \code{alpha}.
#'
#' @references Animashree Anandkumar, Rong Ge, Daniel Hsu, Sham M. Kakade, Matus Telgarsky,
#' \emph{Tensor Decompositions for Learning Latent Variable Models}. Journal of Machine Learning Research 15 (2014) ,
#' pages 2773-2832.
#'
#' @examples
#' #generate a count matrix X(bulk counts) without single cell counts.
#' list1<-samplecounts(N=200,M1=150,Ks=0,Kb=3,Nm=5,Na=0,Nh=20,alpha=c(100,200,300),missing=0)
#' X<-list1$X
#' a<-tensor_eigen(X,a0=600,p=10,K=3,L=100,N=150)
#' a
#'
#' @useDynLib iNMF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
tensor_eigen<-function(X,a0,p,K,L,N){
  m1<-rowMeans(t(t(X)/colSums(X)))
  n<-nrow(X)
  m<-ncol(X)
  P<-p
  Xv <- matrix(X, nc=1)
  C2<-c2cpp(Xv,n,m)
  C2mtx<-matrix(C2,nrow = n)
  m2<-C2mtx-a0*(m1%o%m1)/(a0+1)
  sv <- eigen(m2)
  if (sum(sv$val>0)<P){
    P<-sum(sv$val>0)
  }
  U <-sv$vec[,1:P]
  invD <- diag(1/sqrt(sv$val[1:P]))
  W<-U%*%(invD)
  Wv <- matrix(W, nc=1)
  m3<-m3tildcpp(Xv,Wv,C2,m1,a0,n,m,P)
  evalue<-numeric()
  for (i in 1:K){
    output <- Tspowercpp(m3, P, L, N)
    evalue[i]<-output$lambda.hat
    m3<-output$dtensor
  }
  a<-1/evalue^2
  return(a*a0)
}

Tspowercpp<-function(m3,P,L,N){
  thetaN<-matrix(0,nrow = P,ncol = L)
  rqt<-numeric()
  for (tao in 1:L){
    theta<-Ksphere(P)
    theta<-Npowercpp(N,m3,theta,P)
    thetaN[,tao]<-theta
    rqt[tao]<-lambdacpp(P,m3,theta)
  }
  taost<-which.max(rqt)
  thetaNst<-thetaN[,taost]
  theta.hat<-as.numeric(Npowercpp(N,m3,thetaNst,P))
  lambda.hat<-lambdacpp(P,m3,theta.hat)
  otpd<-as.numeric(lambda.hat*(theta.hat%o%theta.hat%o%theta.hat))
  dtensor<-m3-otpd   ##vector form better change it
  return(list(theta.hat=theta.hat,lambda.hat=lambda.hat,dtensor=dtensor))
}

