#' @title Zero truncated/inflated Nonnegative Matrix Factorization
#'
#' @description Deconvolute the gene count matrix X into a signiture matrix A and a proportion matrix W.
#'
#' @param X A gene count matrix with gene names as row names and sample(bulk) names as column names.
#' If single cell counts,S, are provided, combine single cell count matrix with mixed sample counts. X and S should have the same row names.
#' And X should be in the first few columns of the combined matrix like (X;S).
#' @param K  An integer. This is the number of cell types.
#' @param W2 An indicator matrix. If in matrix X, single cell count matrix is provided, W2 must be provided. The
#' row names of W2 are the cell types and column names are the cell names. If a cell belongs
#' to a cell type, then the column for the cell should have 1 for that cell type and 0 for the rest of the cell types.
#' If no single cell matrix provided, use NULL.
#' @param initial A list of two matrices or a matrix. The initial matrices for deconvoluting signiture matrix A and proportion matrix W.
#' If \code{updateA=1}, initial must be a signature matrix. If no initial information provided, \code{initial=NULL}.
#'
#' @param updateA An integer 0 or 1. If \code{updateA=0}, then A(signature matrix) is updated in every iteration. If \code{updateA=1},
#' then A is fixed. A signature matrix must be provided in \code{initial} option.
#'
#' @param alpha A positive vector. Must have a length of \code{K}.This is the concentration parameter of proportion matrix W.
#' If no information is provided for W,\code{alpha=NULL}.
#'
#' @rdname fixed_point_algorithm
#'
#' @return \code{fixedpoint} returns a list of deconvultion results containing a signiture matrix A,
#' a proportion matrix W and the value of the likelihood for the last iteration.
#'
#' @examples
#' #generate a count matrix X(bulk counts) without single cell counts.
#' list1<-samplecounts(N=10,M1=5,Ks=0,Kb=2,Nm=2,Na=0,Nh=5,alpha=c(1,1),missing="logit")
#' X<-list1$X
#' fit1<-fixedpoint(X,2)
#' fit1$A
#' fit1$W
#' #use the output from fit1 as initial matrices for fit1_Bayes
#' fit1_Bayes<-Bayes_fixedpoint(X, 2,initial=list(fit1$A,fit1$W))
#'
#' #generate a count matrix X(bulk counts) with single cell counts.
#' list2<-samplecounts(N=15,M1=5,Ks=3,P=c(10,15,20),Kb=3,Nm=2,Na=0,Nh=5,alpha=c(1,2,3),missing="logit")
#' X<-list2$X
#' #Check here: the first 5 columns of W are for bulks and the rest of the columns are for single cells.
#' trueW<-list2$W
#' X<-list2$X
#' fit2<-fixedpoint_alpha(X, K=3,W2=trueW[,-c(1:5)], updateA=0,  initial=NULL, alpha=c(1,2,3))
#' #use the output from fit2 to estimate signature matrix A and not change A in fit2_Bayes
#' fit2_Bayes<-Bayes_fixedpoint_alpha(X, 3,updateA=1,initial=fit2$A, alpha=c(1,2,3))
#'
#' @export
#' @useDynLib iNMF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
fixedpoint <- function(X, K,W2=NULL, initial=NULL){

  N <- dim(X)[1]
  M <- dim(X)[2]

  if(!is.null(W2)){
    M1 <- M - dim(W2)[2]
    if(!is.list(initial)){

      Atemp <- X[,(M1+1):M]%*%diag(1/apply(X[,(M1+1):M], 2, sum))
      A  <-   (Atemp%*%t(W2))%*%diag(1/apply(t(W2), 2, sum))
      index <- which(apply(t(W2), 2 , sum) == 0)
      if(length(index)> 0){
        A[, index] <- t(rdirichlet(length(index), rep(1,N)))
      }


      if(!is.null(markergene)){
        for(ii in 1:length(markergene)){
          A[markergene[ii], -celltype[ii]] <-  rep(1e-8, K-1)
        }
      }

      if(is.null(markergene)){
        A  <-  A + 0.02*t(rdirichlet(K, rep(1,N)))
      }
      A <-  A %*% diag(1/apply(A, 2, sum))
      W  <- cbind(t(rdirichlet(M1, rep(1, K))) , W2)
    }else{
      A <-  initial$A
      W <- initial$W
    }
  }else{
    M1 <- M
    if(!is.list(initial)){
      A  <-  t(rdirichlet(K, rep(1,N)))
      W  <- cbind(t(rdirichlet(M1, rep(1, K))) , W2)
    }else{
      A <-  initial$A
      W <- initial$W
    }
  }

  Xv <- matrix(X, nc=1)
  Av <- matrix(A, nc=1)
  Wv <- matrix(t(W), nc=1)
  lambda <- apply(X, 2, sum)
  output <- fitfixedpoint(Xv, lambda, N, M, M1, Av, Wv, K)
  Amat <- matrix(output[1:(N*K)], nc=K)
  Wmat<-t(matrix(output[(N*K+1):(N*K +M*K)], nc=K))
  likevec = output[(N+M)*K+1]
  return(list(A=Amat, W=Wmat, likelihood=likevec))
}

#' @rdname fixed_point_algorithm
#' @return \code{fixedpoint_alpha} returns a list of deconvultion results containing a signiture matrix A,
#' a proportion matrix W and the value of the likelihood for the last iteration.
#' @export
fixedpoint_alpha <- function(X, K,W2=NULL, updateA=0,  initial=NULL, alpha){

  Words<-X
  N <- dim(Words)[1]
  M <- dim(Words)[2]

  if(!is.null(W2)){
    M1 <- M - dim(W2)[2]
    if(!is.list(initial)){
      if(is.matrix(initial)){
        A<-initial
      }else{
        Atemp <- Words[,(M1+1):M]%*%diag(1/apply(Words[,(M1+1):M], 2, sum))
        A  <-   (Atemp%*%t(W2))%*%diag(1/apply(t(W2), 2, sum))
        index <- which(apply(t(W2), 2 , sum) == 0)
        if(length(index)> 0){
          A[, index] <- t(rdirichlet(length(index), rep(1,N)))
        }
        A  <-  A + 0.02*t(rdirichlet(K, rep(1,N)))
        A <-  A %*% diag(1/apply(A, 2, sum))
      }
      W  <- cbind(t(rdirichlet(M1, alpha)) , W2)
    }else{
      A <-  initial$A
      W <- initial$W
    }
  }else{
    M1 <- M
    if(!is.list(initial)){
      if(is.matrix(initial)){
        A<-initial
      }else{
        A  <-  t(rdirichlet(K, rep(1,N)))
      }
      W  <- cbind(t(rdirichlet(M1, alpha)) , W2)
    }else{
      A <-  initial$A
      W <- initial$W
    }
  }

  Xv <- matrix(Words, nc=1)
  Av <- matrix(A, nc=1)
  Wv <- matrix(t(W), nc=1)
  lambda <- apply(Words, 2, sum)

  output <- fitfixedpointAlpha(Xv, lambda, N, M, M1, Av, updateA, Wv, K, alpha)
  Amat <- matrix(output[1:(N*K)], nc=K)
  Wmat<-t(matrix(output[(N*K+1):(N*K +M*K)], nc=K))
  likevec = output[(N+M)*K+1]
  return(list(A=Amat, W=Wmat, likelihood=likevec))
}

#' @rdname fixed_point_algorithm
#' @return \code{Bayes_fixedpoint} returns a list of deconvultion results containing a signiture matrix A,
#' a proportion matrix W and the value of the likelihood for the last iteration.
#' @export
Bayes_fixedpoint <- function(X, K,W2=NULL, initial=NULL){
  Words<-X
  N <- dim(Words)[1]
  M <- dim(Words)[2]
  if(!is.null(W2)){
    M1 <- M - dim(W2)[2]
    if(!is.list(initial)){
      A  <-   Words[, (M1+1):M]%*%t(W2) + 0.02*t(rdirichlet(K, rep(1,N)))
      A <-  A %*% diag(1/apply(A, 2, sum))
      W  <- cbind(t(rdirichlet(M1, rep(1, K))) , W2)
    }else{
      A <-  initial$A
      W <- initial$W
    }
  }else{
    M1 <- M
    if(!is.list(initial)){
      A  <-  t(rdirichlet(K, rep(1,N)))
      W  <- cbind(t(rdirichlet(M1, rep(1, K))) , W2)
    }else{
      A <-  initial$A
      W <- initial$W
    }
  }

  Total <- apply(Words, 2, sum)
  MCcount <-3000
  burnin <- 1001
  L <- MCcount - burnin +1

  OmegaM <- rep(0, length=N*M*L)
  SM <- rep(0, length=N*M*L)
  TauM <- rep(0, length=M*L)
  KappaM <- rep(0, length=M*L)

  EqS <- rep(0, N*M)
  EqSKappa <- rep(0, N*M)
  EqSTau <- rep(0, N*M)
  EqKappa <- rep(0, M)
  EqKappasq<- rep(0, M)
  EqTau<- rep(0, M)
  EqTausq <- rep(0, M)
  EqOmegaKappasq <- rep(0, N*M)
  EqOmegaKappaTau <- rep(0, N*M)
  EqOmegaTausq <- rep(0, N*M)

  output <-MCMCsample2(matrix(Words, nc=1), Total, N, M, M1, K, matrix(A, nc=1), matrix(t(W), 1),  MCcount,burnin, EqS, EqSKappa, EqSTau, EqKappa, EqKappasq, EqTau, EqTausq, EqOmegaKappasq, EqOmegaKappaTau, EqOmegaTausq)
  return(list(A= matrix(output[1:(N*K)], nc=K), W=t(matrix(output[(N*K+1):((N+M)*K)], nc=K)), likelihood=output[(N+M)*K+1]))
}

#' @rdname fixed_point_algorithm
#' @return \code{Bayes_fixedpoint_alpha} returns a list of deconvultion results containing a signiture matrix A,
#' a proportion matrix W and the value of the likelihood for the last iteration.
#' @export
Bayes_fixedpoint_alpha <- function(X, K,W2=NULL,updateA = 0, initial=NULL, alpha){

  Words<-X
  N <- dim(Words)[1]
  M <- dim(Words)[2]

  if(!is.null(W2)){
    M1 <- M - dim(W2)[2]
    if(!is.list(initial)){
      if(is.matrix(initial)){
        A<-initial
      }else{
      A  <-   Words[, (M1+1):M]%*%t(W2) + 0.02*t(rdirichlet(K, rep(1,N)))
      A <-  A %*% diag(1/apply(A, 2, sum))
      }
      W  <- cbind(t(rdirichlet(M1, rep(1, K))) , W2)
    }else{
      A <-  initial$A
      W <- initial$W
    }
  }else{
    M1 <- M
    if(!is.list(initial)){
      if(is.matrix(initial)){
        A<-initial
      }else{
      A  <-  t(rdirichlet(K, rep(1,N)))
      }
      W  <- cbind(t(rdirichlet(M1, rep(1, K))) , W2)
    }else{
      A <-  initial$A
      W <- initial$W
    }
  }

  Total <- apply(Words, 2, sum)
  MCcount <-3000
  burnin <- 1001
  L <- MCcount - burnin +1

  OmegaM <- rep(0, length=N*M*L)
  SM <- rep(0, length=N*M*L)
  TauM <- rep(0, length=M*L)
  KappaM <- rep(0, length=M*L)

  EqS <- rep(0, N*M)
  EqSKappa <- rep(0, N*M)
  EqSTau <- rep(0, N*M)
  EqKappa <- rep(0, M)
  EqKappasq<- rep(0, M)
  EqTau<- rep(0, M)
  EqTausq <- rep(0, M)
  EqOmegaKappasq <- rep(0, N*M)
  EqOmegaKappaTau <- rep(0, N*M)
  EqOmegaTausq <- rep(0, N*M)

  output <- MCMCsample2alpha(matrix(Words, nc=1),Total, N,M, M1, K,matrix(A, nc=1), updateA, matrix(t(W), 1), alpha, MCcount,burnin,EqS, EqSKappa, EqSTau, EqKappa, EqKappasq, EqTau, EqTausq, EqOmegaKappasq, EqOmegaKappaTau,EqOmegaTausq)
  return(list( A= matrix(output[1:(N*K)], nc=K), W=t(matrix(output[(N*K+1):((N+M)*K)], nc=K)), likelihood=output[(N+M)*K+1]))
}

