#' @title Sample mixed or pure count matrix X
#'
#' @description Sample a gene count matrix X(pure,mixed or both) from a signiture matrix A and a proportion matrix W.
#' samplecounts<-function(N,M1,Ks,Kb,Nm,Na,Nh,alpha,missing="logit",mkappa=1,mtau=1.5)
#' @param N An integer. The number of genes. The sum of \code{Nm},\code{Na} and \code{Nh} should be less than \code{N}.
#' @param M1  An integer. The number of bulks.
#' @param Ks An integer. The number of cell types for pure(single cell) samples.
#' @param P  An integer vector. The number of cells from each cell types. The sum of \code{P} is the total number of single cells. And
#' the length of \code{P} must be \code{Ks}. \code{P} is NULL when  \code{Ks}=0.
#' @param Kb An integer. The number of cell types for mixed(bulk) samples.
#' @param Nm A positive integer. The number of marker genes in one cell type.
#' @param Na An integer. The number of anti-marker genes in one cell type.
#' @param Nh A positive integer. The number of house-keeping genes in one cell type.
#' @param alpha A positive vector. Must have a length of \code{Kb}.This is the concentration parameter of proportion matrix W.
#' @param missing A value from 0 to 1 or a character. If missing count rate is fixed, then input
#' an value between 0 and 1(included). Another option is \code{"logit"}. This missing rate is based on A and W.
#' @param mkappa A scalar. Control the missing rate when \code{missing="logit"}.
#' @param mtau A scalar. Control the missing rate when \code{missing="logit"}. The defualt setting of
#' \code{mkappa} and \code{mtau} will generate around 0.1 missing in X matrix.
#'
#' @rdname sample_count_matrix
#'
#' @return \code{sample_count_matrix} returns a list of matrices. \code{X} the count matrix.
#' \code{A} the signature matrix. \code{W} the proportion matrix.
#'
#'
#' @examples
#' #generate a count matrix X(bulk counts) without single cell counts.
#' list1<-samplecounts(N=10,M1=5,Ks=0,Kb=2,Nm=2,Na=0,Nh=5,alpha=c(1,1),missing="logit")
#' X<-list1$X
#' trueA<-list1$A
#' trueW<-list1$W
#' #generate a count matrix X(bulk counts) with single cell counts.
#' list2<-samplecounts(N=15,M1=5,Ks=3,P=c(10,15,20),Kb=3,Nm=2,Na=0,Nh=5,alpha=c(1,2,3),missing="logit")
#' X<-list2$X
#' #Check here: the first 5 columns of W are for bulks and the rest of the columns are for single cells.
#' trueW<-list2$W
#'
#' @export
samplecounts<-function(N,M1,Ks,P,Kb,Nm,Na,Nh,alpha,missing="logit",mkappa=1,mtau=1.5){
  library(MCMCpack)
  library(stringr)
  K <- max(Ks,Kb)

  A <- matrix(0, nc=K, nr=N)
  if (Na==0){
    for(k in 1: K){
      for(i in 1:Nm){
        A[i+(k-1)*Nm,k]=exp(rnorm(1))
      }
    }

    for(i in 1:Nh){
      val = exp(rnorm(1))
      for(k in 1:K){
        A[Nm*K+i, k] = val
      }
    }

    for(i in ( Nm*K + Nh +1):N){
      for(j in 1:K){
        A[i,j] <- exp(rnorm(1))
      }

    }

    tot1 <- apply(A[1:((Nm)*K),], 2, sum)
    A1<- A[1:((Nm)*K),]%*%diag(1/tot1*(1-Nh/N))

    tot2 <- apply(A[((Nm)*K+1):(N), ], 2, sum)
    A2<- A[((Nm)*K+1):(N), ] %*%diag(1/tot2*Nh/N)

    A <- rbind(A1, A2)

  }else{
    for(k in 1: K){
      for(i in 1:Nm){
        A[i+(k-1)*Nm,k]=exp(rnorm(1))
      }
    }

    for(k in 1:K){
      for(i in 1:Na){
        A[i+(Na*(k-1) +Nm*K), -k] = exp(rnorm((K-1)))
      }
    }

    for(i in 1:Nh){
      val = exp(rnorm(1))
      for(k in 1:K){
        A[(Nm+Na)*K+i, k] = val
      }
    }

    for(i in ( (Na+Nm)*K + Nh +1):N){
      for(j in 1:K){
        A[i,j] <- exp(rnorm(1))
      }

    }

    tot1 <- apply(A[1:((Nm+Na)*K),], 2, sum)
    A1<- A[1:((Nm+Na)*K),]%*%diag(1/tot1*(1-Nh/N))

    tot2 <- apply(A[((Nm+Na)*K+1):(N), ], 2, sum)
    A2<- A[((Nm+Na)*K+1):(N), ] %*%diag(1/tot2*Nh/N)

    A <- rbind(A1, A2)

  }

  Depth1 <- rpois(M1, 50*N)

  W1 <- t(rdirichlet(M1, alpha))
  a<-numeric(K)
  a[1]<-1
  if (Ks==0){
    M <- M1
    cn<-str_c("bulk", 1:M)
    W<-W1
    Depth <- Depth1
  }else{
    L<-sum(P)
    M <- M1+L
    cn1<-str_c("bulk", 1:M1)
    cn2<-str_c("Cell", 1:L)
    cn<-c(cn1,cn2)
    Depth2 <- rnbinom(L, size=2,  mu= 2*N)
    Depth <- c(Depth1, Depth2)
    w <- matrix(0,nr=K, nc=L)
    idx<-0
    for (i in 1:Ks){
      w[i,c((1+idx):(P[i]+idx))]<-1
      idx<-idx+P[i]
    }
    W2<-w
    if (Kb<Ks) {
      W1 <- rbind(W1,matrix(0,nrow = (Ks-Kb),ncol = M1))
    }else{
      W1<-W1
      W2<-W2
    }
    W <- cbind(W1, W2)
  }
  trueW  <- W
  trueA <-  A
  WordsbyClass <- array(NA, dim=c(N, M, K))
  for (j in 1:M){
    Z <- rmultinom(Depth[j],1,  W[,j])  ## K by Depth matrix
    Classcounts <- apply(Z, 1, sum)
    for(l in 1:K){
      Dpart = rmultinom(1,Classcounts[l], A[,l])
      WordsbyClass[,j,l] <- Dpart
    }
  }
  Words <- apply(WordsbyClass, c(1,2), sum)
  if (is.character(missing)){
    kappavec <- rnorm(M, mkappa, 0.5)
    tauvec <- rnorm(M, mtau*N, 0.15*N)
    probmat <- A%*%W
    for(i in 1:N){
      for(j in 1:M){
        successprob <- 1/(1+exp(-kappavec[j] - tauvec[j]*probmat[i,j]))
        z<-rbinom(1, 1, successprob)
        Words[i,j] <- z*Words[i,j] + (1-z)*0
      }
    }
  }else{
    for(i in 1:N){
      for(j in 1:M){
        successprob <- 1-missing
        z<-rbinom(1, 1, successprob)
        Words[i,j] <- z*Words[i,j] + (1-z)*0
      }
    }
  }
  X <- Words
  rownames(X)<-str_c("gene", 1:N)
  colnames(X)<-cn
  rownames(A)<-str_c("gene", 1:N)
  colnames(A)<-str_c("CT", 1:K)
  rownames(W)<-str_c("CT", 1:K)
  colnames(W)<-cn
  return(list(A=A,W=W,X=X))
}

