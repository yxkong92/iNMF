#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector Ksphere(const int K) {
  NumericVector theta;
  theta = rnorm(K);
  theta = theta/sqrt(sum(pow(theta,2)));
  return theta;
}
// [[Rcpp::export]]
NumericVector c2cpp(NumericVector X,int n,int m){ //m ncol.X, n nrow.X
  double s;
  NumericVector C2(n*n),cct(n*n),l(m);
  for(int i = 0; i < m; i++){
    s=0;
    for(int j = 0; j < n; j++){
      s=s+X[i*n+j];
      for(int k = 0; k < n; k++){
        if (j==k){
          cct[j*n+k]=X[i*n+j]*X[i*n+k]-X[i*n+k];
        }else{
          cct[j*n+k]=X[i*n+j]*X[i*n+k];
        }
      }
    }
    l[i]=s;
    C2=C2+cct/(l[i]*(l[i]-1.0));
  }
  return C2/m;
}
// [[Rcpp::export]]
NumericVector m3tildcpp(NumericVector X,NumericVector W,NumericVector C2,NumericVector m1,double a0,int n,int m,int p){
  NumericVector m3(p*p*p),sm3(p*p*p);                                   //p ncol.W, n nrow.X , m ncol.X
  NumericVector l(m),xw(p),mw(p);
  double s,s1,parta1,parta2,parta3,parta4,parta5,parta,partb,partc,partd,parte;
  for(int i = 0; i < m; i++){
    s=0;
    for(int j = 0; j < n; j++){
      s=s+X[i*n+j];
    }
    l[i]=s;
  }
  for(int i=0;i<m;i=i+1){
    for(int j=0;j<p;j=j+1){
      s1=0.0;
      for (int k=0;k<n;k=k+1){
        s1=s1+X[i*n+k]*W[j*n+k];
      }
      xw[j]=s1;
    }
    for(int i1=0;i1<p;i1=i1+1){
      for(int i2=0;i2<p;i2=i2+1){
        for(int i3=0;i3<p;i3=i3+1){
          parta1= xw[i1]*xw[i2]*xw[i3];
          parta2=0.0;
          parta3=0.0;
          parta4=0.0;
          parta5=0.0;
          for (int j=0;j<n;j=j+1){
            parta2=parta2+X[i*n+j]*(W[i1*n+j]*W[i2*n+j]*W[i3*n+j]);
            parta3=parta3+X[i*n+j]*(W[i1*n+j]*W[i2*n+j]);
            parta4=parta4+X[i*n+j]*(W[i1*n+j]*W[i3*n+j]);
            parta5=parta5+X[i*n+j]*(W[i2*n+j]*W[i3*n+j]);
          }
          parta3=parta3*xw[i3];
          parta4=parta4*xw[i2];
          parta5=parta5*xw[i1];
          parta=(parta1+2.0*parta2-parta3-parta4-parta5)/l[i]/(l[i]-1.0)/(l[i]-2.0);
          m3(i1+p*i2+p*p*i3)=parta;
        }
      }
    }
    sm3=m3+sm3;
  }
  sm3=sm3/m;
  
  for(int j=0;j<p;j=j+1){
    s1=0.0;
    for (int k=0;k<n;k=k+1){
      s1=s1+m1[k]*W[j*n+k];
    }
    mw[j]=s1;
  }
  
  NumericVector wc(n*p);
  for(int j=0;j<p;j=j+1){
    for(int i=0;i<n;i=i+1){
       s1=0.0;
       for(int k=0;k<n;k=k+1){
        s1=s1+W[j*n+k]*C2[i*n+k];
      }
       wc[n*j+i]=s1;
    } 
  }   
  
  NumericVector wcw(p*p);
  for(int i=0;i<p;i=i+1){
    for(int j=0;j<p;j=j+1){
      s1=0.0;
      for(int k=0;k<n;k=k+1){
        s1=s1+wc[n*i+k]*W[j*n+k];
      }
      wcw[p*i+j]=s1;
    } 
  }           
  
  for(int i1=0;i1<p;i1=i1+1){
    for(int i2=0;i2<p;i2=i2+1){
      for(int i3=0;i3<p;i3=i3+1){
        partb=wcw[i1*p+i2]*mw[i3];
        partc=wcw[i1*p+i3]*mw[i2];
        partd=wcw[i2*p+i3]*mw[i1];
        parte=mw[i1]*mw[i2]*mw[i3];
        sm3(i1+p*i2+p*p*i3) = sm3(i1+p*i2+p*p*i3)-a0*(partb+partc+partd)/(a0+2.0)+2.0*pow(a0,2)*parte/(a0+1.0)/(a0+2.0);
      }
    }
  }          
  
  return(sm3);
}

// [[Rcpp::export]]
NumericVector Npowercpp(int N,NumericVector m3,NumericVector theta,int p){ 
  NumericVector  tta(p*p),otheta(p),stta(p);
  otheta=theta;
  for(int i=0;i<N;i=i+1){
    for(int i1=0;i1<p;i1=i1+1){
      for(int i2=0;i2<p;i2=i2+1){
        for(int i3=0;i3<p;i3=i3+1){
          tta[i3+i2*p]=m3[i1+i2*p+i3*p*p]*otheta[i2]*otheta[i3];
        }
      }
      stta[i1]=sum(tta);
    }
    theta=stta/sqrt(sum(pow(stta,2)));
  }
  return(theta);
}

// [[Rcpp::export]]
double lambdacpp(int p,NumericVector m3,NumericVector theta){
  NumericVector rq(p*p*p);
  for(int i1=0;i1<p;i1=i1+1){
    for(int i2=0;i2<p;i2=i2+1){
      for(int i3=0;i3<p;i3=i3+1){
        rq(i1+i2*p+i3*p*p)=m3(i1+i2*p+i3*p*p)*theta[i1]*theta[i2]*theta[i3];
      }
    }
  }
  return(sum(rq));
}