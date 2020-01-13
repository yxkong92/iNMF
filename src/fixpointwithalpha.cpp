#include <Rcpp.h>
using namespace Rcpp;

double truncatedlikelihood(NumericVector Av, NumericVector Wv, NumericVector Xv, NumericVector lambda, int N, int M, int K, NumericVector alpha){
  int n = N;
  int m = M;
  int kk = K;
  int i,j,k;
  double tot, salpha;
  NumericVector entry(n*m), prob(n*m);
  
  for(j = 0; j < m; j++){
    for(i = 0; i < n; i++){
      if(Xv[j*n + i] == 0){
        entry[j*n + i] = 0.0;
      }else{
        entry[j*n + i] = 1.0;
      }
    }
  }
  
  
  for(i = 0; i < n; i++){
    for(j=0; j < m; j++){
      tot = 0.0;
      for(k = 0; k < kk; k++){
        tot  =  tot + Av[k*n + i ]*Wv[k*m +j];
      }
      //updated
      if(tot < 1e-16){
        tot = 1e-16;
      }
      
      prob[j*n + i] = tot;
    }
  }
  
  tot = 0.0;
  for( j = 0; j < m; j ++){
    for( i = 0; i < n; i++){
      if(entry[j*n +i] > 0){
        tot  = tot +  Xv[j*n + i]*log( lambda[j]*prob[j*n + i] ) -lgamma(Xv[j*n+i]+1) -lambda[j]*prob[j*n+i]- log(1.0 + exp(-lambda[j]*prob[j*n+i]));
      }
    }
  }
  
  for(k = 0; k < kk; k++){
    for(j = 0; j < m; j++){
      if(Wv[m*k+j] < 1e-16){
        Wv[m*k+j] = 1e-16;
      }
    }
  }
  
  salpha = 0.0;
  for( k= 0; k< kk; k++){
    salpha = salpha + alpha[k];
    for(j = 0; j < m; j++){
      tot =  tot - lgamma(alpha[k]) + (alpha[k]-1.0)*log(Wv[m*k+j]);
    }
  }
  
  tot = (double)tot/n/m + lgamma(salpha)/n/m;
  return tot;
}

// [[Rcpp::export]]
Rcpp::NumericVector fitfixedpointAlpha(NumericVector Xv, NumericVector lambda, int N, int M, int M1, NumericVector Av, int updateA,  NumericVector Wv, int K, NumericVector alpha){
  double idiff = 1.0, thres = 1e-3, diffout = 1.0, diff = 1.0, tot, tot1, tot2, toto, oa, ob;
  const int MAXITER = 500;
  int countout = 1, i, j, k, counter;
  int n = N;
  int m = M;
  int kk = K;
  double entry[n*m], prob[n*m], Wnv[m*kk], Anv[n*kk];
  //Rcpp::NumericVector output((n+m)*kk+1);
  Rcpp::NumericVector output((n+m)*kk+1),r2(m*kk), r1(n*kk);
  //NumericVector entry, r2, r1, prob, Wnv, Anv,output;
  
  //Any row should not be completely zero! It should be checked in the R.
  //Preparation start.
  for(j = 0; j < m; j++){
    for(i = 0; i < n; i++){
      if(Xv[j*n + i] == 0.0){
        entry[j*n + i] = 0.0;
      }else{
        entry[j*n + i] = 1.0;
      }
    }
  }
  
  //r1 and r2 are the coefficients stored while Anv and Wnv are the temperary variables.
  for(k =0; k < kk; k++){
    for(i=0; i < n; i++){
      r1[k*n + i] = Av[k*n+i];
      Anv[k*n + i] = 0.0;
    }
    for(j = 0; j < m; j++){
      r2[k*m + j] = Wv[k*m + j];
      Wnv[k*m+j] = 0.0;
    }
  }
  
  //evaluate objective function at the initial
  oa =  truncatedlikelihood(Av, Wv,Xv, lambda, N, M, K, alpha);
  
  //update W and A
  while( (diffout > thres)  & (countout < MAXITER)){
    //update W
    //if M1 == 0, it means that there is no need to updating W
    if(M1 > 0 ){
      counter = 1;
      idiff = 1.0;
      diff = 1.0;
      while( ((idiff > thres) | (diff*10 > thres)) & (counter < MAXITER)){
        // compute prob = AW
        for(i = 0; i < n; i++){
          for(j = 0; j < m; j++){
            tot = 0.0;
            for(k = 0; k < kk; k++){
              tot =  tot + r1[k*n + i]*r2[k*m +j];
            }
            //zero cells are set to be a tiny value
            if(tot < 1e-16){
              tot = 1e-16;
            }
            prob[j*n + i] = tot;
          }
        }
        
        //compute the sums to be used for the numerator and denomenator
        for(k =0; k < kk; k++){
          for(j=0; j < M1; j++){
            tot1 = 0.0;
            tot2 = 0.0;
            for(i = 0; i < n; i++){
              if(entry[n*j +i] > 0){
                tot1 = tot1 + Xv[n*j +i]*r1[n*k +i]/(prob[n*j + i]);
                tot2 = tot2 + r1[n*k +i]*lambda[j]/(1.0 - exp( - lambda[j]*prob[n*j + i]));
              }
            }
            
            if(tot2 < 1e-16){
              tot2 = 1e-16;
            }
            
            Wnv[k*m + j] = (double) (r2[k*m + j] * tot1 + alpha[k] - 1.0)/tot2 ;
          }
        }
        
        
        //normalize columns of W
        toto = 0.0;
        diff = 0.0;
        for(j = 0; j < M1; j++){
          tot = 0.0;
          for(k = 0; k < kk; k++){
            tot =  Wnv[k*m+j] + tot;
            toto = toto + pow(r2[k*m+j],2);
          }
          
          //The following is not needed if things are fine.
          if(tot < 1e-16){
            tot = 1e-16;
          }
          
          for(k = 0; k < kk; k++){
            Wnv[k*m+j] = Wnv[k*m+j]/tot;
            diff = diff + pow(Wnv[k*m+j] - r2[k*m+j], 2);
            r2[k*m + j] = Wnv[k*m+j];
          }
        }
        if(toto < 1e-16){
          toto = 1e-16;
        }
        idiff = sqrt(diff/toto);
        counter++;
      }
    }
    
    //update A
    //A needs to be updated? updateA should be 0
    if(updateA!=1){
      counter = 1;
      idiff = 1.0;
      diff = 1.0;
      while( ((idiff > thres) | (diff*10 > thres)) & (counter < MAXITER)){
        for(i = 0; i < n; i++){
          for(j = 0; j < m; j++){
            tot = 0.0;
            for(k = 0; k < kk; k++){
              tot =  tot + r1[k*n + i]*r2[k*m +j];
            }
            if(tot < 1e-16){
              tot = 1e-16;
            }
            prob[j*n + i] = tot;
          }
        }
        
        for(k = 0; k < kk; k++){
          for(i = 0; i < n; i++){
            tot1 = 0.0;
            tot2 = 0.0;
            for(j = 0; j < m; j++){
              if(entry[n*j + i] > 0){
                tot1 = tot1 +  Xv[n*j +i]* r2[m*k +j]/(prob[n*j + i]);
                tot2 = tot2 +  r2[m*k +j]*lambda[j]/(1.0 - exp( - lambda[j]*prob[n*j + i]));
              }
            }
            
            if(tot2 < 1e-16){
              tot2 = 1e-16;
            }
            
            Anv[k*n + i] = r1[k*n + i] * tot1/tot2;
          }
        }
        
        /*normalize rows of A*/
        diff = 0.0;;
        toto = 0.0;
        for(k = 0; k < kk; k++){
          tot = 0.0;
          for(i = 0; i < n; i++){
            tot = tot + Anv[k*n+i];
            toto = toto + pow(r1[k*n+i],2);
          }
          //updated
          
          if(tot < 1e-16){
            tot = 1e-16;
          }
          
          for(i = 0; i < n; i++){
            Anv[k*n+i] = Anv[k*n+i]/tot;
            diff =  diff + pow(Anv[k*n+i] - r1[k*n+i], 2);
            r1[k*n + i] = Anv[k*n+i];
          }
        }
        
        if(toto < 1e-16){
          toto = 1e-16;
        }
        
        idiff = sqrt(diff/toto);
        counter++;
      }
    }
    /*function evaluation*/
    ob =  truncatedlikelihood(r1,r2, Xv, lambda, N, M, K, alpha );
    diffout = sqrt(pow(oa - ob,2));
    oa = ob;
    countout++;
  }
  
  for(k = 0; k < kk; k++){
    for(i = 0;  i < n; i++){
      output[k*n + i] =  r1[k*n+i];
    }
  }
  
  for(k = 0; k < kk; k++){
    for(j = 0;  j < m; j++){
      output[k*m + j + kk*n] =  r2[k*m+j];
    }
  }
  output[(n+m)*kk] = oa;
  return output;
}