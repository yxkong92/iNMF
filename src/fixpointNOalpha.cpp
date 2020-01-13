#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double truncatedlikelihoodnoalpha(NumericVector Av, NumericVector Wv, NumericVector Xv, NumericVector lambda, int N, int M, int K){
  int n = N;
  int m = M;
  int kk = K;
  int i,j,k;
  double tot;
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
  
  tot = (double) tot/n/m;  
  return tot;
  
}

// [[Rcpp::export]]
Rcpp::NumericVector fitfixedpoint(NumericVector Xv, NumericVector lambda, int N, int M, int M1, NumericVector Av, NumericVector Wv, int K){
  double idiff = 1.0, thres = 5e-3, diffout = 1.0, diff = 1.0, tot, tot1, tot2, toto, oa, ob;
  int countout = 1, i, j, k, counter;
  int n = N;
  int m = M;
  int kk = K;
  double entry[n*m], prob[n*m], Wnv[m*kk], Anv[n*kk];
  NumericVector output((n+m)*kk+1),r2(m*kk), r1(n*kk); 
  
  for(j = 0; j < m; j++){
    for(i = 0; i < n; i++){
      if(Xv[j*n + i] == 0.0){
        entry[j*n + i] = 0.0;  
      }else{
        entry[j*n + i] = 1.0;
      }
    }
  }
  
  for(k =0; k < kk; k++){
    for(i=0; i < n; i++){
      r1[k*n + i] = Av[k*n+i];
    }
    for(j = 0; j < m; j++){
      r2[k*m + j] = Wv[k*m + j];
    }
  }
  
  /*evaluate objective function at the initial */
  
  oa =  truncatedlikelihoodnoalpha(Av, Wv,Xv, lambda, N, M, K );
  
  /*update W and A */
  
  while( (diffout > thres)  & (countout < 1000)){
    
    /*update W*/
    counter = 1; 
    idiff = 1.0;
    while( (idiff > thres) & (counter < 1000)){
      
      for(i = 0; i < n; i++){
        for(j = 0; j < m; j++){
          tot = 0.0;
          for(k = 0; k < kk; k++){
            tot =  tot + r1[k*n + i]*r2[k*m +j];
          }
          //updated
          
          if(tot < 1e-32){
            tot = 1e-32;
          }
          
          prob[j*n + i] = tot; 
        }
      }
      
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
          
          if(tot2 < 1e-32){
            tot2 = 1e-32;
          }
          //   printf("%f", tot2);
          //printf("\n");
          Wnv[k*m + j] = r2[k*m + j] * tot1/tot2;
        }
      }
      
      /*normalize rows of W */              
      toto = 0.0;
      diff = 0.0;
      for(j = 0; j < M1; j++){
        tot = 0.0;
        for(k = 0; k < kk; k++){
          tot = tot + Wnv[k*m+j];
          toto = toto + pow(r2[k*m+j],2);
        }
        
        //updated
        
        if(tot < 1e-32){
          tot = 1e-32;
        }
        
        for(k = 0; k < kk; k++){
          Wnv[k*m+j] = Wnv[k*m+j]/tot;
          diff = diff + pow(Wnv[k*m+j] - r2[k*m+j], 2);
          r2[k*m + j] = Wnv[k*m+j];
        }
      }
      idiff = sqrt(diff/toto);  
      counter++;
    }
    
    
    /*update A*/
    counter = 1; 
    idiff = 1.0;
    while( (idiff > thres) & (counter < 1000)){
      
      
      
      for(i = 0; i < n; i++){
        for(j = 0; j < m; j++){
          tot = 0.0;
          for(k = 0; k < kk; k++){
            tot =  tot + r1[k*n + i]*r2[k*m +j];
          }
          
          //updated
          
          if(tot < 1e-32){
            tot = 1e-32;
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
          
          
          if(tot2 <  1e-32){
            tot2 = 1e-32;
          }
          Anv[k*n + i] = r1[k*n + i] * tot1/tot2;
          
          //      printf("%f", tot2);
          //printf("\n");
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
        
        if(tot < 1e-32){
          tot = 1e-32;
        }
        
        for(i = 0; i < n; i++){
          Anv[k*n+i] = Anv[k*n+i]/tot;
          diff =  diff + pow(Anv[k*n+i] - r1[k*n+i], 2);
          r1[k*n + i] = Anv[k*n+i];
        } 
      }	 
      idiff = sqrt(diff/toto); 
      counter++;	 
    }
    
    /*function evaluation*/
    
    ob =  truncatedlikelihoodnoalpha(r1, r2, Xv, lambda, N, M, K );
    
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