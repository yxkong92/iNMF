#include <Rcpp.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdlib.h>
using namespace Rcpp;


/*Bayes NO alpha*/
/*Mathematical constants computed using Wolfram Alpha */
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

//typedef enum { false, true } bool;

double aterm(int n, double x, double t)
{
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + log(n + 0.5) + 1.5*(MATH_LOG_2_PI-log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }
  return exp(f);
}

double exprnd(double mu)
{
  return -mu * log(1.0 - R::runif(0.0,1.0));
}

double truncgamma()
{
  double c = MATH_PI_2;
  double X, gX;

  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / sqrt(X);

    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }

  return X;
}

double randinvg(double mu)
{
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - sqrt(4*mu*V + mu*mu * V*V) );

  if(R::runif(0.0,1.0) > mu /(mu+out)) {
    out = mu*mu / out;
  }
  return out;
}

double tinvgauss(double z, double t)
{
  double X, u;
  double mu = 1.0/z;

  if(mu > t) {
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma();

      if(log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }
  else {
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }
  return X;
}

double samplepg(double z)
{
  z = fabs(z) * 0.5;
  double t = MATH_2_PI;

  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = log(4) - MATH_LOG_PI - z;
  double logK = log(K);
  double Kt = K * t;
  double w = sqrt(MATH_PI_2);

  double logf1 = logA +  R::pnorm(w*(t*z - 1),0.0,1.0,1,1)+ logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;

  double p_over_q = exp(logf1) + exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q);

  double u, X;

  while(1)
  {
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      X = t + exprnd(1.0)/K;
    }
    else {
      X = tinvgauss(z, t);
    }

    int i = 1;
    double Sn = aterm(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;

    while(1)
    {
      Sn = Sn + asgn * aterm(i, X, t);

      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }

      if(even && (U > Sn)) {
        break;
      }

      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}




void computeomega(NumericVector Kappa, NumericVector Tau, NumericVector probvec, int N,  int M, NumericVector result)
{
  int i;
  int j;

  int n = N;
  int m = M;

  for (i=0; i < n; i++){
    for(j=0; j < m; j++){
      result[j * n + i] = samplepg( Kappa[j]  + Tau[j] * probvec[j *n + i]) ;
    }
  }
}


void sampleKT(NumericVector Omega, NumericVector S, NumericVector probM, double muKappa, double sigmaKappa, double muTau, double sigmaTau, int N, int M, NumericVector output, NumericVector SOmega, NumericVector SOmegaP, NumericVector SOmegaPsq, NumericVector SS, NumericVector SSP){

  int i, j, k;
  int n = N;
  int m = M;

  double s1, s2, s3, s4, s5, v11, v12, v22, denom, chol11, chol12, chol22, a, b, c, d, m1, m2;

  //  double SOmega[n * m],  SOmegaP[ n*m ],  SOmegaPsq [n*m];
  //double SS[ n*m ], SSP[n*m ];


  for(j=0; j <(m * n); j++){

    SOmega[j] = 0;
    SOmegaP[j] = 0;
    SOmegaPsq[j] = 0;
    SS[j] = 0;
    SSP[j] = 0;

  }

  for(j =0; j < m; j++){

    s1 = 0.0;
    s2 = 0.0;
    s3 = 0.0;
    s4 = 0.0;
    s5 = 0.0;

    for( i=0;  i < n; i++){
      k = j*n + i;
      s1 =   s1 + Omega[k];
      s2 =   s2 + Omega[k] * probM[k];
      s3 =   s3 + Omega[k] * pow(probM[k], 2);

      s4 =   s4  + S[k];
      s5 =   s5 +  S[k] * probM[k];
    }

    SOmega[j]  = s1;
    SOmegaP[j] = s2;
    SOmegaPsq[j] = s3;

    SS[j] = s4;
    SSP[j] = s5;

  }

  for(j=0; j < m; j++){
    v11  = SOmega[j] +  1.0 / pow(sigmaKappa, 2);
    v12  =  SOmegaP[j];
    v22  = SOmegaPsq[j] + 1.0/pow(sigmaTau, 2);

    s1 = SS[j] - (double)n / 2.0 + muKappa / pow(sigmaKappa,2);
    s2 = SSP[j] - 0.5  + muTau / pow(sigmaTau, 2);

    denom = v11 * v22 - v12 * v12;
    m1 =( v22 * s1 - v12 * s2) / denom;
    m2 = (v11 * s2 - v12 * s1) / denom;

    chol11 = sqrt(v22 / denom);
    chol12 =  -v12 / denom / chol11 ;
    chol22 = sqrt(v11 / denom - pow(chol12, 2) );

    a = R::rnorm(0.0, 1.0);
    b = R::rnorm(0.0, 1.0);

    c = chol11 * a  + m1;
    d = chol12 * a + chol22 * b +m2;

    output[j] = c;
    output[m +j] = d;

  }

}



void sampleS(NumericVector Words, NumericVector probM, NumericVector S, NumericVector Kappa, NumericVector Tau, NumericVector Total, int N, int M, NumericVector output){

  int n = N;
  int m = M;

  int i, j, k, l;
  double tot, success, successprob;

  for( j = 0; j < m ; j++){
    for( i =0; i < n ; i++){
      output[n*j + i] =S[n*j +i];
    }
  }

  for( j = 0; j < m ; j++){
    for( i =0; i < n ; i++){
      k = j * n + i;
      if(Words[k] > 0){
        output[k] = 1;
      }else{
        tot = 0;
        for( l = 0; l < n; l++){
          if(l != i){
            tot =  tot + output[j * n + l ] *  probM[ j* n + l];
          }
        }
        success = Kappa[j] + Tau[j] * probM[k] + Total[j]*log(tot / (probM[k] + tot)) ;
        successprob = 1.0 / (1.0 + exp(-success));
        output[k]  = R::rbinom(1, successprob);
      }
    }
  }
}


void test( NumericVector Av, NumericVector Wv, NumericVector Xv,  NumericVector Total,  NumericVector eqS, NumericVector EqSKappa, NumericVector EqSTau, NumericVector EqKappa, NumericVector EqTau, NumericVector EqOmegaKappasq, NumericVector EqOmegaKappaTau, NumericVector EqOmegaTausq, int N, int M, int M1, int K, NumericVector output){

  int m = M;
  int n = N;
  int kk = K;

  int i, j, k;
  int count = 1;

  double tot = 0.0, tot2=0.0;
  double s1, s2;
  double diff1 = 1.0;

  //double  prob[n * m], A[ n * m ],  D[m];
  // double r1[n * kk], r2[m*kk];

  //Test
  Rcpp::NumericVector A(n * m),D(m),r2(m*kk), r1(n*kk),prob(n * m);

  while( (diff1 > 5e-4) & (count < 500)){
    /* probability update */
    for(i=0; i < n; i++){
      for(j=0; j < m; j++){
        tot = 0;
        for(k=0; k < kk; k++){
          tot = tot + Av[n * k + i] * Wv[m * k + j];
        }
        /* added*/
        if(tot < 1e-32){
          tot = 1e-32;
        }
        //
        prob[ j * n +i] = tot;
      }
    }

    for(j = 0; j < m; j++){
      for(i=0; i < n; i++){
        A[j * n + i] =  eqS[j * n + i] * Xv[j *n +i] / prob[j * n + i];
      }
    }

    for(j = 0; j < m; j++){
      tot = 0;
      for(i = 0; i < n; i++){
        tot = tot + eqS[j * n +i] * prob[ j * n + i];
      }

      /* added*/
      if(tot < 1e-32){
        tot = 1e-32;
      }
      //

      D[j] =  Total[j] / tot;
    }

    for(k = 0; k < kk; k++){
      for(i=0; i < n; i++){
        tot = 0.0;
        tot2 = 0.0;
        for(j=0; j < m; j++){
          tot = tot + ( A[j*n +i] + EqSTau[j*n+i] - 0.5*EqTau[j] -EqOmegaKappaTau[j*n+i] - EqOmegaTausq[j*n+i]*prob[j*n+i])*Wv[ k*m +j];
          tot2 = tot2 +  eqS[j*n+i] * D[j]*Wv[k*m + j];
        }

        /* added*/
        if(tot2 < 1e-32){
          tot2 = 1e-32;
        }
        //

        r1[n*k+i] = Av[n*k+i]*tot/tot2;
      }
    }

    for(k=0; k < kk; k++){
      tot = 0;
      for(i = 0; i < n; i++){
        tot  = tot + r1[n*k + i];
      }

      for(i = 0; i < n ; i++){
        r1[n*k +i] =  r1[n*k + i] / tot;
      }
    }

    s1  = 0.0;
    s2 = 0.0;
    for(k = 0; k < kk; k++){
      for(i =0; i < n;  i++){
        s1 = s1 +  pow(Av[ n * k + i] - r1[n * k + i], 2);
        s2 = s2 +  pow(Av[n*k + i], 2);
        Av[n*k + i] = r1[n*k +i];
      }
    }

    diff1 = sqrt(s1 / s2);

    count++;
  }



  count = 1;
  diff1 = 1;
  while( (diff1 > 5e-4) & (count < 500)){

    for(i=0; i < n; i++){
      for(j=0; j < m; j++){
        tot = 0;
        for(k=0; k < kk; k++){
          tot = tot + Av[n * k + i] * Wv[m * k + j];
        }
        /* added*/
        if(tot < 1e-32){
          tot = 1e-32;
        }
        //
        prob[ j * n +i] = tot;
      }
    }


    for(j = 0; j < m; j++){
      for(i=0; i < n; i++){
        A[j * n + i] =  eqS[j * n + i] * Xv[j *n +i] / prob[j * n + i];
      }
    }

    for(j = 0; j < m; j++){
      tot = 0;
      for(i = 0; i < n; i++){
        tot = tot + eqS[j * n +i] * prob[ j * n + i];
      }
      /* added*/
      if(tot < 1e-32){
        tot = 1e-32;
      }
      //


      D[j] =  Total[j] / tot;
    }

    for(k = 0; k < kk; k++){
      for(j=0; j < (M1); j++){
        tot = 0.0;
        tot2 = 0.0;
        for(i=0; i < n; i++){
          tot = tot + ( A[j*n +i] +  EqSTau[j*n+i] - 0.5*EqTau[j] -EqOmegaKappaTau[j*n+i] - EqOmegaTausq[j*n+i]*prob[j*n+i])*Av[ k*n +i];
          tot2 = tot2 +  eqS[j*n+i] * D[j]*Av[k*n + i];
        }
        /* added*/
        if(tot2 < 1e-32){
          tot2 = 1e-32;
        }
        //


        r2[m*k+j] = Wv[m*k+j]*tot/tot2;
      }
    }
    for(j=0; j < (M1); j++){
      tot = 0;
      for(k = 0; k < kk; k++){
        tot  = tot + r2[m*k + j];
      }

      for(k= 0; k < kk ; k++){
        r2[m*k +j] =  r2[m*k + j] / tot;
      }
    }

    s1  = 0.0;
    s2 = 0.0;
    for(k = 0; k < kk; k++){
      for(j=0; j < (M1);  j++){
        s1 = s1 +  pow(Wv[ m * k + j] - r2[m * k + j], 2);
        s2 = s2 +  pow(Wv[m *k + j], 2);
        Wv[ k*m +j ] = r2[m*k + j];
      }
    }

    diff1 = sqrt(s1 / s2);

    count++;
  }

  /* store output*/
  for(k =0; k < kk;  k++){
    for(i =0; i < n;  i++){
      output[n*k +i] =r1[n*k + i];
    }
    for(j =0; j < M1;  j++){
      output[n*kk +m*k+j] =r2[m*k + j];
    }
  }
}


double lfunction(NumericVector Av,NumericVector Wv, NumericVector Xv, NumericVector Total, NumericVector EqS, NumericVector EqSKappa, NumericVector EqSTau, NumericVector EqKappa, NumericVector EqTau, NumericVector EqOmegaKappasq,NumericVector EqOmegaKappaTau,NumericVector EqOmegaTausq, int N, int M, int K, double sdKappa, double sdTau){

  int n = N, m = M, kk=K;
  double tot, a, Term3[m];//, prob[n*m];
  int i, j, k;

  NumericVector prob(n * m);

  for(j = 0; j < m; j++){
    for(i = 0; i < n; i++){
      tot =  0.0;
      for(k = 0; k < kk; k++){
        tot = tot + Av[k*n + i]*Wv[k*m+j];
      }
      prob[j*n + i] = tot;
    }
  }

  for(j = 0; j < m; j++){
    Term3[j] = 0.0;
    for(i=0; i < n; i++){
      Term3[j] = Term3[j] + EqS[n*j + i]*prob[n*j +i];
    }
  }

  tot = 0.0;
  for(j = 0; j < m; j++){
    for(i = 0; i < n; i++){
      k = j*n +i;
      tot = tot + EqS[k]*Xv[k]*log(prob[k]) + EqSKappa[k] - 0.5*EqKappa[j] +  EqSTau[k]*prob[k] - 0.5* EqTau[j]*prob[k] - 0.5*EqOmegaKappasq[k] -EqOmegaKappaTau[k]*prob[k] - 0.5*EqOmegaTausq[k]*pow(prob[k], 2.0);
    }
    tot = tot - Total[j]*(n+log(Term3[j]));
  }

  a = (double) tot /n / m  - (double) (log(sdKappa) + log(sdTau))/ n;
  return a;
}






// [[Rcpp::export]]
Rcpp::NumericVector MCMCsample2(NumericVector Xv, NumericVector Total, int N, int M, int M1, int K, NumericVector Av, NumericVector Wv, int Mccount, int burnin, NumericVector EqS, NumericVector EqSKappa, NumericVector EqSTau, NumericVector EqKappa,  NumericVector EqKappasq, NumericVector EqTau, NumericVector EqTausq, NumericVector EqOmegaKappasq, NumericVector EqOmegaKappaTau, NumericVector EqOmegaTausq){
  double thres=5e-3,  diff1=1.0,  tot=0.0;
  double sdKappa  = 100.0, sdTau = 100.0, muKappa = 0.0, muTau = 0.0;
  int n = N, m = M, kk = K, L =  Mccount - burnin + 1, miter, i, j, k, outer=1;
  //  double prob[n*m], Omega[n*m], S[n*m], Kappa[m], Tau[m], KT[2*m], nS[n*m];
  //NumericVector EqS(n*m), EqSKappa(n*m), EqSTau(n*m), EqKappa(m),EqKappasq(m), EqTau(m), EqTausq(m), EqOmegaKappasq(n*m), EqOmegaKappaTau(n*m), EqOmegaTausq(n*m);
  // double output[(n+m)*kk] ;
  double ob, oa;

  NumericVector prob(n * m),Omega(n * m),S(n * m),nS(n * m),output((n+m)*kk),output2((n+m)*kk +1);

  NumericVector Kappa(m),Tau(m),KT(2*m);

  //Sample KT
  NumericVector SOmega(n * m),SOmegaP(n * m),SOmegaPsq(n * m),SS(n * m),SSP(n * m);

  for(k = 0; k < kk; k++){
    for(i = 0; i < n; i++){
      output[n*k + i] = Av[n*k+i];
    }
    for(j = 0; j < m; j++){
      output[m*k + j + n*kk] = Wv[m*k + j];
    }
  }


  /*initialize Kappa and Tau, S*/
  for(j=0; j < m; j++){
    Kappa[j] = R::rnorm(muKappa, sdKappa);
    Tau[j] = R::rnorm(muTau, sdTau);
    for(i =0; i < n; i++){
      if(Xv[j*n+i] > 0){
        S[j*n+i] =1;
      }else{
        S[j*n+i] = R::rbinom(1,0.5);
      }
    }
  }

  for(i=0; i < n*m; i++){
    Omega[i] = 0.0;
    nS[i] = 0.0;
  }
  for(i=0; i < 2*m; i++) KT[i] = 0.0;

  oa = 10000000;
  while((diff1 > thres)  & (outer  < 50)){
    /*Gibbs sampling*/
    for(j = 0; j < m; j++){
      for(i = 0; i < n; i++){
        tot =  0.0;
        for(k = 0; k < kk; k++){
          tot = tot + Av[k*n + i]*Wv[k*m+j];
        }
        prob[j*n + i] = tot;
      }
    }

    for(miter = 0; miter < burnin-1; miter++){
      computeomega(Kappa, Tau, prob, N,  M, Omega);
      sampleKT(Omega, S, prob, muKappa, sdKappa,  muTau, sdTau, N, M, KT, SOmega, SOmegaP, SOmegaPsq, SS, SSP);
      for(j = 0; j < m; j ++){
        Kappa[j] = KT[j];
        Tau[j] = KT[m+j];
      }
      sampleS(Xv, prob, S, Kappa, Tau, Total, N, M, nS);
      for(i = 0 ;  i < n*m; i++) S[i] = nS[i];
    }

    for(j = 0; j < m; j++){
      EqKappa[j] = 0.0;
      EqTau[j] = 0.0;
      EqKappasq[j] = 0.0;
      EqTausq[j] = 0.0;
      for(i = 0; i < n; i++){
        k = j*n + i;
        EqS[k] = 0.0;
        EqSKappa[k] = 0.0;
        EqSTau[k] = 0.0;
        EqOmegaKappasq[k] = 0.0;
        EqOmegaKappaTau[k] = 0.0;
        EqOmegaTausq[k] = 0.0;
      }
    }

    for(miter = 0; miter < L; miter++){
      computeomega(Kappa, Tau, prob, N,  M, Omega );
      sampleKT(Omega, S, prob, muKappa, sdKappa,  muTau, sdTau, N, M, KT, SOmega, SOmegaP, SOmegaPsq, SS, SSP);

      for(j = 0; j < m; j ++){
        Kappa[j] = KT[j];
        Tau[j] = KT[m+j];
      }
      sampleS(Xv, prob, S, Kappa, Tau, Total, N, M, nS);
      for(i = 0 ;  i < n*m; i++) S[i] = nS[i];

      /*compute expectations*/
      for(j = 0; j < m; j++){
        EqKappa[j] = EqKappa[j] + Kappa[j];
        EqTau[j] = EqTau[j] + Tau[j];
        EqKappasq[j] = EqKappasq[j] + pow(Kappa[j], 2.0);
        EqTausq[j] = EqTausq[j] + pow(Tau[j], 2.0);
        for(i = 0; i <n; i++){
          k = j*n + i;
          EqS[k] = EqS[k] + S[k];
          EqSKappa[k] = EqSKappa[k] + S[k]*Kappa[j];
          EqSTau[k] = EqSTau[k] + S[k]*Tau[j];
          EqOmegaKappasq[k] = EqOmegaKappasq[k] + Omega[k]*pow(Kappa[j], 2.0);
          EqOmegaKappaTau[k] = EqOmegaKappaTau[k] + Omega[k]*Kappa[j]*Tau[j];
          EqOmegaTausq[k] = EqOmegaTausq[k] + Omega[k]*pow(Tau[j], 2.0);
        }
      }
    }

    for(j = 0; j < m; j++){
      EqKappa[j] = (double) EqKappa[j]/L;
      EqTau[j] = (double) EqTau[j]/L;
      EqKappasq[j] = (double) EqKappasq[j]/L;
      EqTausq[j] = (double) EqTausq[j]/L;
      for(i = 0; i < n; i++){
        k = j*n + i;
        EqS[k] = (double)EqS[k]/L;
        EqSKappa[k] =  (double) EqSKappa[k]/L;
        EqSTau[k] = (double) EqSTau[k]/L;
        EqOmegaKappasq[k] = (double) EqOmegaKappasq[k]/L;
        EqOmegaKappaTau[k] = (double) EqOmegaKappaTau[k]/L;
        EqOmegaTausq[k] = (double) EqOmegaTausq[k]/L;
      }
    }

    /*update the parameter based on the Gibbs samples */

    muKappa = 0.0;
    muTau = 0.0;
    sdKappa = 0.0;
    sdTau = 0.0;
    for(j = 0; j < m; j++){
      muKappa = muKappa +  EqKappa[j];
      muTau = muTau + EqTau[j];
      sdKappa = sdKappa + EqKappasq[j];
      sdTau = sdTau + EqTausq[j];
    }
    muKappa = (double) muKappa/m;
    muTau =  (double) muTau/m;
    sdKappa = sqrt((double)sdKappa/m - pow(muKappa, 2.0));
    sdTau = sqrt( (double)sdTau/m- pow(muTau, 2.0));


    //compute update
    test(Av, Wv,  Xv, Total, EqS, EqSKappa, EqSTau, EqKappa, EqTau,EqOmegaKappasq,EqOmegaKappaTau,EqOmegaTausq, N, M,M1, K, output);


    /*update*/

    for(k = 0; k < kk; k++){
      for(i=0; i < n; i++){
        Av[n*k + i] = output[n*k + i];
      }
      for(j = 0; j < M1; j++){
        Wv[m*k +j] = output[n*kk + m*k + j ];
      }
    }

    ob = lfunction(Av, Wv, Xv, Total,  EqS, EqSKappa, EqSTau, EqKappa, EqTau, EqOmegaKappasq, EqOmegaKappaTau, EqOmegaTausq, N, M,K, sdKappa, sdTau);

    if(outer == 1){
      diff1 = 1.0;
    }else{
      diff1 =  sqrt(pow(oa-ob, 2))/sqrt(pow(oa, 2));
    }


    oa = ob;
    outer ++;
  }

  for(k = 0; k < kk; k++){
    for(i=0; i < n; i++){
      output2[k*n + i] = output[k*n+i];
    }
    for(j=0; j < m; j++){
      output2[k*m + j + kk*n] = output[k*m+j + kk*n];
    }
  }
  output2[(m+n)*kk] = oa;
  return output2;
}





