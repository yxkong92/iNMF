#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector Ksphere(const int K) {
  NumericVector theta;
  theta = rnorm(K);
  theta = theta/sqrt(sum(pow(theta,2)));
  return theta;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
Ksphere(2)
*/
