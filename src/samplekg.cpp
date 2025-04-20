#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//----------------------------------------------------------------
// Mathematical constants
//----------------------------------------------------------------
#define MATH_PI    3.1415926535897932384626433832795028841971693993751
#define MATH_PI2   (MATH_PI * MATH_PI)
#define TRUNC      0.050239

//----------------------------------------------------------------
// Function: a_coef
//
// Computes the coefficient a_n(x) 
//----------------------------------------------------------------
double a_coef(int n, double xval) {
  if(xval > TRUNC) {
    double relval = std::pow(MATH_PI, 2) * std::pow(n + 1, 2);
    return 4.0 * relval * std::exp(-2.0 * relval * xval);
  } else {
    double relval = std::pow(n + 1, 2);
    if(std::isnan(relval)) Rcpp::stop("NaN encountered in a_coef");
    if(n % 2 == 0) { // even
      return std::pow(2.0 * MATH_PI, -0.5) / 4.0 *
        std::pow(xval, -2.5) * relval *
        std::exp(-relval / (8.0 * xval));
    } else { // odd
      return std::pow(2.0 * MATH_PI, -0.5) *
        std::pow(xval, -1.5) *
        std::exp(- (n * n) / (8.0 * xval));
    }
  }
}



//----------------------------------------------------------------
// Function: pinvgauss
//
// inverse Gaussian CDF (GIG with p = -1/2)
//----------------------------------------------------------------
double pinvgauss(double x, double a, double b) {
  double term1 = R::pnorm(std::sqrt(b/x) * (x * std::sqrt(a/b) - 1.0), 0.0, 1.0, 1, 0);
  double term2 = std::exp(2.0 * std::sqrt(a * b)) *
    R::pnorm(-std::sqrt(b/x) * (x * std::sqrt(a/b) + 1.0), 0.0, 1.0, 1, 0);
  return term1 + term2;
}

// Exported version
// [[Rcpp::export]]
double pinvgauss_export(double x, double a, double b) {
  return pinvgauss(x, a, b);
}

//----------------------------------------------------------------
// Function: pgig_onehalf
//
// computes pgig(x, 1/2, a, b)
//----------------------------------------------------------------
double pgig_onehalf(double x, double a, double b) {
  double k_one_half = std::sqrt((MATH_PI / 2.0) / std::sqrt(a * b)) *
    std::exp(-std::sqrt(a * b));
  double w = 1.0 / (1.0 + std::pow(a * b, -0.5));
  double k_three_halves = k_one_half / w;
  double gammainc_half = std::sqrt(MATH_PI) * std::erfc(std::sqrt(x));
  double gammainc_neg_half = -2.0 * gammainc_half +
    2.0 / std::sqrt(b / (2.0 * x)) *
    std::exp(-b / (2.0 * x));
  double ival = std::pow(a * b, -0.25) / std::sqrt(2.0) / k_three_halves * gammainc_half +
    (1.0 - w) * std::pow(a * b, 0.25) / std::pow(2.0, 1.5) / k_one_half * gammainc_neg_half;
  return w * pinvgauss(x, a, b) +
    (1.0 - w) * (1.0 - pinvgauss(1.0 / x, b, a)) -
    std::exp(-a * x / 2.0) * ival;
}

//----------------------------------------------------------------
// Function: pgig
//
// Computes the GIG CDF (with p -3/2).
//----------------------------------------------------------------
double pgig_negonehalf(double x, double a, double b) {
  if(x == 0.0) return 0.0;
  return 1.0 - pgig_onehalf(1.0 / x, b, a);
}

// Exported version
// [[Rcpp::export]]
double pgig_onehalf_export(double x, double a, double b) {
  return pgig_onehalf(x, a, b);
}
// Exported version
// [[Rcpp::export]]
double pgig_negonehalf_export(double x, double a, double b) {
  return pgig_negonehalf(x, a, b);
}

//----------------------------------------------------------------
// Function: log_sum_exp
//
// Numerically stable computation of log(exp(a)+exp(b)).
//----------------------------------------------------------------
double log_sum_exp(double a, double b) {
  double maxab = std::max(a, b);
  double neg_diffab = -std::fabs(a - b);
  return maxab + std::log(1.0 + std::exp(neg_diffab));
}

// Forward declaration for rgig_half
double rgig_negonehalf(double a, double b);

// Forward declaration
double rkg_one_zero();

//----------------------------------------------------------------
// Function: rkg_one
//
// Samples one value from the target distribution.
//----------------------------------------------------------------
double rkg_one(double z) {
  if(std::fabs(z) < 1e-6) {
    return rkg_one_zero();
  }
  z = std::fabs(z);
  double logp = std::log(z + 2.0) - z / 2.0 +
    std::log(pgig_negonehalf(TRUNC, z * z, 0.25));
  double k = 2.0 * MATH_PI2 + (z * z) / 2.0;
  double logq = std::log(4.0 * MATH_PI2) - k * TRUNC - std::log(k);
  double q_over_p_plus_q = (z < 200.0) ? std::exp(logq - log_sum_exp(logp, logq))
    : 0.0;

  double propval = 0.0;
  while(true) {
    if(R::runif(0.0, 1.0) < q_over_p_plus_q) {
      // Truncated exponential branch
      double gamma_sample = R::rgamma(1.0, 1.0); // Gamma(shape=1, scale=1)
      propval = TRUNC + gamma_sample / k;
      if(std::isnan(propval) || propval <= 0.0)
        Rcpp::stop("Invalid propval in rkg_one (right branch)");
    } else {
      propval = 1.0;
      while(propval > TRUNC) {
        propval = rgig_negonehalf(z * z, 0.25);
      }
      if(std::isnan(propval) || propval <= 0.0)
        Rcpp::stop("Invalid propval in rkg_one (left branch)");
    }

    double s = a_coef(0, propval);
    if(std::isnan(s)) Rcpp::stop("s is NaN in rkg_one");
    double y = R::runif(0.0, 1.0) * s;
    int n = 0;
    while(true) {
      n++;
      if(n % 2 == 1) { // odd term
        s -= a_coef(n, propval);
        if(y <= s) break;
      } else {         // even term
        s += a_coef(n, propval);
        if(y > s) break;
      }
    }
    if(y <= s) break;
  }
  return propval;
}

//----------------------------------------------------------------
// Function: rkg_one_zero
//
// Special-case sampler for z near zero, due to pGIG with a=0 becomes inverse gamma
//----------------------------------------------------------------
double rkg_one_zero() {
  double p = 0.347094406263; // this is TRUNC dependent!!
  double q = 0.741907333966;
  double q_over_p_plus_q = q / (p + q);
  double propval = 0.0;
  while(true) {
    if(R::runif(0.0, 1.0) < q_over_p_plus_q) {
      double k = 2.0 * MATH_PI2;
      double gamma_sample = R::rgamma(1.0, 1.0);
      propval = TRUNC + gamma_sample / k;
      if(std::isnan(propval) || propval <= 0.0)
        Rcpp::stop("Invalid propval in rkg_one_zero (right branch)");
    } else {
      propval = 1.0;
      // Gamma with shape=1.5 and scale=8.0
      while(propval > TRUNC) {
        double gamma_sample = R::rgamma(1.5, 8.0);
        propval = 1.0 / gamma_sample;
      }
      if(std::isnan(propval) || propval <= 0.0)
        Rcpp::stop("Invalid propval in rkg_one_zero (left branch)");
    }

    double s = a_coef(0, propval);
    if(std::isnan(s)) Rcpp::stop("s is NaN in rkg_one_zero");
    double y = R::runif(0.0, 1.0) * s;
    int n = 0;
    while(true) {
      n++;
      if(n % 2 == 1) {
        s -= a_coef(n, propval);
        if(y <= s) break;
      } else {
        s += a_coef(n, propval);
        if(y > s) break;
      }
    }
    if(y <= s) break;
  }
  return propval;
}

//----------------------------------------------------------------
// Function: rkg_b
//
// Returns the sum of 'b' independent samples from rkg_one.
//----------------------------------------------------------------
double rkg_b(int b, double z) {
  double retval = 0.0;
  for(int i = 0; i < b; i++) {
    retval += rkg_one(z);
  }
  return retval;
}

//----------------------------------------------------------------
// Function: rkg
//
// Vectorized version: takes two NumericVectors b and c (of equal length or 1)
// and returns a NumericVector of corresponding rkg_b samples.
//----------------------------------------------------------------
//' Sample Kolmogorov-Gamma random variables
//'
//' A random variable \eqn{X} follows Kolmogorov-Gamma(b,c) distribution, in short KG(b,c), if
//' \deqn{
//'  X \stackrel{d}{=} \dfrac{1}{2\pi^2}\sum_{k=1}^\infty \dfrac{\epsilon_k}{k^2 + c^2/(4\pi^2)}, \quad \epsilon_k\stackrel{iid}{\sim} Gamma(b,1)
//' }
//' where \eqn{\stackrel{d}{=}} denotes equality in distribution. 
//' The random variate generation is based on alternating series method, a fast and exact method (without infinite sum truncation) implemented in cpp. 
//' This function only supports integer b, which is sufficient for cobin and micobin regression models.  
//'
//' @param n The number of samples. 
//' @param b First parameter, positive integer (1,2,...). Length must be 1 or n.
//' @param c Second parameter, real, associated with tilting. Length must be 1 or n.
//' @return It returns n independent Kolmogorov-Gamma(\code{b[i]},\code{c[i]}) samples. If input b or c is scalar, it is assumed to be length n vector with same entries.  
//' @examples
//' \dontrun{
//' rkgcpp(1000, 1, 2)
//' rkgcpp(1000, 1, rnorm(1000))
//' rkgcpp(1000, rep(c(1,2),500), rnorm(1000))
//' }
//' @export
// [[Rcpp::export]]
NumericVector rkgcpp(int n, NumericVector b, NumericVector c) {
  // If b has length 1, replicate its value to length n.
  if (b.size() == 1) {
    NumericVector b_rep(n, b[0]);
    b = b_rep;
  }
  // If c has length 1, replicate its value to length n.
  if (c.size() == 1) {
    NumericVector c_rep(n, c[0]);
    c = c_rep;
  }
  if (b.size() != n)
    Rcpp::stop("rkgcpp error: Length of b must equal 1 or n.");
  if (c.size() != n)
    Rcpp::stop("rkgcpp error: Length of c must equal 1 or n.");

  // Pre-check b vector and convert to integer vector.
  std::vector<int> b_int(n);
  for (int i = 0; i < n; i++) {
    double b_val = b[i];
    if (b_val <= 0)
      Rcpp::stop("rkgcpp error: All elements of b must be positive.");
    if (std::fabs(b_val - std::round(b_val)) > 1e-10)
      Rcpp::stop("rkgcpp error: Non-integer value found in b.");
    b_int[i] = static_cast<int>(std::round(b_val));
  }

  // Now compute the output using the pre-validated b_int.
  NumericVector ret(n);
  for (int i = 0; i < n; i++) {
    ret[i] = rkg_b(b_int[i], c[i]);
  }
  return ret;
}

// //----------------------------------------------------------------
// // Function: rkg_vector
// //
// // Vectorized version: takes two NumericVectors b and z (of equal length)
// // and returns a NumericVector of corresponding rkg_b samples.
// //----------------------------------------------------------------
// // [[Rcpp::export]]
// NumericVector rkg_vector(NumericVector b, NumericVector z) {
//   int n = b.size();
//   if(z.size() != n) {
//     Rcpp::stop("Lengths of b and z must be equal");
//   }
//   NumericVector ret(n);
//   for(int i = 0; i < n; i++) {
//     ret[i] = rkg_b(static_cast<int>(b[i]), z[i]);
//   }
//   return ret;
// }

//----------------------------------------------------------------
// Function: rinvgauss
//
// Samples from the inverse Gaussian distribution.
//----------------------------------------------------------------
double rinvgauss(double mu, double lambda) {
  double u = R::rnorm(0.0, 1.0);
  double y = u * u;
  double x = mu + (mu * mu * y) / (2.0 * lambda) -
    (mu / (2.0 * lambda)) *
    std::sqrt(4.0 * mu * lambda * y + mu * mu * y * y);
  if(R::runif(0.0, 1.0) <= mu / (mu + x)) {
    return x;
  } else {
    return (mu * mu) / x;
  }
}

//----------------------------------------------------------------
// Function: rinvgauss_mulambda_export
//
// Returns a NumericVector of n inverse Gaussian samples.
// Exported for use in R.
//----------------------------------------------------------------
// [[Rcpp::export]]
NumericVector rinvgauss_mulambda_export(int n, double mu, double lambda) {
  NumericVector ret(n);
  for(int i = 0; i < n; i++) {
    ret[i] = rinvgauss(mu, lambda);
  }
  return ret;
}

//----------------------------------------------------------------
// Function: rgig_half
//
// Generates one sample from the generalized inverse Gaussian
// parameter p = -3/2.
//
// Based on Alg 2 of https://arxiv.org/pdf/2401.00749v2
//
//----------------------------------------------------------------
double rgig_negonehalf(double a, double b) {
  double sqrt_ab = std::sqrt(a * b);
  double w = sqrt_ab / (1.0 + sqrt_ab);
  double exp_contribution = R::rgamma(1.0, 2.0 / b); // Gamma(shape=1, scale=2/b)
  double x = 0.0;
  if(R::runif(0.0, 1.0) < w) {
    x = rinvgauss(std::sqrt(a / b), a);
  } else {
    x = 1.0 / rinvgauss(std::sqrt(b / a), b);
  }
  return 1.0 / (x + exp_contribution);
}

//----------------------------------------------------------------
// Function: rgig_negonehalf_export
//
// Returns a NumericVector of n samples from rgig_negonehalf(a, b).
// Exported for use in R.
//----------------------------------------------------------------
// [[Rcpp::export]]
NumericVector rgig_negonehalf_export(int n, double a, double b) {
  NumericVector ret(n);
  for(int i = 0; i < n; i++) {
    ret[i] = rgig_negonehalf(a, b);
  }
  return ret;
}
