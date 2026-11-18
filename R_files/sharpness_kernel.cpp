#include <Rcpp.h>
using namespace Rcpp;
// Exact multinomial bootstrap of standardized Gaussian bin sums.
// [[Rcpp::export]]
NumericMatrix gaussian_sharpness(int s, int B, int repetitions) {
  RNGScope scope;
  NumericMatrix result(repetitions, 3);
  std::vector<double> w(s);
  for (int r = 0; r < repetitions; ++r) {
    double mean = 0, sumsq = 0;
    for (int j = 0; j < s; ++j) {
      w[j] = R::rnorm(0, 1); mean += w[j]; sumsq += w[j] * w[j];
    }
    mean /= s;
    double ideal = std::sqrt(sumsq / s - mean * mean);
    double bootstrap_mean = 0, bootstrap_ss = 0;
    for (int b = 0; b < B; ++b) {
      double value = 0;
      for (int j = 0; j < s; ++j) value += w[(int)(R::unif_rand() * s)];
      value /= std::sqrt((double)s);
      double delta = value - bootstrap_mean;
      bootstrap_mean += delta / (b + 1);
      bootstrap_ss += delta * (value - bootstrap_mean);
    }
    double estimated = std::sqrt(bootstrap_ss / (B - 1));
    result(r, 0) = ideal - 1;
    result(r, 1) = estimated - ideal;
    result(r, 2) = estimated - 1;
    if (r % 100 == 0) checkUserInterrupt();
  }
  return result;
}
