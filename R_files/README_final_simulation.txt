FINAL SIMULATION STUDY

Requirements: R, Rcpp, RcppArmadillo, a C++ compiler. Only base R plotting and
parallel are otherwise needed. Drivers limit BLAS/OpenMP threads and parallelize
independent data sets, using deterministic L'Ecuyer-CMRG streams and checkpoints.
The macOS driver can link the Fortran runtime bundled with R when Makeconf names
a missing external toolchain; it does not install system software.

Main experiment (the reported run):
Rscript final_simulation.R --phase=statistics --reps=200 --cores=2 --B=200,500 --run=final_simulation

After all parallel studies have stopped, collect serial timings:
Rscript final_simulation.R --phase=timing --reps=200 --cores=1 --B=200,500 --timing-reps=5 --run=final_simulation
Rscript check_final_simulation.R final_simulation
Rscript plot_final_simulation.R final_simulation

Change --run when changing the statistical configuration. A saved configuration
cannot silently be mixed with a new one. Increase --reps for a more precise study.
Timing includes the full-data point fit, sufficient-statistic construction or
kernel preprocessing, all resampling, and the final SD/interval calculations.
It excludes data generation and compilation. One warm-up run is discarded and
the median of five subsequent serial runs is used. Timings evaluate the null
statistic once; the Monte Carlo experiment shares computations across alternatives.
For the U-statistic, runtime files additionally report the initial full-data fit
and cache construction, and the remaining resampling/aggregation time. The
remaining BLB time includes its per-bag kernel construction. Medians of separate
components need not add exactly to the median total. These timing fields do not
change the statistical experiment.

Models:
  Log mean: Gamma(shape=.005, rate=.005), multiplied by exp(theta).
  Capped reciprocal mean: Gamma(shape=.05, rate=.05), divided by theta;
    statistic min(1000, 1/mean(X)).
  Pairwise Euclidean distance: Z ~ .98 N_3(0,I)+.02 N_3(10 e1,I), X=a Z;
    unbiased order-two U-statistic of ||X_i-X_j||. Dimension three prevents the
    one-dimensional sorting shortcut. Exact U variance is computed using the
    Hoeffding decomposition and deterministic Gaussian quadrature (40- versus
    60-node rules agree to 1e-9). Both ordinary and CBB use exact cached kernels;
    CBB timing streams pairs directly into an s-by-s cache.
    BLB uses the weighted multinomial implementation, evaluating the exact
    quadratic form from its cached within-bag kernel.
  Regularized precision: X ~ N_20(0, diag(tau,1,...,1)); statistic
    e1' [ (1/n) sum_i (X_i-Xbar)(X_i-Xbar)' + .1 I ]^-1 e1.
    The population target is 1/(tau+.1). A Schur complement evaluates the five
    first-coordinate scale alternatives exactly, without repeated eigenanalysis.
    Sampling SD is estimated from 20,000 independent Wishart(n-1,I)/n draws;
    reference draws and their Monte Carlo standard errors are saved.

Sample sizes are 10,000/50,000 for Gamma, 1,000/4,000 for U-statistics, and
2,000/8,000 for precision. Rounded bootstrap budgets are B=200/500. The smaller
U-statistic sizes reflect the genuine quadratic exact ordinary-bootstrap cost.
Each experiment has 200 independent data sets and local effect indices d=0:4.
The maximum Bernoulli Monte Carlo standard error is 3.54 percentage points;
at probability .05 it is 1.54 points. Paired method differences are also saved.

Methods:
  Ordinary bootstrap; m/n bootstrap with m=ceil(log n), ceil(sqrt n);
  CBB with s approximately log n, sqrt n, and exactly B;
  BLB with b=ceil(log n), ceil(sqrt n), K=ceil(n/b) disjoint balanced bags.
CBB counts are rounded to a nearest divisor on the logarithmic scale so bins
are exactly equal. All n values are multiples of B. Actual settings are saved.

BLB draws Multinomial(n; 1/b_j,...,1/b_j) within each bag of actual size b_j.
Each bag has B inner draws, so total BLB fits equal K*B, not B. Within-bag
bootstrap errors center at the bag's plugin parameter; the reported SD averages
bag SDs, and interval error endpoints average bag quantiles. Uncentered bag fits
are never pooled. For a U-statistic the bag plugin is the V-statistic, not U_b.
The root-m/n method rescales centered errors by sqrt(m/n); CBB needs no rescaling.
All methods use unstudentized 95% basic intervals. At B=200 the root ranks are
5 and 196; at B=500 they are 12 and 488. High apparent power accompanied by
inflated null rejection is explicitly a calibration failure, not a power benefit.

BLB's b=log n is a deliberately aggressive small-support stress regime. Theorem 3
of Kleiner et al. additionally needs b of order at least sqrt n, regular expansion
conditions, and excludes finite-inner-bootstrap Monte Carlo error. K=n/b alone
does not imply its higher-order correctness conclusion.

Tuning driver (gamma examples only; same conventions):
Rscript final_simulation.R --phase=tuning --cases=log_gamma,inverse_gamma --n-gamma=50000 --B=500 --reps=100 --cores=2 --run=final_tuning
Rscript final_simulation.R --phase=tuning-timing --cases=log_gamma,inverse_gamma --n-gamma=50000 --B=500 --reps=100 --cores=1 --timing-reps=5 --run=final_tuning
Rscript check_final_tuning.R final_tuning
Rscript plot_final_tuning.R final_tuning

Source/reference: Kleiner et al. (2014), A scalable bootstrap for massive data,
JRSS B 76:795-816, DOI 10.1111/rssb.12050; Algorithm 1 and Theorem 3.
No failed or unfavorable method runs are removed. The validation script checks
all records, actual count budgets, equal CBB bins, and paired comparisons.
