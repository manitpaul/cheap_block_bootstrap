# Cheap Block Bootstrap

Reproducibility code for **Cheap Block Bootstrap**, by Manit Paul, Abir Sarkar,
Raaz Dwivedi, and Lester Mackey.

The cheap block bootstrap (CBB) partitions observations into bins, stores
reusable summaries, and applies shared bootstrap weights within each bin. This
repository compares CBB with the ordinary bootstrap, the *m*-out-of-*n*
bootstrap, and the bag of little bootstraps (BLB). It contains the analysis
code, saved simulation outputs and aggregate clinical results, and the figures
used in the manuscript. Manuscript drafts and participant-level data are not
included.

## Requirements

- R and a C++ toolchain compatible with R.
- R packages `Rcpp` and `RcppArmadillo`:

  ```r
  install.packages(c("Rcpp", "RcppArmadillo"))
  ```

The remaining dependencies are distributed with R. The parallel simulation
driver uses forked workers on macOS/Linux. Use `--cores=1` on Windows. Session
information from the reported runs is saved with each experiment.

## Reproduce figures from the supplied results

From the repository root, run:

```sh
Rscript R_files/run_final_studies.R --phase=figures
```

This validates the saved results and regenerates all eleven figures in
`PDF_images/`. It does not require participant-level data or rerunning the
Monte Carlo experiments. The driver checks that the Rcpp packages are installed,
although figure generation itself uses base R.

## Rerun the experiments

```sh
# Statistical experiments, with two worker processes:
Rscript R_files/run_final_studies.R --phase=statistics --cores=2

# Measure runtimes after other computational jobs have stopped:
Rscript R_files/run_final_studies.R --phase=timing --cores=1

# Validate results and regenerate figures:
Rscript R_files/run_final_studies.R --phase=figures
```

Alternatively, `Rscript R_files/run_final_studies.R --cores=2` runs all phases
in order. Full reruns are computationally substantial. Random seeds and
configurations are recorded, and deterministic checkpoints allow interrupted
statistical runs to resume. Checkpoints are not distributed here. Use a new
`--run` directory when changing an experiment's configuration; the individual
drivers reject incompatible saved configurations.

Wall-clock timings depend on hardware, operating system, numerical libraries,
and system load. The reported timings include estimator evaluation,
preprocessing, resampling, and interval/standard-error calculations, and exclude
data generation and compilation. Main simulation timings use a warm-up followed
by five serial measurements, summarized by their median.

## Experiments and conventions

The main simulations use 200 independent data sets and bootstrap budgets
`B = 200, 500`:

| Model | Sample sizes | Statistic |
| --- | --- | --- |
| Gamma, shape = rate = 0.005 | 10,000; 50,000 | Logarithm of the sample mean |
| Gamma, shape = rate = 0.05 | 10,000; 50,000 | Capped reciprocal sample mean |
| Three-dimensional Gaussian mixture | 1,000; 4,000 | Mean pairwise Euclidean distance, an order-two U-statistic |
| Twenty-dimensional Gaussian model | 2,000; 8,000 | Regularized precision-matrix entry |

The methods are ordinary bootstrap; *m*-out-of-*n* bootstrap with logarithmic
and square-root resample sizes; CBB with logarithmic, square-root, and `B` bins;
and BLB with logarithmic and square-root subset sizes. Actual integer settings
are saved with the results. Tuning experiments vary the compression/subsample
parameter in the two Gamma models. A separate Gaussian experiment examines the
two terms in the standard-error bound.

BLB uses `K = ceiling(n/b)` balanced disjoint bags and **B inner draws in each
bag**, for `K*B` weighted fits. Bag-level standard deviations and error quantiles
are averaged separately. The logarithmic subset size is an aggressive stress
test; `K = n/b` alone does not establish BLB's higher-order accuracy conditions.
All methods use centered errors and unstudentized basic intervals. The
*m*-out-of-*n* errors are scaled by `sqrt(m/n)`; CBB needs no such scaling.
Uncalibrated increases in rejection are not interpreted as power improvements.

For the U-statistic, weighted bootstrap fits use exact kernel calculations,
including cached within-bag kernels for BLB. Both preprocessing and resampling
are included in the comparisons. No failed or unfavorable methods are removed.

Detailed configurations and commands are in
[the simulation notes](R_files/README_final_simulation.txt).

## International Stroke Trial

The clinical analysis estimates the prognosis-adjusted aspirin-minus-control
risk difference for death or dependency at six months, using the original
observed-minus-expected adjustment across four fixed prognosis strata. The
complete-case analysis includes 19,285 participants.

The corrected version-2 data are available from
[Edinburgh DataShare, DOI 10.7488/ds/104](https://doi.org/10.7488/ds/104).
For statistical or timing reruns, the clinical driver downloads the source if
it is absent and stores it in the ignored directory
`Real_data/raw/IST_corrected.csv`. An existing file can be supplied to the
individual driver with `--data=/absolute/path/IST_corrected.csv`. The source URL
and checksum are recorded in `results/final_ist/data_source.csv`. Consult the
data provider for the data's terms of use.

The public repository contains aggregate counts, bootstrap summaries, and
Monte Carlo outputs, not individual clinical records. The synthetic experiment
generates new allocations and outcomes using the empirical distribution of
the recorded prognostic score. Its simulated outcomes are not observed trial
outcomes. Generic row-resampling timings are accompanied by an optimized
exact-cell ordinary-bootstrap benchmark. See
[the clinical analysis notes](R_files/README_ist_final.md) for the estimand,
resampling conventions, limitations, and phase-specific commands.

## Repository layout

| Location | Contents |
| --- | --- |
| `R_files/` | R drivers, plotting/validation code, C++ kernels, and detailed notes |
| `R_files/results/final_simulation/` | Main simulation configurations, records, references, summaries, and timings |
| `R_files/results/final_tuning/` | Tuning experiment outputs |
| `R_files/results/final_sharpness/` | Gaussian sharpness experiment outputs |
| `results/final_ist/` | Aggregate clinical and synthetic-experiment results |
| `PDF_images/` | Generated manuscript figures |

## Contact

Manit Paul: [paulman@wharton.upenn.edu](mailto:paulman@wharton.upenn.edu).
Please use GitHub issues for reproducibility questions or bug reports.

AI assistance was used in developing and reviewing portions of the code and
documentation. Statistical comparisons should be interpreted with the
calibration checks, Monte Carlo uncertainty, and implementation qualifications
documented above.
