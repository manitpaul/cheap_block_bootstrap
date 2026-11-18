# Final International Stroke Trial study

From this `R_files` directory, run:

```sh
Rscript real_data_ist_final.R --phase=observed --cores=2
Rscript real_data_ist_final.R --phase=synthetic --cores=2 --reps=300
Rscript real_data_ist_final.R --phase=timing --cores=1
Rscript real_data_ist_final.R --phase=plot --cores=1
```

Alternatively, call the script by its full path from any working directory.
Run timing only after other compute jobs have finished. `--phase=all` also runs
these phases in sequence. Statistical and timing phases require Rcpp and a
working C++ compiler. The
`--phase=plot` option reads saved results and the estimate in `design.csv`;
it needs only base R, without compiling C++ or accessing the patient data.
The configuration and required saved inputs are checked before regeneration.
For the statistical and timing phases, the
public corrected version-2 data default to
`../Real_data/raw/IST_corrected.csv`, relative to this folder, or a path supplied
as `--data=/absolute/path/IST_corrected.csv`. Missing data are downloaded from
Edinburgh DataShare. Provenance is recorded in the study's `data_source.csv` and
the project's `Real_data/data_sources.csv`. Individual records are
never exported by this script. Results go to `../results/final_ist/`, and figures
to `../PDF_images/`.

The target is the original IST prognosis-adjusted aspirin-minus-control risk
difference. The four strata of recorded `EXPDD` have cutpoints 0.40, 0.55, 0.70.
The estimate equals `n/(n0*n1) * sum_h(a_h - t_h*d_h/n_h)`. An empty stratum
contributes zero. No ridge fit, tuning parameter, or new risk model is used.

All bootstraps stratify by aspirin allocation. Original arm masses are
preserved. An m-out-of-n resample allocates m between arms proportionally (nearest
integer), restores original arm masses, and multiplies its centered error by
sqrt(m/n). Cheap bootstrap allocates s bins between arms proportionally, uses
balanced bins within each arm, resamples the two sets of bins independently, and
restores each arm's original mass. BLB uses K=ceiling(n/b) disjoint, balanced
subsets within each arm; a bag combines one subset from each arm. Each bag has B
independent multinomial draws with the full original arm sample sizes. Errors
are centered at the arm-standardized bag estimator. The bag standard deviations
and lower/upper error order statistics are averaged separately. B is the inner
draw count, so BLB performs K*B estimates. The standard error is the average of
bag standard deviations, not the square root of an average variance.

Small requested sizes are stress tests. In particular, log(n) bags are outside
the usual large-subset recommendations for BLB. Treatment stratification makes
all estimates defined and avoids silently discarding one-arm resamples. This
is a conditional-on-arm-sizes adaptation of the general construction; it is not
a direct application of an iid equal-bin theorem to the actual trial's
minimisation randomization scheme.

Accuracy runs use exact categorical sampling identities, since the statistic
only depends on 16 treatment/stratum/outcome cell counts. These preserve exactly
the distributions of row resampling, random partitions, and disjoint subsets.
Separate serial timings use generic row resampling and bag-weight computation,
including cell encoding and partition construction. The figures also report
the much faster specialized exact-cell ordinary bootstrap. Thus generic
speedups must not be interpreted as superiority to every optimized ordinary
implementation.

Synthetic data sample recorded prognostic probabilities with replacement,
generate fresh independent Bernoulli(1/2) allocations, and generate outcomes
with logit probability logit(q)+(A-1/2)*h. q is the recorded EXPDD value. Its
distribution is empirical; the resulting outcomes are synthetic. The known
population target is the empirical average of expit(logit(q)+h/2) minus
expit(logit(q)-h/2). Independent reference data sets estimate finite-n sampling
standard deviations. Calibration and power refer only to this synthetic model.

Seeds use separate L'Ecuyer streams per data set. `checkpoints/` allows identical
runs to resume. Configuration or input-data changes are rejected if an existing
run is found. Select a new directory with `--run=new_name` for a changed study;
also update the configuration version when changing the DGP or algorithm.
Configurations, references, summaries,
replicate-level results, timings, source checksum, and session information are
saved. Run `check_ist_final.R` after completion for numerical invariants.
