#!/usr/bin/env Rscript
# Reproduce the reported experiments, then measure all runtimes serially.
args <- commandArgs(TRUE)
opt <- function(key, default) {
  x <- args[startsWith(args, paste0('--', key, '='))]
  if (length(x)) sub(paste0('^--', key, '='), '', tail(x, 1)) else default
}
script_dir <- dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1])))
cores <- as.integer(opt('cores', '2'))
phase <- opt('phase', 'all')
stopifnot(is.finite(cores), cores >= 1L,
          phase %in% c('all', 'statistics', 'timing', 'figures'))
if (.Platform$OS.type != 'unix' && cores > 1L)
  stop('The main simulation uses forked workers: use macOS/Linux, or --cores=1.')
needed <- c('Rcpp', 'RcppArmadillo')
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('Install the following R packages: ', paste(missing, collapse=', '))
rscript <- file.path(R.home('bin'), 'Rscript')
run <- function(file, arguments=character()) {
  message('\nRunning ', file, ' ', paste(arguments, collapse=' '))
  status <- system2(rscript, c(shQuote(file.path(script_dir, file)), shQuote(arguments)))
  if (!identical(status, 0L)) stop(file, ' failed with exit status ', status)
}
main <- c('--reps=200', '--B=200,500', '--run=final_simulation')
tuning <- c('--cases=log_gamma,inverse_gamma', '--n-gamma=50000',
            '--B=500', '--reps=100', '--run=final_tuning')
workers <- paste0('--cores=', cores)
if (phase %in% c('all', 'statistics')) {
  run('final_simulation.R', c('--phase=statistics', main, workers))
  run('final_simulation.R', c('--phase=tuning', tuning, workers))
  run('real_data_ist_final.R', c('--phase=observed', workers))
  run('real_data_ist_final.R', c('--phase=synthetic', '--reps=300', workers))
  run('final_sharpness.R', as.character(cores))
}
if (phase %in% c('all', 'timing')) {
  message('Serial timing phase: stop other computational jobs before running this phase.')
  run('final_simulation.R', c('--phase=timing', main, '--cores=1', '--timing-reps=5'))
  run('final_simulation.R', c('--phase=tuning-timing', tuning, '--cores=1', '--timing-reps=5'))
  run('real_data_ist_final.R', c('--phase=timing', '--cores=1'))
}
if (phase %in% c('all', 'figures')) {
  run('check_final_simulation.R', 'final_simulation')
  run('check_final_tuning.R', 'final_tuning')
  run('plot_final_simulation.R', 'final_simulation')
  run('plot_final_tuning.R', 'final_tuning')
  run('real_data_ist_final.R', c('--phase=plot', '--cores=1'))
  run('check_ist_final.R')
  run('final_sharpness.R', '--plot-only')
}
message('\nRequested final-study phases completed.')
