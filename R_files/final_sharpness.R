#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
plot_only <- '--plot-only' %in% args
worker_args <- setdiff(args, '--plot-only')
cores <- if (length(worker_args)) as.integer(worker_args[1]) else 2L
script_dir <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1])))
out <- file.path(script_dir, 'results', 'final_sharpness')
fig <- file.path(dirname(script_dir), 'PDF_images')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
dir.create(fig, showWarnings=FALSE)
if (!plot_only) {
Rcpp::sourceCpp(file.path(script_dir, 'sharpness_kernel.cpp'), cacheDir=file.path(script_dir, 'cache'))
grid <- expand.grid(s=c(50L, 200L, 1000L), B=c(50L, 200L, 1000L))
RNGkind("L'Ecuyer-CMRG"); set.seed(20261122)
streams <- vector('list', nrow(grid))
for (i in seq_len(nrow(grid))) {
  streams[[i]] <- .Random.seed
  .Random.seed <- parallel::nextRNGStream(.Random.seed)
}
worker <- function(i) {
  assign('.Random.seed', streams[[i]], envir=.GlobalEnv)
  e <- gaussian_sharpness(grid$s[i], grid$B[i], 1000L)
  write.csv(e, file.path(out, sprintf('errors_s%d_B%d.csv', grid$s[i], grid$B[i])), row.names=FALSE)
  data.frame(grid[i, ], repetitions=1000L, compression=sqrt(mean(e[,1]^2)),
    monte_carlo=sqrt(mean(e[,2]^2)), combined=sqrt(mean(e[,3]^2)))
}
rows <- parallel::mclapply(seq_len(nrow(grid)), worker, mc.cores=cores, mc.set.seed=FALSE)
stopifnot(!any(vapply(rows, inherits, FALSE, 'try-error')))
z <- do.call(rbind, rows)
write.csv(z, file.path(out, 'summary.csv'), row.names=FALSE)
} else z <- read.csv(file.path(out, 'summary.csv'))
pdf(file.path(fig, 'cbb_sharpness.pdf'), width=9, height=3.7, pointsize=12, family='Helvetica', useDingbats=FALSE)
par(mfrow=c(1,3), mar=c(4.1,4.5,2.4,.8), mgp=c(2.8,.8,0), las=1, cex=.95)
cols <- c('#0072B2', '#D55E00', '#117733')
for (h in 1:3) {
  value <- c('compression','monte_carlo','combined')[h]
  horizontal <- if (h==1) 's' else 'B'
  plot(NA, xlim=c(45,1100), ylim=range(z[[value]])*c(.85,1.2), log='xy',
    xlab=if(h==1) 'Number of bins s' else 'Bootstrap draws B', ylab='Relative RMSE',
    main=c('Compression error','Monte Carlo error','Combined error')[h], xaxt='n')
  axis(1, at=c(50,200,1000), labels=c('50','200','1,000'))
  for (j in 1:3) {
    v <- c(50,200,1000)[j]
    a <- if(h==1) z[z$B==v,] else z[z$s==v,]
    a <- a[order(a[[horizontal]]),]
    lines(a[[horizontal]],a[[value]],type='b',pch=15+j,col=cols[j],lwd=1.3)
    if(h==3) lines(a$B,sqrt((1/a$s+1/a$B)/2),col=cols[j],lty=2)
  }
  if(h<3) lines(c(50,200,1000),1/sqrt(2*c(50,200,1000)),lty=2)
  legend('bottomleft',paste(if(h==1)'B =' else 's =',c(50,200,1000)),col=cols,
    pch=16:18,lty=1,bty='n',cex=.95)
}
dev.off()
if (!plot_only) writeLines(capture.output(sessionInfo()), file.path(out,'sessionInfo.txt'))
cat('Exact Gaussian sharpness experiment and figure completed.\n')
