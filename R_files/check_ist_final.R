#!/usr/bin/env Rscript
f<-grep('^--file=',commandArgs(),value=TRUE);script_dir<-dirname(normalizePath(sub('^--file=','',f[1])));root<-dirname(script_dir);out<-file.path(root,'results','final_ist')
cfg<-readRDS(file.path(out,'config.rds'));raw<-read.csv(file.path(out,'synthetic_replicates.csv'));obs<-read.csv(file.path(out,'observed_runs.csv'));tune<-read.csv(file.path(out,'observed_tuning.csv'));tm<-read.csv(file.path(out,'runtime_serial.csv'));tt<-read.csv(file.path(out,'tuning_runtime_serial.csv'))
stopifnot(nrow(raw)==cfg$R*length(cfg$B)*8L*3L,nrow(obs)==cfg$observed_runs*8L,nrow(tune)==cfg$tuning_runs*length(cfg$tuning)*3L)
for(x in list(raw,obs,tune)){stopifnot(all(is.finite(x$se)),all(x$se>0),all(is.finite(x$lower)),all(is.finite(x$upper)),all(x$lower<=x$upper),all(x$inner_fits==x$K*x$B))}
stopifnot(nrow(tm)==7L*9L,nrow(tt)==5L*3L*length(cfg$tuning),all(is.finite(tm$seconds)),all(tm$seconds>0),all(is.finite(tt$seconds)),all(tt$seconds>0),all(raw$rejection==(raw$lower>0|raw$upper<0)),all(raw$coverage==(raw$lower<=raw$target&raw$upper>=raw$target)))
stopifnot(all(table(raw$method,raw$B,raw$delta)==cfg$R))
design<-read.csv(file.path(out,'design.csv'));stopifnot(design$n==19285L,abs(design$estimate+0.01367058)<1e-8)
tr<-read.csv(file.path(out,'sampling_truth.csv'));stopifnot(tr$target[1]==0,all(diff(tr$target)<0),all(tr$sigma>0))
stopifnot(length(readRDS(file.path(out,'observed_reference.rds')))==cfg$reference_B,nrow(readRDS(file.path(out,'sampling_reference.rds')))==cfg$reference_R)
figures<-file.path(root,'PDF_images',c('ist_observed_final.pdf','ist_tuning_final.pdf','ist_synthetic_final.pdf'));stopifnot(all(file.exists(figures)),all(file.info(figures)$size>1000))
cat(sprintf('PASS: %d synthetic records, %d observed records, %d tuning records; all fits valid and finite.\n',nrow(raw),nrow(obs),nrow(tune)))
