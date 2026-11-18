#!/usr/bin/env Rscript
# Integrity checks on every tuning replicate; no favorable-result selection.
args<-commandArgs(TRUE)
script_dir<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
run<-if(length(args))args[1]else'final_tuning';out<-file.path(script_dir,'results',run)
x<-read.csv(file.path(out,'tuning_replicates.csv'));cfg<-readRDS(file.path(out,'config.rds'))
stopifnot(all(is.finite(x$relative_error)),all(!is.na(x$rejection)),all(!is.na(x$coverage)))
stopifnot(!anyDuplicated(x[c('case','rep','delta','method','size')]))
counts<-aggregate(rep~case+delta+method+size,x,length)
stopifnot(all(counts$rep==cfg$reps),setequal(unique(x$case),cfg$cases),setequal(unique(x$delta),cfg$deltas))
stopifnot(setequal(unique(x$method),c('ordinary','mn','cbb','blb')))
cb<-x[x$method=='cbb',];stopifnot(all(cb$n%%cb$size==0))
settings<-unique(x[c('case','n','B','method','size','requested')])
settings$K<-ifelse(settings$method=='blb',ceiling(settings$n/settings$size),1)
settings$total_draws<-settings$K*settings$B
settings$minimum_bag_size<-ifelse(settings$method=='blb',floor(settings$n/settings$K),NA_integer_)
settings$maximum_bag_size<-ifelse(settings$method=='blb',ceiling(settings$n/settings$K),NA_integer_)
stopifnot(nrow(x)==nrow(settings)*cfg$reps*length(cfg$deltas))
write.csv(settings,file.path(out,'tuning_method_settings.csv'),row.names=FALSE)
o<-x[x$method=='ordinary',c('case','n','B','rep','delta','rejection','relative_error')]
q<-merge(x,o,by=c('case','n','B','rep','delta'),suffixes=c('','_ordinary'))
q$rejection_diff<-as.numeric(q$rejection)-as.numeric(q$rejection_ordinary)
q$mse_diff<-q$relative_error^2-q$relative_error_ordinary^2
p<-do.call(rbind,lapply(split(q,interaction(q$case,q$delta,q$method,q$size,drop=TRUE)),function(z)data.frame(z[1,c('case','n','B','delta','method','size')],reps=nrow(z),rejection_difference=mean(z$rejection_diff),rejection_difference_mcse=sd(z$rejection_diff)/sqrt(nrow(z)),mse_difference=mean(z$mse_diff),mse_difference_mcse=sd(z$mse_diff)/sqrt(nrow(z)))))
write.csv(p,file.path(out,'tuning_paired_comparisons.csv'),row.names=FALSE)
writeLines(c(paste('Validated rows:',nrow(x)),paste('Replicates per tuning setting:',cfg$reps),'All settings and replicates retained; equal CBB bins; actual BLB budgets recorded.'),file.path(out,'validation.txt'))
cat('Validated',nrow(x),'tuning rows\n')

# Confirm complete timing coverage before joining statistics to runtime figures.
if(file.exists(file.path(out,'tuning_runtime_summary.csv'))){
 tt<-read.csv(file.path(out,'tuning_runtime_serial.csv'));ts<-read.csv(file.path(out,'tuning_runtime_summary.csv'));timing_keys<-c('case','n','B','method','size')
 stopifnot(!anyDuplicated(ts[timing_keys]),!anyDuplicated(tt[c(timing_keys,'repeat_id')]),all(is.finite(tt$seconds)),all(tt$seconds>0),all(ts$q25>0),all(ts$q25<=ts$seconds),all(ts$seconds<=ts$q75))
 actual<-unique(x[timing_keys]);stopifnot(nrow(ts)==nrow(actual),nrow(merge(actual,ts,by=timing_keys))==nrow(actual))
 cat('Validated complete serial timing coverage:',nrow(ts),'settings\n')
}
