#!/usr/bin/env Rscript
args<-commandArgs(TRUE);script_dir<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
run<-if(length(args))args[1]else'final_simulation';out<-file.path(script_dir,'results',run)
x<-read.csv(file.path(out,'replicates.csv'));cfg<-readRDS(file.path(out,'config.rds'))
ns<-function(ca)if(ca%in%c('log_gamma','inverse_gamma'))cfg$n_gamma else if(ca=='u_distance')cfg$n_u else cfg$n_precision
expected_datasets<-cfg$reps*sum(vapply(cfg$cases,function(ca)length(ns(ca)),0L))
stopifnot(nrow(x)==expected_datasets*length(cfg$B)*8*length(cfg$deltas))
stopifnot(all(is.finite(x$estimate)),all(is.finite(x$se)),all(x$se>=0),all(is.finite(x$truth)),all(x$truth>0),all(is.finite(x$relative_error)),all(x$upper>=x$lower))
counts<-aggregate(rep~case+n+B+delta+method,x,length);stopifnot(all(counts$rep==cfg$reps))
stopifnot(!anyDuplicated(x[c('case','n','B','rep','delta','method')]))
bb<-x[startsWith(x$method,'blb'),];stopifnot(all(bb$K==ceiling(bb$n/bb$size)),all(bb$total_draws==bb$K*bb$B))
cb<-x[startsWith(x$method,'cbb'),];stopifnot(all(cb$n%%cb$size==0),all(cb$size[cb$method=='cbb_B']==cb$B[cb$method=='cbb_B']))
ordinary<-x[x$method=='ordinary',c('case','n','B','rep','delta','rejection','relative_error')]
q<-merge(x,ordinary,by=c('case','n','B','rep','delta'),suffixes=c('','_ordinary'));q$rejection_diff<-as.numeric(q$rejection)-as.numeric(q$rejection_ordinary);q$squared_error_diff<-q$relative_error^2-q$relative_error_ordinary^2
p<-do.call(rbind,lapply(split(q,interaction(q$case,q$n,q$B,q$delta,q$method,drop=TRUE)),function(z)data.frame(z[1,c('case','n','B','delta','method')],reps=nrow(z),rejection_difference=mean(z$rejection_diff),rejection_difference_mcse=sd(z$rejection_diff)/sqrt(nrow(z)),mse_difference=mean(z$squared_error_diff),mse_difference_mcse=sd(z$squared_error_diff)/sqrt(nrow(z)))))
write.csv(p,file.path(out,'paired_comparisons.csv'),row.names=FALSE)
settings<-unique(x[c('case','n','B','method','size','K','total_draws')]);settings$minimum_bag_size<-ifelse(startsWith(settings$method,'blb'),floor(settings$n/settings$K),NA_integer_);settings$maximum_bag_size<-ifelse(startsWith(settings$method,'blb'),ceiling(settings$n/settings$K),NA_integer_)
write.csv(settings,file.path(out,'method_settings.csv'),row.names=FALSE)
writeLines(c(paste('Validated rows:',nrow(x)),paste('Datasets per model/sample-size:',cfg$reps),'All estimates and standard errors finite; all methods retained.','Paired comparisons use common data, but independent bootstrap random numbers by method.','BLB computes B draws in every one of K near-equal disjoint bags.'),file.path(out,'validation.txt'))
cat('Validated',nrow(x),'rows\n')

# Confirm complete timing coverage before joining statistics to runtime figures.
if(file.exists(file.path(out,'runtime_serial_summary.csv'))){
 tt<-read.csv(file.path(out,'runtime_serial.csv'));ts<-read.csv(file.path(out,'runtime_serial_summary.csv'));timing_keys<-c('case','n','B','method')
 stopifnot(!anyDuplicated(ts[timing_keys]),!anyDuplicated(tt[c(timing_keys,'repeat_id')]),all(is.finite(tt$seconds)),all(tt$seconds>0),all(ts$q25>0),all(ts$q25<=ts$seconds),all(ts$seconds<=ts$q75))
 actual<-unique(x[timing_keys]);stopifnot(nrow(ts)==nrow(actual),nrow(merge(actual,ts,by=timing_keys))==nrow(actual))
 cat('Validated complete serial timing coverage:',nrow(ts),'settings\n')
}
