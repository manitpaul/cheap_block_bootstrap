#!/usr/bin/env Rscript
# Final IST study. Each bootstrap conditions on its data set's treatment-arm sizes.
# The original trial's prognosis-adjusted risk-difference score replaces ridge.
args <- commandArgs(trailingOnly=TRUE)
opt <- function(k,d){z<-args[startsWith(args,paste0('--',k,'='))];if(length(z))sub(paste0('^--',k,'='),'',tail(z,1))else d}
filearg <- grep('^--file=',commandArgs(),value=TRUE)
script_dir <- dirname(normalizePath(sub('^--file=','',filearg[1])))
root <- dirname(script_dir); project <- root
out <- file.path(root,'results',opt('run','final_ist'));dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(root,'PDF_images'),showWarnings=FALSE)
cores <- as.integer(opt('cores','2'));R <- as.integer(opt('reps','300'));phase <- opt('phase','all')
stopifnot(is.finite(cores),cores>=1L,is.finite(R),R>=1L,phase%in%c('all','observed','synthetic','timing','plot'))
cfg <- list(version='IST-final-score-v1',R=R,B=c(200L,1000L),observed_B=1000L,
 observed_runs=30L,tuning_runs=15L,reference_B=50000L,reference_R=50000L,seed=2026111801L,
 tuning=c(10L,20L,50L,100L,200L,500L,1000L,2000L,5000L),arm_stratified=TRUE,
 blb='K=ceil(n/b) disjoint balanced bags in each arm; B inner draws per bag')
config_path<-file.path(out,'config.rds')
if(phase=='plot'&&!file.exists(config_path))stop('Saved configuration is missing: ',config_path)
if(file.exists(config_path)&&!identical(readRDS(config_path),cfg))stop('Stored configuration differs. Use --run=new_name and a new version for changed experiments; existing checkpoints are not reused.')
if(phase!='plot')saveRDS(cfg,config_path)
streams <- function(N,seed){RNGkind("L'Ecuyer-CMRG");set.seed(seed);x<-vector('list',N);for(i in seq_len(N)){x[[i]]<-.Random.seed;.Random.seed<-parallel::nextRNGStream(.Random.seed)};x}
setstream <- function(s)assign('.Random.seed',s,envir=.GlobalEnv)
parallel_map <- function(ids,fun){
 if(cores<=1)return(lapply(ids,fun))
 if(.Platform$OS.type=='unix')a<-parallel::mclapply(ids,fun,mc.cores=cores,mc.preschedule=TRUE,mc.set.seed=FALSE)
 else{cl<-parallel::makePSOCKcluster(cores);on.exit(parallel::stopCluster(cl));parallel::clusterExport(cl,ls(.GlobalEnv),envir=.GlobalEnv);parallel::clusterEvalQ(cl,Rcpp::sourceCpp(file.path(script_dir,'ist_kernel_final.cpp'),cacheDir=file.path(root,'tmp','rcpp')));a<-parallel::parLapply(cl,ids,fun)}
 if(any(vapply(a,inherits,logical(1),'try-error')))stop('Worker failed; no fits are discarded.');a
}
labels <- c(ordinary='Ordinary',mn_log='m/n: m=log n',mn_sqrt='m/n: m=sqrt(n)',
 cbb_log='CBB: s=log n',cbb_sqrt='CBB: s=sqrt(n)',cbb_B='CBB: s=B',
 blb_log='BLB: b=log n',blb_sqrt='BLB: b=sqrt(n)')
cols <- c(ordinary='#222222',mn_log='#AA4499',mn_sqrt='#CC79A7',cbb_log='#117733',
 cbb_sqrt='#44AA99',cbb_B='#0072B2',blb_log='#997700',blb_sqrt='#D55E00')
pchs <- c(ordinary=16,mn_log=17,mn_sqrt=15,cbb_log=1,cbb_sqrt=2,cbb_B=18,blb_log=0,blb_sqrt=5)
ltys <- c(ordinary=1,mn_log=2,mn_sqrt=3,cbb_log=2,cbb_sqrt=3,cbb_B=1,blb_log=2,blb_sqrt=3)
if(phase=='plot'){
 # Regeneration uses saved computational outputs only: no patient data or C++.
 required<-c('design.csv','observed_runs.csv','observed_tuning.csv','runtime_serial.csv',
             'tuning_runtime_serial.csv','synthetic_replicates.csv','observed_reference.rds')
 missing<-required[!file.exists(file.path(out,required))]
 if(length(missing))stop('Saved plot inputs are missing from ',out,': ',paste(missing,collapse=', '))
 design<-read.csv(file.path(out,'design.csv'));theta<-design$estimate
 if(length(theta)!=1L||!is.finite(theta))stop('Saved design.csv must contain one finite estimate.')
 source(file.path(script_dir,'plot_ist_final.R'))
 quit(save='no',status=0L)
}
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
if(!requireNamespace('Rcpp',quietly=TRUE))stop('Install Rcpp to run this study.')
Rcpp::sourceCpp(file.path(script_dir,'ist_kernel_final.cpp'),cacheDir=file.path(root,'tmp','rcpp'))
RNGkind("L'Ecuyer-CMRG")
grid <- function(n,B)data.frame(method=names(labels),kind=c('ordinary','mn','mn','cbb','cbb','cbb','blb','blb'),size=as.integer(c(n,ceiling(log(n)),ceiling(sqrt(n)),ceiling(log(n)),ceiling(sqrt(n)),B,ceiling(log(n)),ceiling(sqrt(n)))))
counts_from <- function(a,h,y)matrix(tabulate(1L+a*8L+h*2L+y,16L),nrow=2,byrow=TRUE)
data_path<-opt('data',file.path(project,'Real_data','raw','IST_corrected.csv'))
data_url<-'https://datashare.ed.ac.uk/server/api/core/bitstreams/eca52b56-6d1b-47b3-814a-0b7089f8c870/content'
if(!file.exists(data_path)){dir.create(dirname(data_path),recursive=TRUE,showWarnings=FALSE);download.file(data_url,data_path,mode='wb',method='libcurl')}
z <- read.csv(data_path)
source_path<-file.path(out,'data_source.csv');source_md5<-unname(tools::md5sum(data_path))
if(file.exists(source_path)&&!identical(as.character(read.csv(source_path)$md5),source_md5))stop('Input data differ from the stored study. Use --run=new_name; existing checkpoints are not reused.')
write.csv(data.frame(dataset='International Stroke Trial corrected version 2',doi='10.7488/ds/104',url=data_url,md5=source_md5),source_path,row.names=FALSE)
missing <- data.frame(arm=c('Aspirin','Control'),excluded=c(sum(z$RXASP=='Y'&!z$OCCODE%in%1:4),sum(z$RXASP!='Y'&!z$OCCODE%in%1:4)))
write.csv(missing,file.path(out,'missingness.csv'),row.names=FALSE)
z <- z[z$OCCODE%in%1:4,];a <- as.integer(z$RXASP=='Y');y<-as.integer(z$OCCODE%in%c(1,2));q<-z$EXPDD
stopifnot(all(is.finite(q)),all(q>0&q<1));h<-findInterval(q,c(.4,.55,.7));n<-length(y)
ct<-counts_from(a,h,y);theta<-ist_estimate(ct)
write.csv(data.frame(n=n,aspirin=sum(a),control=sum(1-a),adverse=sum(y),estimate=theta,unadjusted=mean(y[a==1])-mean(y[a==0])),file.path(out,'design.csv'),row.names=FALSE)
write.csv(aggregate(cbind(n=rep(1,n),aspirin=a,adverse=y,aspirin_adverse=a*y),list(stratum=h+1),sum),file.path(out,'strata.csv'),row.names=FALSE)
write.csv(grid(n,1000),file.path(out,'methods.csv'),row.names=FALSE)
meta<-data.frame(method=names(labels),size=grid(n,1000)$size)
meta$K<-ifelse(grepl('blb',meta$method),ceiling(n/meta$size),1)
meta$min_bag_arm0<-ifelse(grepl('blb',meta$method),floor(sum(1-a)/meta$K),NA)
meta$max_bag_arm0<-ifelse(grepl('blb',meta$method),ceiling(sum(1-a)/meta$K),NA)
meta$min_bag_arm1<-ifelse(grepl('blb',meta$method),floor(sum(a)/meta$K),NA)
meta$max_bag_arm1<-ifelse(grepl('blb',meta$method),ceiling(sum(a)/meta$K),NA)
write.csv(meta,file.path(out,'bag_sizes.csv'),row.names=FALSE)
summarize_one <- function(v,th,B,method,kind,size,run){data.frame(run=run,B=B,method=method,kind=kind,size=size,se=v['se'],lower=th-v['qhi'],upper=th-v['qlo'],K=v['K'],inner_fits=v['inner_fits'],row.names=NULL)}
if(phase%in%c('all','observed')){
 set.seed(cfg$seed);ref<-ist_reference(ct,cfg$reference_B);saveRDS(ref,file.path(out,'observed_reference.rds'))
 ss<-streams(cfg$observed_runs,cfg$seed+10L);mg<-grid(n,cfg$observed_B)
 worker<-function(r){setstream(ss[[r]]);do.call(rbind,lapply(seq_len(nrow(mg)),function(j){m<-mg[j,];v<-ist_resample(ct,cfg$observed_B,m$kind,m$size);summarize_one(v,theta,cfg$observed_B,m$method,m$kind,m$size,r)}))}
 ob<-do.call(rbind,parallel_map(seq_len(cfg$observed_runs),worker));ob$reference_se<-sd(ref);ob$relative_error<-ob$se/sd(ref)-1
 write.csv(ob,file.path(out,'observed_runs.csv'),row.names=FALSE)
 ss<-streams(cfg$tuning_runs,cfg$seed+20L)
 worker<-function(r){setstream(ss[[r]]);ans<-list();k<-0
 for(kind in c('mn','cbb','blb'))for(size in cfg$tuning){k<-k+1;v<-ist_resample(ct,1000,kind,size);ans[[k]]<-summarize_one(v,theta,1000,kind,kind,size,r)};do.call(rbind,ans)}
 tu<-do.call(rbind,parallel_map(seq_len(cfg$tuning_runs),worker));tu$reference_se<-sd(ref);tu$relative_error<-tu$se/sd(ref)-1
 referr<-sort(ref-theta);lo<-floor((length(ref)+1)*.025+1e-10);hi<-lo+ceiling((length(ref)+1)*.95-1e-10);refci<-theta-referr[c(hi,lo)]
 tu$endpoint_error<-sqrt(((tu$lower-refci[1])^2+(tu$upper-refci[2])^2)/2)/sd(ref)
 write.csv(tu,file.path(out,'observed_tuning.csv'),row.names=FALSE)
 cat('Observed analysis complete\n');flush.console()
}
if(phase%in%c('all','synthetic')){
 # q is the recorded prognostic probability; no new outcome model is fitted.
 qh<-tapply(q,h,mean);wh<-tabulate(h+1L,4L)/n;tau0<-sqrt(4*sum(wh*qh*(1-qh)))
 delta<-c(0,2,4);hs<--delta*tau0/(sqrt(n)*mean(q*(1-q)))
 truth<-vapply(hs,function(v)mean(plogis(qlogis(q)+v/2)-plogis(qlogis(q)-v/2)),0)
 sample_counts<-function(){ii<-sample.int(n,n,replace=TRUE);aa<-rbinom(n,1,.5);u<-runif(n);lapply(hs,function(v)counts_from(aa,h[ii],as.integer(u<plogis(qlogis(q[ii])+(aa-.5)*v))))}
 # Exact categorical sampling is equivalent to generating all n individual rows.
 # This makes a high-precision, independent finite-n reference inexpensive.
 set.seed(cfg$seed+30L)
 ref<-vapply(hs,function(v){p<-matrix(0,2,8);for(arm in 0:1){prob<-plogis(qlogis(q)+(arm-.5)*v);for(hh in 0:3){p[arm+1,2*hh+2]<-mean(prob*(h==hh));p[arm+1,2*hh+1]<-mean((1-prob)*(h==hh))}};ist_population_reference(p,n,cfg$reference_R)},numeric(cfg$reference_R));sig<-apply(ref,2,sd)
 saveRDS(ref,file.path(out,'sampling_reference.rds'))
 mc<-vapply(seq_len(ncol(ref)),function(j){v<-var(ref[,j]);sqrt(max(0,mean((ref[,j]-mean(ref[,j]))^4)-v^2)/(4*nrow(ref)*v))},0)
 write.csv(data.frame(delta=delta,logit_shift=hs,target=truth,sigma=sig,sigma_mcse=mc,tau0=tau0),file.path(out,'sampling_truth.csv'),row.names=FALSE)
 ss<-streams(R,cfg$seed+40L);cp<-file.path(out,'checkpoints');dir.create(cp,showWarnings=FALSE)
 worker<-function(r){fn<-file.path(cp,sprintf('task_%04d.rds',r));if(file.exists(fn))return(fn);setstream(ss[[r]]);cts<-sample_counts();ans<-list();k<-0L
 for(B in cfg$B){mg<-grid(n,B);for(j in seq_len(nrow(mg))){m<-mg[j,];for(e in seq_along(delta)){k<-k+1L;th<-ist_estimate(cts[[e]]);v<-ist_resample(cts[[e]],B,m$kind,m$size);one<-summarize_one(v,th,B,m$method,m$kind,m$size,r);one$delta<-delta[e];one$target<-truth[e];one$sigma<-sig[e];one$relative_error<-one$se/sig[e]-1;one$rejection<-one$lower>0||one$upper<0;one$coverage<-one$lower<=truth[e]&&one$upper>=truth[e];ans[[k]]<-one}}}
 saveRDS(do.call(rbind,ans),fn);fn}
 paths<-character(R);start<-proc.time()[3]
 for(first in seq(1,R,by=20)){ids<-first:min(R,first+19);paths[ids]<-unlist(parallel_map(ids,worker));cat(sprintf('Synthetic IST %d/%d (%.1f s)\n',max(ids),R,proc.time()[3]-start));flush.console()}
 raw<-do.call(rbind,lapply(paths,readRDS));raw$sigma<-sig[match(raw$delta,delta)];raw$relative_error<-raw$se/raw$sigma-1
 write.csv(raw,file.path(out,'synthetic_replicates.csv'),row.names=FALSE)
}
if(phase%in%c('all','timing')){
 # Serial elapsed time, randomized method order; one warm-up then seven repeats.
 # Input encoding, point estimator, partitions/bags, draws, SD and quantiles included.
 set.seed(cfg$seed+50L);mg<-rbind(grid(n,1000),data.frame(method='ordinary_exact_cells',kind='ordinary',size=n));ans<-list();k<-0L
 bench<-function(kind,size,generic=TRUE){start<-proc.time()[3];calls<-0L
 repeat{cc<-counts_from(a,h,y);th<-ist_estimate(cc);v<-ist_resample(cc,1000,kind,size,generic);ci<-th-c(v['qhi'],v['qlo']);calls<-calls+1L;elapsed<-proc.time()[3]-start;if(elapsed>=.05||calls>=200L)break}
 c(seconds=unname(elapsed/calls),calls=calls,se=unname(v['se']))}
 for(r in 0:7)for(j in sample.int(nrow(mg))){m<-mg[j,];v<-bench(m$kind,m$size,m$method!='ordinary_exact_cells')
 if(r>0){k<-k+1L;ans[[k]]<-data.frame(run=r,method=m$method,seconds=v['seconds'],batch_calls=v['calls'],se=v['se'],row.names=NULL)}}
 write.csv(do.call(rbind,ans),file.path(out,'runtime_serial.csv'),row.names=FALSE)
 ans<-list();k<-0L
 for(r in 0:5)for(kind in sample(c('mn','cbb','blb')))for(size in sample(cfg$tuning)){v<-bench(kind,size,TRUE)
 if(r>0){k<-k+1L;ans[[k]]<-data.frame(run=r,kind=kind,size=size,seconds=v['seconds'],batch_calls=v['calls'])}}
 write.csv(do.call(rbind,ans),file.path(out,'tuning_runtime_serial.csv'),row.names=FALSE)
 cat('Serial timing complete\n');flush.console()
}
writeLines(c(capture.output(str(cfg)),capture.output(sessionInfo()),paste('source MD5',tools::md5sum(data_path))),file.path(out,'sessionInfo.txt'))
if(phase=='all')source(file.path(script_dir,'plot_ist_final.R'))
