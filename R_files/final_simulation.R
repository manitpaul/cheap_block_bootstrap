#!/usr/bin/env Rscript
# Reproducible final study. No result is filtered for favorability.
args <- commandArgs(TRUE)
opt <- function(k,default) {x<-args[startsWith(args,paste0('--',k,'='))];if(length(x))sub(paste0('^--',k,'='),'',tail(x,1)) else default}
script_dir <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
# Some macOS R binaries ship their Fortran runtime but retain an unavailable
# toolchain path in Makeconf. Link the bundled runtime without installing tools.
if(Sys.info()[['sysname']]=='Darwin' && !dir.exists('/opt/gfortran/lib')){
 lib<-file.path(R.home('lib'),'libgfortran.5.dylib');quad<-file.path(R.home('lib'),'libquadmath.0.dylib')
 if(file.exists(lib)&&file.exists(quad)){mv<-file.path(script_dir,'cache','Makevars');dir.create(dirname(mv),showWarnings=FALSE);writeLines(paste('FLIBS =',shQuote(lib),shQuote(quad)),mv);Sys.setenv(R_MAKEVARS_USER=mv)}}
# Process-local compilation avoids cache collisions between concurrent drivers.
Rcpp::sourceCpp(file.path(script_dir,'final_simulation_kernel.cpp'),cacheDir=file.path(tempdir(),'cbb_cpp'))
cores<-as.integer(opt('cores','2'));reps<-as.integer(opt('reps','200'));phase<-opt('phase','all')
Bs<-as.integer(strsplit(opt('B','200,500'),',')[[1]]);deltas<-0:4
cases<-strsplit(opt('cases','log_gamma,inverse_gamma,u_distance,precision'),',')[[1]]
ns_gamma<-as.integer(strsplit(opt('n-gamma','10000,50000'),',')[[1]])
ns_u<-as.integer(strsplit(opt('n-u','1000,4000'),',')[[1]])
ns_precision<-as.integer(strsplit(opt('n-precision','2000,8000'),',')[[1]])
run<-opt('run','final_simulation');out<-file.path(script_dir,'results',run)
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(file.path(out,'checkpoints'),showWarnings=FALSE)
wall<-function()as.numeric(Sys.time())
wait_for_timing_window<-function(){while(file.exists(file.path(script_dir,'PAUSE_MONTE_CARLO')))Sys.sleep(1)}
method_labels<-c(ordinary='Ordinary',mn_log='m/n: m = log n',mn_sqrt='m/n: m = sqrt(n)',cbb_log='CBB: s = log n',cbb_sqrt='CBB: s = sqrt(n)',cbb_B='CBB: s = B',blb_log='BLB: b = log n',blb_sqrt='BLB: b = sqrt(n)')
method_cols<-c(ordinary='#222222',mn_log='#AA4499',mn_sqrt='#CC79A7',cbb_log='#117733',cbb_sqrt='#44AA99',cbb_B='#0072B2',blb_log='#997700',blb_sqrt='#D55E00')
case_labels<-c(log_gamma='Log Gamma mean',inverse_gamma='Capped reciprocal Gamma mean',u_distance='Pairwise Euclidean distance',precision='Regularized precision')
get_ns<-function(ca)if(ca %in% c('log_gamma','inverse_gamma'))ns_gamma else if(ca=='u_distance')ns_u else ns_precision
sizes<-function(n,k){g<-as.integer(floor((0:(n-1))*k/n));list(g=g,sizes=tabulate(g+1,nbins=k))}
nearest_divisor<-function(n,target){lo<-seq_len(floor(sqrt(n)));lo<-lo[n%%lo==0];ds<-sort(unique(c(lo,n/lo)));ds[which.min(abs(log(ds/target)))]}
grid_methods<-function(n,B){stopifnot(n%%min(n,B)==0);data.frame(method=names(method_labels),kind=c('ordinary',rep('mn',2),rep('cbb',3),rep('blb',2)),size=as.integer(c(n,ceiling(log(n)),ceiling(sqrt(n)),nearest_divisor(n,log(n)),nearest_divisor(n,sqrt(n)),min(n,B),ceiling(log(n)),ceiling(sqrt(n)))))}
summarize_errors<-function(err){if(is.null(dim(err)))err<-matrix(err,ncol=1);B<-nrow(err);lo<-floor((B+1)*.025+1e-10);hi<-lo+ceiling((B+1)*.95-1e-10);stopifnot(lo>=1,hi<=B,all(is.finite(err)));rbind(se=apply(err,2,sd),lower=apply(err,2,function(x)sort(x)[lo]),upper=apply(err,2,function(x)sort(x)[hi]))}
gh_rule<-function(k){J<-matrix(0,k,k);J[cbind(1:(k-1),2:k)]<-sqrt(1:(k-1));e<-eigen(J+t(J),symmetric=TRUE);list(x=e$values,w=e$vectors[1,]^2)}
u_constants<-function(k=40){p<-.02;L<-10;c0<-4/sqrt(pi);cl<-(L+2/L)*(2*pnorm(L/sqrt(2))-1)+2/sqrt(pi)*exp(-L^2/4);theta<-(p^2+(1-p)^2)*c0+2*p*(1-p)*cl
 gh<-gh_rule(k);ii<-expand.grid(1:k,1:k,1:k);X<-cbind(gh$x[ii[,1]],gh$x[ii[,2]],gh$x[ii[,3]]);ww<-gh$w[ii[,1]]*gh$w[ii[,2]]*gh$w[ii[,3]]
 f<-function(r)ifelse(r<1e-7,2*sqrt(2/pi),sqrt(2/pi)*exp(-r*r/2)+(r+1/r)*(2*pnorm(r)-1));eq2<-0
 for(g in 0:1){Z<-X;Z[,1]<-Z[,1]+g*L;q<-(1-p)*f(sqrt(rowSums(Z^2)));Z[,1]<-Z[,1]-L;q<-q+p*f(sqrt(rowSums(Z^2)));eq2<-eq2+c(1-p,p)[g+1]*sum(ww*q^2)}
 c(theta=theta,zeta1=eq2-theta^2,varh=6+2*p*(1-p)*L^2-theta^2)}
uc<-u_constants();stopifnot(max(abs(uc-u_constants(60)))<1e-9)
truth_gamma<-function(ca,n,theta=1){kap<-if(ca=='log_gamma').005 else .05;a<-n*kap;if(ca=='log_gamma')return(sqrt(trigamma(a)))
 cut<-a*theta/1000;m1<-1000*pgamma(cut,a)+theta*a/(a-1)*pgamma(cut,a-1,lower.tail=FALSE);m2<-1e6*pgamma(cut,a)+(theta*a)^2/((a-1)*(a-2))*pgamma(cut,a-2,lower.tail=FALSE);sqrt(m2-m1*m1)}
parameters<-function(ca,n){if(ca=='log_gamma'){theta<-deltas/sqrt(n*.005);return(list(theta=theta,alt=theta,null=0,truth=rep(truth_gamma(ca,n),5)))}
 if(ca=='inverse_gamma'){theta<-1+deltas/sqrt(n*.05);return(list(theta=theta,alt=theta,null=1,truth=vapply(theta,function(th)truth_gamma(ca,n,th),0)))}
 if(ca=='u_distance'){sd0<-sqrt(4*(n-2)*uc['zeta1']/(n*(n-1))+2*uc['varh']/(n*(n-1)));alt<-1+deltas*sd0/uc['theta'];return(list(theta=uc['theta']*alt,alt=alt,null=uc['theta'],truth=sd0*alt))}
 theta<-1/1.1+deltas*sqrt(2)/(1.1^2*sqrt(n));list(theta=theta,alt=1/theta-.1,null=1/1.1,truth=NULL)}
make_data<-function(ca,n){if(ca=='log_gamma')return(matrix(rgamma(n,.005,rate=.005),ncol=1));if(ca=='inverse_gamma')return(matrix(rgamma(n,.05,rate=.05),ncol=1));if(ca=='precision')return(matrix(rnorm(n*20),ncol=20));z<-matrix(rnorm(3*n),ncol=3);z[,1]<-z[,1]+10*rbinom(n,1,.02);z}
features<-function(x,ca)if(ca=='precision')precision_features(x) else cbind(mass=1,x[,1])
run_method<-function(x,ca,met,B,pa,h=NULL,profile=FALSE){profile_start<-if(profile)wall()else NA_real_;initial_seconds<-NA_real_;n<-nrow(x);kind<-if(ca=='log_gamma')0L else if(ca=='inverse_gamma')1L else 2L;p<-if(ca=='precision')20L else 1L
 if(ca!='u_distance'){ff<-features(x,ca);estimate<-feature_estimate(ff,kind,p,pa$alt)
  if(met$kind=='blb'){K<-ceiling(n/met$size);ss<-sizes(n,K)$sizes;res<-blb_feature_summary(ff,ss,n,B,kind,p,pa$alt);rownames(res)<-c('se','lower','upper')}
  else {cache<-if(met$kind=='cbb')rowsum(ff,sizes(n,met$size)$g,reorder=FALSE) else ff;stars<-draw_features(cache,met$size,B,kind,p,pa$alt);err<-sweep(stars,2,estimate)*if(met$kind=='mn')sqrt(met$size/n) else 1;res<-summarize_errors(err)}
 }else{cached_a<-NULL
  if(is.null(h)){if(met$kind=='ordinary'){h<-distance_kernel(x);est0<-sum(h)/(n*(n-1))}else{ss0<-if(met$kind=='cbb')sizes(n,met$size)else list(g=integer(0));cache0<-distance_bin_cache(x,ss0$g,if(met$kind=='cbb')met$size else 0L);est0<-cache0$estimate;cached_a<-cache0$table}}
  else est0<-sum(h)/(n*(n-1));estimate<-est0*pa$alt
  if(profile)initial_seconds<-wall()-profile_start
  if(met$kind=='blb'){K<-ceiling(n/met$size);bags<-draw_blb_u(x,sizes(n,K)$sizes,n,B);rr<-Reduce('+',lapply(bags,summarize_errors))/K;res<-rr %o% pa$alt;dim(res)<-c(3,length(pa$alt));rownames(res)<-c('se','lower','upper')}
  else {if(met$kind=='mn')star<-draw_u_small(x,met$size,B) else if(met$kind=='ordinary')star<-draw_u(h,rep.int(1L,n),n,B) else {ss<-sizes(n,met$size);a<-if(is.null(cached_a))compress_kernel(h,ss$g,met$size)else cached_a;star<-draw_u(a,ss$sizes,met$size,B)}
  err<-(star-est0)*if(met$kind=='mn')sqrt(met$size/n) else 1;res<-summarize_errors(err %o% pa$alt)}
 }
 list(estimate=estimate,se=res['se',],lower=estimate-res['upper',],upper=estimate-res['lower',],initial_seconds=initial_seconds,remainder_seconds=if(profile&&ca=='u_distance')wall()-profile_start-initial_seconds else NA_real_)}
RNGkind("L'Ecuyer-CMRG");set.seed(20261118)
streams_for<-function(n,seed){set.seed(seed);s<-vector('list',n);for(i in seq_len(n)){s[[i]]<-.Random.seed;.Random.seed<-parallel::nextRNGStream(.Random.seed)};s}
par_run<-function(ix,fun){if(cores==1)return(lapply(ix,fun));z<-parallel::mclapply(ix,fun,mc.cores=cores,mc.preschedule=TRUE,mc.set.seed=FALSE);if(any(vapply(z,inherits,FALSE,'try-error')))stop('worker failure; no failed runs discarded');z}
precision_reference<-function(n,R=20000){path<-file.path(out,paste0('precision_reference_n',n,'.rds'));if(file.exists(path))return(readRDS(path));set.seed(999000+n);pa<-parameters('precision',n);values<-matrix(NA_real_,R,5)
 for(r in seq_len(R)){S<-rWishart(1,n-1,diag(20))[,,1]/n;q<-S[1,1]-sum(S[1,-1]*solve(S[-1,-1]+diag(.1,19),S[-1,1]));values[r,]<-1/(pa$alt*q+.1)}
 z<-list(n=n,reps=R,truth=apply(values,2,sd),means=colMeans(values),mcse=apply(values,2,sd)/sqrt(2*(R-1)),values=values);saveRDS(z,path);z}
if('precision'%in%cases){refs<-setNames(lapply(ns_precision,precision_reference),ns_precision)}else refs<-list()
cfg<-list(seed=if(startsWith(phase,'tuning'))20261120L else 20261118L,reps=reps,B=Bs,cases=cases,n_gamma=ns_gamma,n_u=ns_u,n_precision=ns_precision,deltas=deltas,blb_inner='B per bag',uc=uc,version=2L,precision_alternative='diag(tau,1,...,1)')
if(file.exists(file.path(out,'config.rds'))&&!identical(readRDS(file.path(out,'config.rds')),cfg))stop('Existing run has a different configuration; choose a new --run.')
saveRDS(cfg,file.path(out,'config.rds'))
if(phase%in%c('all','pilot','statistics')){
 grid<-do.call(rbind,lapply(cases,function(ca)expand.grid(rep=seq_len(reps),n=get_ns(ca),case=ca,stringsAsFactors=FALSE)));streams<-streams_for(nrow(grid),20261118L)
 work<-function(i){path<-file.path(out,'checkpoints',sprintf('task_%06d.rds',i));if(file.exists(path))return(path);assign('.Random.seed',streams[[i]],.GlobalEnv);g<-grid[i,];n<-g$n;ca<-g$case;x<-make_data(ca,n);x<-x[sample.int(n),,drop=FALSE];pa<-parameters(ca,n);if(ca=='precision')pa$truth<-refs[[as.character(n)]]$truth
 h<-if(ca=='u_distance')distance_kernel(x) else NULL;ans<-list();z<-0
 for(B in Bs){methods<-grid_methods(n,B);for(j in sample.int(nrow(methods))){met<-methods[j,];start<-wall();a<-run_method(x,ca,met,B,pa,h);elapsed<-wall()-start;z<-z+1
 ans[[z]]<-data.frame(case=ca,n=n,rep=g$rep,B=B,delta=deltas,method=met$method,size=met$size,K=if(met$kind=='blb')ceiling(n/met$size) else 1,total_draws=if(met$kind=='blb')B*ceiling(n/met$size) else B,theta=pa$theta,estimate=a$estimate,se=a$se,truth=pa$truth,relative_error=a$se/pa$truth-1,lower=a$lower,upper=a$upper,rejection=pa$null<a$lower|pa$null>a$upper,coverage=a$lower<=pa$theta&a$upper>=pa$theta,elapsed_parallel=elapsed)} }
 saveRDS(do.call(rbind,ans),path);path}
 start<-wall();paths<-character(nrow(grid));for(first in seq(1,nrow(grid),by=cores*2)){wait_for_timing_window();ix<-first:min(nrow(grid),first+cores*2-1);paths[ix]<-unlist(par_run(ix,work));cat(sprintf('%s %d/%d %.1fs\n',run,max(ix),nrow(grid),wall()-start));flush.console()}
 raw<-do.call(rbind,lapply(paths,readRDS));write.csv(raw,file.path(out,'replicates.csv'),row.names=FALSE)
 summary<-do.call(rbind,lapply(split(raw,interaction(raw$case,raw$n,raw$B,raw$delta,raw$method,drop=TRUE)),function(z){r<-sqrt(mean(z$relative_error^2));data.frame(z[1,c('case','n','B','delta','method','size','K','total_draws')],reps=nrow(z),rejection=mean(z$rejection),rejection_mcse=sqrt(mean(z$rejection)*(1-mean(z$rejection))/nrow(z)),coverage=mean(z$coverage),relative_bias=mean(z$relative_error),relative_rmse=r,rmse_mcse=sd(z$relative_error^2)/(2*r*sqrt(nrow(z))),mean_length=mean(z$upper-z$lower))}));write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE)
 cat('Statistics complete\n');flush.console()
}
if(phase%in%c('all','timing')){set.seed(20261119);ans<-list();z<-0;times<-as.integer(opt('timing-reps','5'));for(ca in cases)for(n in get_ns(ca)){x<-make_data(ca,n);pa<-parameters(ca,n);pa$alt<-pa$alt[1]
 for(r in 0:times)for(B in Bs){methods<-grid_methods(n,B);for(j in sample.int(nrow(methods))){met<-methods[j,];start<-wall();a<-run_method(x,ca,met,B,pa,profile=TRUE);elapsed<-wall()-start;if(r>0){z<-z+1;ans[[z]]<-data.frame(case=ca,n=n,B=B,method=met$method,repeat_id=r,seconds=elapsed,initial_seconds=a$initial_seconds,remainder_seconds=a$remainder_seconds)}}};cat('Timing',ca,n,'complete\n');flush.console()}
 timing<-do.call(rbind,ans);write.csv(timing,file.path(out,'runtime_serial.csv'),row.names=FALSE);ss<-do.call(rbind,lapply(split(timing,interaction(timing$case,timing$n,timing$B,timing$method,drop=TRUE)),function(z)data.frame(z[1,c('case','n','B','method')],seconds=median(z$seconds),q25=unname(quantile(z$seconds,.25)),q75=unname(quantile(z$seconds,.75)),initial_seconds=median(z$initial_seconds),remainder_seconds=median(z$remainder_seconds))));write.csv(ss,file.path(out,'runtime_serial_summary.csv'),row.names=FALSE)}
if(phase=='tuning'){
 stopifnot(all(cases%in%c('log_gamma','inverse_gamma')));n<-max(ns_gamma);B<-max(Bs)
 controls<-sort(unique(as.integer(c(ceiling(log(n)),25,50,100,ceiling(sqrt(n)),500,1000))))
 tuning_methods<-rbind(data.frame(method='ordinary',kind='ordinary',size=n,requested=n),do.call(rbind,lapply(c('mn','cbb','blb'),function(k)data.frame(method=k,kind=k,size=if(k=='cbb')vapply(controls,function(s)nearest_divisor(n,s),0) else controls,requested=controls))))
 tuning_methods<-tuning_methods[!duplicated(tuning_methods[c('kind','size')]),]
 grid<-expand.grid(rep=seq_len(reps),case=cases,stringsAsFactors=FALSE);streams<-streams_for(nrow(grid),20261120L)
 work<-function(i){path<-file.path(out,'checkpoints',sprintf('tuning_%06d.rds',i));if(file.exists(path))return(path);assign('.Random.seed',streams[[i]],.GlobalEnv);g<-grid[i,];ca<-g$case;x<-make_data(ca,n);x<-x[sample.int(n),,drop=FALSE];pa<-parameters(ca,n);ans<-list()
 for(j in sample.int(nrow(tuning_methods))){met<-tuning_methods[j,];a<-run_method(x,ca,met,B,pa);ans[[j]]<-data.frame(case=ca,n=n,B=B,rep=g$rep,delta=deltas,method=met$method,size=met$size,requested=met$requested,relative_error=a$se/pa$truth-1,rejection=pa$null<a$lower|pa$null>a$upper,coverage=a$lower<=pa$theta&a$upper>=pa$theta)}
 saveRDS(do.call(rbind,ans),path);path}
 paths<-character(nrow(grid));start<-wall();for(first in seq(1,nrow(grid),by=cores*2)){wait_for_timing_window();ix<-first:min(nrow(grid),first+cores*2-1);paths[ix]<-unlist(par_run(ix,work));cat(sprintf('Tuning %d/%d %.1fs\n',max(ix),nrow(grid),wall()-start));flush.console()}
 raw<-do.call(rbind,lapply(paths,readRDS));write.csv(raw,file.path(out,'tuning_replicates.csv'),row.names=FALSE)
 ss<-do.call(rbind,lapply(split(raw,interaction(raw$case,raw$delta,raw$method,raw$size,drop=TRUE)),function(z){r<-sqrt(mean(z$relative_error^2));data.frame(z[1,c('case','n','B','delta','method','size','requested')],reps=nrow(z),rejection=mean(z$rejection),coverage=mean(z$coverage),relative_rmse=r,rmse_mcse=sd(z$relative_error^2)/(2*r*sqrt(nrow(z))),relative_bias=mean(z$relative_error))}));write.csv(ss,file.path(out,'tuning_summary.csv'),row.names=FALSE)
}
if(phase=='tuning-timing'){
 n<-max(ns_gamma);B<-max(Bs);controls<-sort(unique(as.integer(c(ceiling(log(n)),25,50,100,ceiling(sqrt(n)),500,1000))))
 mm<-rbind(data.frame(method='ordinary',kind='ordinary',size=n),do.call(rbind,lapply(c('mn','cbb','blb'),function(k)data.frame(method=k,kind=k,size=if(k=='cbb')vapply(controls,function(s)nearest_divisor(n,s),0) else controls))));mm<-mm[!duplicated(mm[c('kind','size')]),]
 set.seed(20261121);ans<-list();k<-0;for(ca in cases){x<-make_data(ca,n);pa<-parameters(ca,n);pa$alt<-pa$alt[1];for(r in 0:as.integer(opt('timing-reps','5')))for(j in sample.int(nrow(mm))){met<-mm[j,];start<-wall();invisible(run_method(x,ca,met,B,pa));elapsed<-wall()-start;if(r>0){k<-k+1;ans[[k]]<-data.frame(case=ca,n=n,B=B,method=met$method,size=met$size,repeat_id=r,seconds=elapsed)}}}
 tt<-do.call(rbind,ans);write.csv(tt,file.path(out,'tuning_runtime_serial.csv'),row.names=FALSE);ts<-do.call(rbind,lapply(split(tt,interaction(tt$case,tt$method,tt$size,drop=TRUE)),function(z)data.frame(z[1,c('case','n','B','method','size')],seconds=median(z$seconds),q25=unname(quantile(z$seconds,.25)),q75=unname(quantile(z$seconds,.75)))));write.csv(ts,file.path(out,'tuning_runtime_summary.csv'),row.names=FALSE)
}
session_text<-c(paste('Phase',phase),paste('Workers',cores),capture.output(sessionInfo()))
writeLines(session_text,file.path(out,paste0('sessionInfo_',phase,'.txt')))
if(!file.exists(file.path(out,'sessionInfo.txt')))writeLines(session_text,file.path(out,'sessionInfo.txt'))
