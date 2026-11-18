# Sourced by real_data_ist_final.R after all computation/timing has completed.
ob<-read.csv(file.path(out,'observed_runs.csv'));tu<-read.csv(file.path(out,'observed_tuning.csv'))
tm<-read.csv(file.path(out,'runtime_serial.csv'));tt<-read.csv(file.path(out,'tuning_runtime_serial.csv'))
raw<-read.csv(file.path(out,'synthetic_replicates.csv'));ref<-readRDS(file.path(out,'observed_reference.rds'))
reference_se<-sd(ref);referr<-sort(ref-theta);lo<-floor((length(ref)+1)*.025+1e-10);hi<-lo+ceiling((length(ref)+1)*.95-1e-10);refci<-theta-referr[c(hi,lo)]
tu$endpoint_error<-sqrt(((tu$lower-refci[1])^2+(tu$upper-refci[2])^2)/2)/reference_se
ob$endpoint_error<-sqrt(((ob$lower-refci[1])^2+(ob$upper-refci[2])^2)/2)/reference_se
wilson<-function(x){N<-length(x);p<-mean(x);zz<-qnorm(.975);den<-1+zz^2/N;c((p+zz^2/(2*N)-zz*sqrt(p*(1-p)/N+zz^2/(4*N^2)))/den,(p+zz^2/(2*N)+zz*sqrt(p*(1-p)/N+zz^2/(4*N^2)))/den)}
groups<-split(raw,interaction(raw$method,raw$B,raw$delta,drop=TRUE))
ss<-do.call(rbind,lapply(groups,function(x){v<-sqrt(mean(x$relative_error^2));mc<-sd(x$relative_error^2)/(2*v*sqrt(nrow(x)));w<-wilson(x$rejection);cov<-wilson(x$coverage);data.frame(method=x$method[1],B=x$B[1],delta=x$delta[1],reps=nrow(x),rejection=mean(x$rejection),rejection_low=w[1],rejection_high=w[2],coverage=mean(x$coverage),coverage_low=cov[1],coverage_high=cov[2],rrmse=v,rrmse_low=pmax(0,v-1.96*mc),rrmse_high=v+1.96*mc,relative_bias=mean(x$relative_error))}))
write.csv(ss,file.path(out,'synthetic_summary.csv'),row.names=FALSE)
obsum<-do.call(rbind,lapply(names(labels),function(m){x<-ob[ob$method==m,];t<-tm[tm$method==m,];data.frame(method=m,mean_se=mean(x$se),lower=median(x$lower),upper=median(x$upper),se_rrmse=sqrt(mean(x$relative_error^2)),ci_rrmse=sqrt(mean(x$endpoint_error^2)),runtime=median(t$seconds),runtime_q25=quantile(t$seconds,.25),runtime_q75=quantile(t$seconds,.75))}))
write.csv(obsum,file.path(out,'observed_summary.csv'),row.names=FALSE)
tusum<-do.call(rbind,lapply(split(tu,interaction(tu$kind,tu$size,drop=TRUE)),function(x){tim<-tt[tt$kind==x$kind[1]&tt$size==x$size[1],];data.frame(kind=x$kind[1],size=x$size[1],se_rrmse=sqrt(mean(x$relative_error^2)),ci_rrmse=sqrt(mean(x$endpoint_error^2)),runtime=median(tim$seconds),mean_se=mean(x$se),mean_lower=mean(x$lower),mean_upper=mean(x$upper))}))
write.csv(tusum,file.path(out,'tuning_summary.csv'),row.names=FALSE)
figure<-function(name,w=11.5,h=8.0){pdf(file.path(root,'PDF_images',name),width=w,height=h,pointsize=16,family='Helvetica',useDingbats=FALSE);par(mgp=c(2.5,.7,0),tcl=-.25,las=1,cex=.95)}
short<-c('Ordinary','m/n: m=10','m/n: m=139','CBB: s=10','CBB: s=139','CBB: s=1000','BLB: b~10','BLB: b~139')
ord<-names(labels);yy<-8:1
hbase<-function(xlim,title,xlab,log=''){plot(NA,xlim=xlim,ylim=c(.5,8.5),yaxt='n',ylab='',xlab=xlab,main=title,log=log);axis(2,at=yy,labels=short,las=1,cex.axis=1);abline(h=yy,col='#eeeeee')}
figure('ist_observed_final.pdf',12,8.2);par(mfrow=c(2,2),mar=c(4.2,7.2,2.8,1),cex=.95)
hbase(range(100*c(obsum$lower,obsum$upper,refci)) + c(-.15,.15),'A  Observed 95% basic intervals','Aspirin risk difference (percentage points)')
abline(v=0,lty=2,col='gray60');abline(v=100*refci,lty=3,col='gray60')
for(j in 1:8){segments(100*obsum$lower[j],yy[j],100*obsum$upper[j],yy[j],col=cols[j],lwd=2);points(100*theta,yy[j],col=cols[j],pch=pchs[j])}
hbase(range(obsum$se_rrmse)*100*c(.8,1.2),'B  Observed SE agreement','RMS relative SE disagreement (%)',log='x')
for(j in 1:8)points(100*obsum$se_rrmse[j],yy[j],pch=pchs[j],col=cols[j],cex=1.1)
exact_time<-median(tm$seconds[tm$method=='ordinary_exact_cells']);write.csv(data.frame(generic_seconds=obsum$runtime[1],exact_cell_seconds=exact_time),file.path(out,'exact_cell_timing.csv'),row.names=FALSE)
hbase(range(c(pmax(.0001,obsum$runtime_q25),obsum$runtime_q75,exact_time))*c(.7,1.5),'C  Observed-data computation','Serial time (seconds; log scale)',log='x')
abline(v=exact_time,lty=3,col='gray45');for(j in 1:8){segments(pmax(.0001,obsum$runtime_q25[j]),yy[j],obsum$runtime_q75[j],yy[j],col=cols[j],lwd=2);points(obsum$runtime[j],yy[j],pch=pchs[j],col=cols[j])}
hbase(range(obsum$ci_rrmse)*100*c(.8,1.2),'D  Interval-endpoint agreement','Endpoint RMS / reference SE (%)',log='x')
for(j in 1:8)points(100*obsum$ci_rrmse[j],yy[j],pch=pchs[j],col=cols[j],cex=1.1)
dev.off()

figure('ist_synthetic_final.pdf',12,8.5)
layout(matrix(c(1,2,3,4,5,5),3,2,byrow=TRUE),heights=c(1,1,.3));par(mar=c(4.3,5.1,2.5,1),cex=.95)
offset<-seq(-.28,.28,length.out=8)
null<-ss[ss$delta==0,]
plot(NA,xlim=c(.6,2.4),ylim=c(0,max(null$rejection_high)*1.05),xaxt='n',xlab='Bootstrap draws B',ylab='Null rejection probability',main='A  Synthetic null calibration')
axis(1,1:2,c(200,1000));abline(h=.05,lty=2,col='gray50')
for(j in 1:8){x<-null[null$method==ord[j],];x<-x[order(x$B),];xx<-1:2+offset[j];segments(xx,x$rejection_low,xx,x$rejection_high,col=cols[j]);lines(xx,x$rejection,col=cols[j],pch=pchs[j],lty=ltys[j],type='b')}
plot(NA,xlim=c(0,4),ylim=c(0,1),xlab='Standardized effect d',ylab='Rejection probability',main='B  Synthetic power (B = 1000)')
abline(h=.05,lty=2,col='gray50');for(j in 1:8){x<-ss[ss$method==ord[j]&ss$B==1000,];x<-x[order(x$delta),];lines(x$delta,x$rejection,col=cols[j],pch=pchs[j],lty=ltys[j],type='b',lwd=1.5)}
plot(NA,xlim=c(.6,2.4),ylim=range(c(null$rrmse_low,null$rrmse_high))*c(.85,1.1),xaxt='n',xlab='Bootstrap draws B',ylab='Relative RMSE of\nstandard-error estimate',main='C  Sampling SE accuracy',log='y')
axis(1,1:2,c(200,1000));for(j in 1:8){x<-null[null$method==ord[j],];x<-x[order(x$B),];xx<-1:2+offset[j];segments(xx,x$rrmse_low,xx,x$rrmse_high,col=cols[j]);lines(xx,x$rrmse,col=cols[j],pch=pchs[j],lty=ltys[j],type='b')}
cov<-ss[ss$B==1000,];plot(NA,xlim=c(0,4),ylim=c(min(cov$coverage_low)*.98,1),xlab='Standardized effect d',ylab='Coverage probability',main='D  Coverage (B = 1000)');abline(h=.95,lty=2,col='gray50')
for(j in 1:8){x<-cov[cov$method==ord[j],];x<-x[order(x$delta),];lines(x$delta,x$coverage,col=cols[j],pch=pchs[j],lty=ltys[j],type='b',lwd=1.5)}
par(mar=c(0,0,0,0));plot.new();legend('center',legend=labels,col=cols,pch=pchs,lty=ltys,ncol=4,bty='n',cex=1)
dev.off()

figure('ist_tuning_final.pdf',11.5,5.1);layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),heights=c(1,.32));par(mar=c(4.5,5.1,2.5,1),cex=.95)
tc<-c(mn='#AA4499',cbb='#0072B2',blb='#D55E00');tp<-c(mn=17,cbb=18,blb=5);tl<-c(mn=2,cbb=1,blb=3)
plot(NA,xlim=range(tusum$size),ylim=range(100*c(tusum$se_rrmse,obsum$se_rrmse[1]))*c(.8,1.1),log='xy',xlab='m (m/n), s (CBB), or b (BLB)',ylab='RMS relative SE\ndisagreement (%)',main='A  Observed tuning (B = 1000)')
abline(h=100*obsum$se_rrmse[1],lty=2,col='#222222')
for(kind in names(tc)){x<-tusum[tusum$kind==kind,];x<-x[order(x$size),];lines(x$size,100*x$se_rrmse,type='b',pch=tp[kind],lty=tl[kind],col=tc[kind],lwd=1.6)}
plot(NA,xlim=range(c(tusum$runtime,obsum$runtime[1],exact_time))*c(.7,1.2),ylim=range(100*c(tusum$se_rrmse,obsum$se_rrmse[1]))*c(.8,1.1),log='xy',xlab='Serial elapsed time (seconds)',ylab='RMS relative SE\ndisagreement (%)',main='B  Accuracy and computation')
for(kind in names(tc)){x<-tusum[tusum$kind==kind,];x<-x[order(x$size),];lines(x$runtime,100*x$se_rrmse,type='b',pch=tp[kind],lty=tl[kind],col=tc[kind],lwd=1.6)}
points(obsum$runtime[1],100*obsum$se_rrmse[1],pch=16,cex=1.3);points(exact_time,100*obsum$se_rrmse[1],pch=6,cex=1.3)
par(mar=c(0,0,0,0));plot.new();legend('center',legend=c('m/n','Ordinary (row)','CBB','Ordinary (exact cells)','BLB',''),col=c(tc['mn'],'#222222',tc['cbb'],'#222222',tc['blb'],NA),pch=c(tp['mn'],16,tp['cbb'],6,tp['blb'],NA),lty=c(tl['mn'],NA,tl['cbb'],NA,tl['blb'],NA),ncol=3,bty='n',cex=1)
dev.off()
getobs<-function(m)obsum[obsum$method==m,]
oo<-getobs('ordinary');cb<-getobs('cbb_B');mm<-getobs('mn_sqrt');bb<-getobs('blb_sqrt');small<-getobs('blb_log')
observed_text<-sprintf(paste0(
 'The reference standard error is $%.3f$ percentage points. Ordinary bootstrap and cheap bootstrap with $s=B=1000$ give mean standard errors $%.3f$ and $%.3f$ percentage points, respectively; their median basic intervals are $[%.2f,%.2f]$ and $[%.2f,%.2f]$ percentage points. Thus the cheap method closely reproduces the ordinary-bootstrap uncertainty estimate and the same qualitative conclusion. Its observed SE disagreement is $%.1f\\%%$, compared with $%.1f\\%%$ for ordinary bootstrap. Cheap bootstrap takes $%.4f$ seconds per generic serial run, compared with $%.4f$ for ordinary bootstrap, a $%.1f$-fold reduction.\n\n',
 'The $m=139$ and $b\\approx139$ alternatives are also competitive in agreement, with relative discrepancies $%.1f\\%%$ and $%.1f\\%%$. Compression must be sufficient: all three methods deteriorate at their smallest size, and the $b\\approx10$ BLB standard error is only $%.1f\\%%$ of the reference on average. The exact-cell ordinary bootstrap takes $%.4f$ seconds, faster than cheap bootstrap with $s=B=1000$ for this particular estimator. The clinical example therefore illustrates accurate general-purpose compression and its tuning trade-off, while also showing why an available exact sufficient-summary implementation should be used as a benchmark.'),
 100*reference_se,100*oo$mean_se,100*cb$mean_se,100*oo$lower,100*oo$upper,100*cb$lower,100*cb$upper,100*cb$se_rrmse,100*oo$se_rrmse,cb$runtime,oo$runtime,oo$runtime/cb$runtime,100*mm$se_rrmse,100*bb$se_rrmse,100*small$mean_se/reference_se,exact_time)
writeLines(observed_text,file.path(root,'ist_results_final.tex'))
pick<-function(m,d=0)ss[ss$method==m&ss$B==1000&ss$delta==d,]
oo<-pick('ordinary');cb<-pick('cbb_B');bb<-pick('blb_log');mml<-pick('mn_log');cbl<-pick('cbb_log')
op<-pick('ordinary',2);cp<-pick('cbb_B',2);mp<-pick('mn_sqrt',2);bp<-pick('blb_sqrt',2)
synthetic_text<-sprintf(paste0(
 'At $B=1000$, ordinary bootstrap and cheap bootstrap with $s=B$ reject $%.1f\\%%$ and $%.1f\\%%$ of null samples and have powers $%.1f\\%%$ and $%.1f\\%%$ at $d=2$. Their null SE relative RMSEs are $%.1f\\%%$ and $%.1f\\%%$. These differences should be read with the Monte Carlo intervals: with 300 data sets, the standard error of a rejection frequency near one half is about 2.9 percentage points. The $m=139$ and $b\\approx139$ methods have powers $%.1f\\%%$ and $%.1f\\%%$ at this alternative, so the application does not establish uniform accuracy or power dominance over either competitor.\n\n',
 'The smallest sizes expose a different limitation. At $B=1000$, null rejection is $%.1f\\%%$ for $m=10$, $%.1f\\%%$ for $s=10$, and $%.1f\\%%$ for $b\\approx10$. Their apparently larger rejection frequencies under alternatives must therefore be assessed against their null miscalibration. Increasing the number of bootstrap draws reduces Monte Carlo variability but cannot repair inadequate subset support or excessive block compression. The observed-data intervals and synthetic experiment together support $s$ large enough to preserve inferential accuracy, with $s=B$ providing an effective choice in this example.'),
 100*oo$rejection,100*cb$rejection,100*op$rejection,100*cp$rejection,100*oo$rrmse,100*cb$rrmse,100*mp$rejection,100*bp$rejection,100*mml$rejection,100*cbl$rejection,100*bb$rejection)
writeLines(synthetic_text,file.path(root,'ist_synthetic_results_final.tex'))
cat('IST figures and summaries written\n')
