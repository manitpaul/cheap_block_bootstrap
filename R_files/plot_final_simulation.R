#!/usr/bin/env Rscript
args<-commandArgs(TRUE);statistics_only<-'--statistics-only'%in%args;args<-args[args!='--statistics-only'];f<-grep('^--file=',commandArgs(),value=TRUE);script_dir<-dirname(normalizePath(sub('^--file=','',f[1])))
run<-if(length(args))args[1] else 'final_simulation';out<-file.path(script_dir,'results',run);fig<-file.path(dirname(script_dir),'PDF_images');dir.create(fig,showWarnings=FALSE)
s<-read.csv(file.path(out,'summary.csv'));if(!statistics_only){t<-read.csv(file.path(out,'runtime_serial_summary.csv'));z<-merge(s,t,by=c('case','n','B','method'))}
methods<-c('ordinary','mn_log','mn_sqrt','cbb_log','cbb_sqrt','cbb_B','blb_log','blb_sqrt')
labels<-c('Ordinary','m/n: m = log n','m/n: m = sqrt(n)','CBB: s = log n','CBB: s = sqrt(n)','CBB: s = B','BLB: b = log n','BLB: b = sqrt(n)')
cols<-c('#222222','#AA4499','#CC79A7','#117733','#44AA99','#0072B2','#997700','#D55E00');ltys<-c(1,2,3,2,3,1,2,3);pchs<-c(16,17,15,1,2,18,0,5)
cases<-c('log_gamma','inverse_gamma','u_distance','precision');titles<-c('Log Gamma mean','Capped reciprocal Gamma mean','Pairwise Euclidean distance','Regularized precision')
fmt<-function(x)format(x,big.mark=',',scientific=FALSE,trim=TRUE)
interval<-function(x0,y0,x1,y1,col){span<-sqrt((grconvertX(x1,'user','inches')-grconvertX(x0,'user','inches'))^2+(grconvertY(y1,'user','inches')-grconvertY(y0,'user','inches'))^2);j<-span>=.001;if(any(j))arrows(x0[j],y0[j],x1[j],y1[j],angle=90,code=3,length=.025,col=col);if(any(!j))segments(x0[!j],y0[!j],x1[!j],y1[!j],col=col)}
openplot<-function(name){pdf(file.path(fig,name),width=9,height=7.8,pointsize=12,family='Helvetica',useDingbats=FALSE);layout(matrix(c(1,2,3,4,5,5),nrow=3,byrow=TRUE),heights=c(1,1,.38));par(mar=c(3.8,5.1,3.8,1.1),oma=c(0,0,0,0),mgp=c(2.8,.75,0),tcl=-.25,las=1,cex=.95,cex.main=1.1,cex.lab=1,cex.axis=1,font.main=1,xpd=FALSE)}
endplot<-function(note=NULL){par(mar=c(0,0,0,0));plot.new();legend(.5,if(is.null(note)).8 else .95,xjust=.5,yjust=1,legend=labels,col=cols,lty=ltys,pch=pchs,ncol=3,bty='n',cex=.95,x.intersp=.7,y.intersp=1.15);if(!is.null(note))text(.5,.08,note,cex=.95);dev.off()}
linepanel<-function(data,x,y,title,xlab,ylab,log='',ylim=NULL,xat=NULL,ref=NULL,ci=FALSE){
 if(is.null(ylim)){ylim<-range(data[[y]],finite=TRUE)*c(.75,1.3);if(y=='relative_rmse')ylim<-range(c(ylim,pmax(data[[y]]-1.96*data$rmse_mcse,data[[y]]*.1),data[[y]]+1.96*data$rmse_mcse))}
 xlim<-if(x=='seconds')range(data$q25,data$q75)*c(.8,1.2)else range(data[[x]])*c(.85,1.15)
 plot(xlim,ylim,type='n',log=log,xlab=xlab,ylab=ylab,main='',xaxt=if(is.null(xat))'s'else'n',yaxt=if(grepl('y',log))'n'else's')
 mtext(title,side=3,line=1.1,cex=1.05,font=1)
 if(grepl('y',log)){yt<-axTicks(2);axis(2,at=yt,labels=vapply(yt,function(v)format(v,scientific=FALSE,trim=TRUE),''))}
 if(!is.null(xat))axis(1,at=xat,labels=fmt(xat));if(!is.null(ref))abline(h=ref,lty=2,col='gray65')
 for(j in seq_along(methods)){
  a<-data[data$method==methods[j],];a<-a[order(a[[x]]),];ec<-adjustcolor(cols[j],alpha.f=.6)
  if(ci){zz<-qnorm(.975);den<-1+zz^2/a$reps;ctr<-(a[[y]]+zz^2/(2*a$reps))/den;half<-zz*sqrt(a[[y]]*(1-a[[y]])/a$reps+zz^2/(4*a$reps^2))/den;interval(a[[x]],ctr-half,a[[x]],ctr+half,ec)}
  if(y=='relative_rmse'&&'rmse_mcse'%in%names(a))interval(a[[x]],pmax(a[[y]]-1.96*a$rmse_mcse,.Machine$double.eps),a[[x]],a[[y]]+1.96*a$rmse_mcse,ec)
  if(x=='seconds'&&all(c('q25','q75')%in%names(a)))interval(a$q25,a[[y]],a$q75,a[[y]],ec)
  lines(a[[x]],a[[y]],col=cols[j],lty=ltys[j],pch=pchs[j],type='b',lwd=1.5,cex=.8)
 }
}
openplot('final_calibration.pdf')
for(i in seq_along(cases)){a<-s[s$case==cases[i]&s$n==max(s$n[s$case==cases[i]])&s$delta==0,];linepanel(a,'B','rejection',paste0(titles[i],'\nn = ',fmt(a$n[1])),'Bootstrap draws B','Null rejection probability',ylim=c(0,1),xat=sort(unique(a$B)),ref=.05,ci=TRUE)}
endplot()
openplot('final_power.pdf')
for(i in seq_along(cases)){a<-s[s$case==cases[i]&s$n==max(s$n[s$case==cases[i]])&s$B==max(s$B),];linepanel(a,'delta','rejection',paste0(titles[i],'\nn = ',fmt(a$n[1]),', B = ',a$B[1]),'Local effect d','Rejection probability',ylim=c(0,1),xat=0:4,ref=.05)}
endplot()
openplot('final_se_rmse.pdf')
for(i in seq_along(cases)){a<-s[s$case==cases[i]&s$B==max(s$B)&s$delta==0,];linepanel(a,'n','relative_rmse',paste0(titles[i],'\nB = ',a$B[1]),'Sample size n','Relative RMSE of\nstandard-error estimate',log='xy',xat=sort(unique(a$n)))}
endplot()
if(statistics_only){cat('Three statistical figures written; runtime figures await isolated timings.\n');quit(save='no')}
openplot('final_runtime_rmse.pdf')
for(i in seq_along(cases)){a<-z[z$case==cases[i]&z$n==max(z$n[z$case==cases[i]])&z$delta==0,];linepanel(a,'seconds','relative_rmse',paste0(titles[i],'\nn = ',fmt(a$n[1])),'Elapsed time (seconds)','Relative RMSE of\nstandard-error estimate',log='xy')}
endplot('Connected points: B = 200 and 500. Horizontal bars: timing IQR.')
openplot('final_runtime_power.pdf')
for(i in seq_along(cases)){a<-z[z$case==cases[i]&z$n==max(z$n[z$case==cases[i]])&z$delta==2,];linepanel(a,'seconds','rejection',paste0(titles[i],'\nn = ',fmt(a$n[1]),', d = 2'),'Elapsed time (seconds)','Rejection probability',log='x',ylim=c(0,1),ref=.05,ci=TRUE)}
endplot('Connected points: B = 200 and 500. Interpret power alongside calibration.')
cat('Five main figures written.\n')
