#!/usr/bin/env Rscript
args<-commandArgs(TRUE);script_dir<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
run<-if(length(args))args[1]else'final_tuning';out<-file.path(script_dir,'results',run);fig<-file.path(dirname(script_dir),'PDF_images')
s<-read.csv(file.path(out,'tuning_summary.csv'));t<-read.csv(file.path(out,'tuning_runtime_summary.csv'));z<-merge(s,t,by=c('case','n','B','method','size'))
cols<-c(ordinary='#222222',mn='#AA4499',cbb='#0072B2',blb='#D55E00');pchs<-c(ordinary=16,mn=17,cbb=18,blb=5);lty<-c(ordinary=1,mn=2,cbb=1,blb=3)
interval<-function(x0,y0,x1,y1,col){span<-sqrt((grconvertX(x1,'user','inches')-grconvertX(x0,'user','inches'))^2+(grconvertY(y1,'user','inches')-grconvertY(y0,'user','inches'))^2);j<-span>=.001;if(any(j))arrows(x0[j],y0[j],x1[j],y1[j],angle=90,code=3,length=.025,col=col);if(any(!j))segments(x0[!j],y0[!j],x1[!j],y1[!j],col=col)}
draw<-function(file,delta,y,ylab){pdf(file.path(fig,file),width=9,height=4.8,pointsize=12,family='Helvetica',useDingbats=FALSE);layout(matrix(c(1,2,3,3),2,byrow=TRUE),heights=c(1,.27));par(mar=c(4.1,5.1,3.8,1),oma=c(0,0,0,0),mgp=c(2.8,.8,0),tcl=-.25,las=1,cex=.95,cex.main=1.1,cex.lab=1,cex.axis=1,font.main=1,xpd=FALSE)
 for(ca in c('log_gamma','inverse_gamma')){a<-z[z$case==ca&z$delta==delta,];yl<-if(y=='rejection')c(0,1)else range(a[[y]])*c(.65,1.3);plot(range(a$seconds)*c(.8,1.25),yl,type='n',log=if(y=='rejection')'x'else'xy',xlab='Elapsed time (seconds)',ylab=ylab,yaxt=if(y=='relative_rmse')'n'else's',main='');mtext(paste0(if(ca=='log_gamma')'Log Gamma mean'else'Capped reciprocal Gamma mean','\nn = ',format(a$n[1],big.mark=','),', B = ',a$B[1],if(delta>0)paste0(', d = ',delta)else''),side=3,line=1.35,cex=1.05,font=1)
  if(y=='relative_rmse'){yt<-axTicks(2);axis(2,at=yt,labels=vapply(yt,function(v)format(v,scientific=FALSE,trim=TRUE),''))}
  if(y=='rejection')abline(h=.05,lty=2,col='gray70');for(m in names(cols)){b<-a[a$method==m,];b<-b[order(b$size),];ec<-adjustcolor(cols[m],alpha.f=.55)
   interval(b$q25,b[[y]],b$q75,b[[y]],ec)
   if(y=='relative_rmse')interval(b$seconds,pmax(b[[y]]-1.96*b$rmse_mcse,.Machine$double.eps),b$seconds,b[[y]]+1.96*b$rmse_mcse,ec)
   else {zz<-qnorm(.975);den<-1+zz^2/b$reps;ctr<-(b[[y]]+zz^2/(2*b$reps))/den;half<-zz*sqrt(b[[y]]*(1-b[[y]])/b$reps+zz^2/(4*b$reps^2))/den;interval(b$seconds,ctr-half,b$seconds,ctr+half,ec)}
   lines(b$seconds,b[[y]],type='b',col=cols[m],pch=pchs[m],lty=lty[m],lwd=1.6,cex=.9);if(m!='ordinary'){ii<-unique(c(1,nrow(b)));label_pos<-switch(m,mn=c(1,3),cbb=c(3,1),blb=c(2,2));if(y=='rejection')label_pos[b[[y]][ii]<.12]<-3;text(b$seconds[ii],b[[y]][ii],labels=b$size[ii],pos=label_pos[seq_along(ii)],cex=.95,col=cols[m])}
   if(y=='rejection'){nul<-z[z$case==ca&z$delta==0&z$method==m,];nul<-nul[order(nul$size),];lines(nul$seconds,nul$rejection,col=cols[m],lty=4,lwd=.8);points(nul$seconds,nul$rejection,col=cols[m],pch=pchs[m],cex=.55)}}
  if(y=='rejection')mtext('Thin dot-dashed curves: null rejection',side=3,line=.15,cex=.95)}
 par(mar=c(0,0,0,0));plot.new();legend('center',legend=c('Ordinary','m/n (vary m)','CBB (vary s)','BLB (vary b; K = ceiling(n/b))'),col=cols,pch=pchs,lty=lty,ncol=2,bty='n',cex=.95,y.intersp=1.2);dev.off()}
draw('final_tuning_rmse.pdf',0,'relative_rmse','Relative RMSE of\nstandard-error estimate')
draw('final_tuning_power.pdf',2,'rejection','Rejection probability')
