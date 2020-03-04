volcanoexp = function(table,plotname,fileprefix,pthreshold,fc_varname,fcp,fcn,pval_varname,xlims,ylims,minfold=0) {
  
  tablj     = na.omit(table)
  if (minfold == 0) {
    tabsp     = tablj[tablj[,fc_varname] > 0 & tablj[,pval_varname] < pthreshold,]
    tabsn     = tablj[tablj[,fc_varname] < 0 & tablj[,pval_varname] < pthreshold,]
  } else {
    tabsp     = tablj[tablj[,fc_varname] > minfold  & tablj[,pval_varname] < pthreshold,]
    tabsn     = tablj[tablj[,fc_varname] < -minfold & tablj[,pval_varname] < pthreshold,]
  }
  signifpos = nrow(tabsp)
  signifneg = nrow(tabsn)
  signiftot = sum(signifpos+signifneg)
  
  # limits (if "self", evaluate each time)
  if (xlims == "self") {
    xlims = ceiling(max(abs(range(tablj[,fc_varname]))))
    xlims = c(-xlims,xlims)
  }
  if (ylims == "self") {
    ylims = ceiling(max(abs(range(log(tablj[,pval_varname],10)))))
    ylims = c(0,-ylims)
  }
  
  pdf(file=paste(fileprefix,".",plotname,"_volcano.pdf",sep=""),height=4.5,width=5)
  
  # plot volcano
  plot(x=tablj[,fc_varname],
       y=log(tablj[,pval_varname],10),
       col=alpha("slategray",0.6),
       main=plotname,cex=0.5,
       xlab=fc_varname,ylab="log10(pval)",ylim=ylims,xlim=xlims,
       sub=paste("total p<",pthreshold," N=",signiftot,"/",nrow(tablj),
                 " | pos:",fcp," neg:",fcn,sep=""))
  text(paste("dif>",minfold," & p<",pthreshold,": N=",signifpos,sep=""),y=ylims[2],x=xlims[2]/2,col="red",cex=0.8)
  text(paste("dif<-",minfold," & p<",pthreshold,": N=",signifneg,sep=""),y=ylims[2],x=xlims[1]/2,col="red",cex=0.8)
  if (minfold == 0) {
    abline(v=0,lty=2)
  } else {
    abline(v=c(-minfold,0,minfold),lty=2, col=c("blue","black","blue"))
  }
  abline(h=log(pthreshold,10),lty=2,col="blue")
  
  # plot distribution of FC values (diagnostic)
  ecdfi = ecdf(abs(tablj[,fc_varname]))
  plot(x=seq(0,xlims[2],length.out=1000), y=ecdfi(v = seq(0,xlims[2],length.out=1000)),
       col="slategray",type="l",
       xlab=fc_varname,main=paste("CDF abs",fc_varname),
       ylab="Fraction",
       xlim = c(0,xlims[2]),ylim=c(0,1))
  if (minfold == 0) {
    abline(v=0,lty=2)
  } else {
    abline(v=c(-minfold,0,minfold),lty=2, col=c("blue","black","blue"))
  }
  abline(h=0.5,lty=2)
  
  # plot distribution of FC values (below and above pval threshold)
  ecdfi = ecdf(abs(tablj[tablj[,pval_varname]<pthreshold,fc_varname]))
  plot(x=seq(0,xlims[2],length.out=1000), y=ecdfi(v = seq(0,xlims[2],length.out=1000)),
       col="blue",type="l",
       xlab=fc_varname,main=paste("CDF abs",fc_varname,"signif and non-signif"),
       ylab="Fraction",
       xlim = c(0,xlims[2]),ylim=c(0,1))
  ecdfi = ecdf(abs(tablj[tablj[,pval_varname]>pthreshold,fc_varname]))
  lines(x=seq(0,xlims[2],length.out=1000), y=ecdfi(v = seq(0,xlims[2],length.out=1000)),col="slategray3")
  legend("bottomright",legend = c("p<threshold","p>threshold"),
         col=c("blue","slategray3"),lty = 1,bty = "n",cex=0.8)
  if (minfold == 0) {
    abline(v=0,lty=2)
  } else {
    abline(v=c(-minfold,0,minfold),lty=2, col=c("blue","black","blue"))
  }
  abline(h=0.5,lty=2)
  
  
  # plot distribution of pvalues (diagnostic)
  ecdfi = ecdf(log(tablj[,pval_varname],10))
  plot(x=seq(ylims[2],ylims[1],length.out=1000), y=ecdfi(v = seq(ylims[2],ylims[1],length.out=1000)),
       col="slategray",type="l",
       xlab="log10(p)",main=paste("CDF",pval_varname),
       ylab="Fraction",
       xlim = ylims,ylim=c(0,1))
  abline(v=log(pthreshold,10),lty=2,col="blue")
  abline(h=0.5,lty=2)
  
  # plot distribution of p values (below and above FC threshold)
  ecdfi = ecdf(log(tablj[abs(tablj[,fc_varname])>minfold,pval_varname],10))
  plot(x=seq(ylims[2],ylims[1],length.out=1000), y=ecdfi(v = seq(ylims[2],ylims[1],length.out=1000)),
       col="blue",type="l",
       xlab="log10(p)",main=paste("CDF",pval_varname,"signif and non-signif"),
       ylab="Fraction",
       xlim = ylims,ylim=c(0,1))
  if (minfold == 0) {
    abline(v=0,lty=2)
  } else {
    ecdfi = ecdf(log(tablj[abs(tablj[,fc_varname])<minfold,pval_varname],10))
    lines(x=seq(ylims[2],ylims[1],length.out=1000), y=ecdfi(v = seq(ylims[2],ylims[1],length.out=1000)),col="slategray3")
  }
  legend("topright",legend = c("dif>|threshold|","dif<|threshold|"),
         col=c("blue","slategray3"),lty = 1,bty = "n",cex=0.8)
  abline(v=log(pthreshold,10),lty=2,col="blue")
  abline(h=0.5,lty=2)
  
  
  dev.off()
  
  write.table(tabsp,file=paste(fileprefix,".",plotname,"_difpos.csv",sep=""),
              quote = F)
  write.table(tabsn,file=paste(fileprefix,".",plotname,"_difneg.csv",sep=""),
              quote = F)
  
  signiflists = list(genes_sp=rownames(tabsp),
                     genes_sn=rownames(tabsn))
  return(signiflists)
  
}
