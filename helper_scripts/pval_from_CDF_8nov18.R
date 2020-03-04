## defines two functions: pval calculation and friendly volcano plots
# 
# matrix_frequencies = taf
# matrix_tpms = tac
# name_group_i = coi
# name_group_j = coj
# samples_group_i = cli
# samples_group_j = clj
# pval_correction = "fdr"
# ecdf_area = 1000

pval_from_CDF = function(matrix_frequencies,matrix_tpms,
                         name_group_i, name_group_j,
                         samples_group_i, samples_group_j,
                         comparison_code,
                         pval_correction="BH",ecdf_area=1000,
                         dpsi_cutoff=0) {
  
  # rename variables
  taf = matrix_frequencies
  tac = matrix_tpms
  coi = name_group_i
  coj = name_group_j
  cli = samples_group_i
  clj = samples_group_j
  con = comparison_code
  arv = ecdf_area
  
  # find pval from ECDF
  pval_from_CDFvicinity = function(vector,start,final,val_eval) {
    
    if (is.na(val_eval)) {
      pval = NA
    } else {
      if (abs(val_eval) > dpsi_cutoff) {
        vecti  = vector[start:final]
        vecti  = vecti[!is.na(vecti)]
        vecti  = c(0,vecti,1)
        f.ecdf = ecdf(vecti)
        # f.ecdf(val_eval)
        pval = (1-f.ecdf(val_eval)) / 2
      } else {
        pval = 1
      }
    }
    return(pval)
  
  }
  
  # find closer value in vector
  closestVal.fun = function(vector,value) {
    which.min(abs(vector - value)) 
  }
  
  # across-condition absolute frequency difference, and TPMs
  dis.ac = data.frame(
    event  = rownames(taf),
    freq_diff  = rowMeans(taf[,cli],na.rm = T)        - rowMeans(taf[,clj],na.rm = T),
    freq_diffA = abs(rowMeans(taf[,cli],na.rm = T)    - rowMeans(taf[,clj],na.rm = T)),
    mean_expr  = 1/2 * (rowMeans(log10(tac[,cli]),na.rm = T) + rowMeans(log10(tac[,clj]),na.rm = T))
  )
  
  # inter-replicate distributions of diff_edifreq and mean_expression
  # within groups i and j (separately)
  dis.ir = data.frame()
  for (clx in list(cli,clj)) {
    for (clx.i in 1:length(clx)) {
      for (clx.j in 1:length(clx)) {
        if (clx.i < clx.j) {
          dis.ir = rbind(dis.ir,data.frame(
            event      = rownames(taf),
            freq_diff  = taf[,clx[clx.i]] - taf[,clx[clx.j]],        # difference in frequency for comparison i-j
            freq_diffA = abs(taf[,clx[clx.i]] - taf[,clx[clx.j]]),   # absolute difference (for CDF)
            mean_expr  = 1/2 * (log10(tac[,clx[clx.i]]) + log10(tac[,clx[clx.j]])) # mean TPM count in the event
          )
          )
        }
      }
    }
  }
  
  dis.ir   = na.omit(dis.ir)
  # dis.ir.a = aggregate(dis.ir,by=list(dis.ir$event),FUN=mean) # mean aggregate
  dis.ir.a = dis.ir
  # dis.ir.a = dis.ir.a[dis.ir.a$freq_diffA!= 0, ]
  dis.ir.a = dis.ir.a[base::order(-dis.ir.a$mean_expr),]
  
  
  # find position of expression-ordered replicates that
  # is closest to expression of query (and its vicinity)
  dis.ac$ref.closest = sapply(dis.ac$mean_expr,   function(x) closestVal.fun(dis.ir.a$mean_expr,x))
  dis.ac$ref.start   = sapply(dis.ac$ref.closest, function(x) max(x-(arv/2),0))
  dis.ac$ref.final   = sapply(dis.ac$ref.closest, function(x) min(x+(arv/2),nrow(dis.ir.a)))
  
  # calculate pvalue from CDF of within-replicate editing
  # frequency distribution, using only frequencies with 
  # corresponding expression in the vicinity of 
  dis.ac$freq_diffA[is.nan(dis.ac$freq_diffA)] = NA
  dis.ac$freq_diff [is.nan(dis.ac$freq_diff)]  = NA
  dis.ac$pval = apply(dis.ac[,c("freq_diffA","ref.start","ref.final")], 1, function(x)
    pval_from_CDFvicinity(vector   = dis.ir.a$freq_diffA,
                          start    = x[2], 
                          final    = x[3],
                          val_eval = x[1])
  )
  
  # prepare output
  dis.ac             = subset(dis.ac, select=c("event","freq_diff","mean_expr","pval"))
  dis.ac$pval.corr   = p.adjust(dis.ac$pval, method = pval_correction)
  dis.ac$comparison  = con
  dis.ac$group.i     = coi
  dis.ac$group.j     = coj
  dis.ac$mean_freq.i = rowMeans(taf[,cli],na.rm = T)
  dis.ac$mean_freq.j = rowMeans(taf[,clj],na.rm = T)

  plotfun = function() {
    plot(y=dis.ac$freq_diff, x=dis.ac$mean_expr,
         col=c("slategray3",rainbow(5,v = 0.7)[2:4])[as.factor(as.numeric(dis.ac$pval<0.1) + 
                                                                 as.numeric(dis.ac$pval<0.05) + 
                                                                 as.numeric(dis.ac$pval<0.01))],
         xlab="expression",ylab="freq_diff",cex=0.6,
         main=paste("freq_diff vs expr & signif"))
    legend("topright",legend = c(paste("p<0.1 n=",sum(as.numeric(dis.ac$pval<0.1),na.rm = T),sep=""),
                                 paste("p<0.05 n=",sum(as.numeric(dis.ac$pval<0.05),na.rm = T),sep=""),
                                 paste("p<0.01 n=",sum(as.numeric(dis.ac$pval<0.01),na.rm = T),sep=""))
           ,col=rainbow(5,v = 0.7)[2:4],pch=1,cex=0.6)
    abline(h=0,lty=2)
  }
    
  return(list(result = dis.ac, plot = plotfun, intrareplicate = dis.ir.a))
  
}


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
       xlim = c(0,xlims[2]))
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
