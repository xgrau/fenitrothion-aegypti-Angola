dimredfun = function(
  matriu,outputname,varname,cols_are,rows_are,isbidi,cols_dist_method,
  rows_dist_method,clus_method,
  printpdfpca=T,printpdfheatmaps=T,hm_height=8,hm_width=8) {
  mi = matriu
  
  library(gplots)
  library(ape)
  library(factoextra)
  
  # Colors
  if (isbidi == T) {
    hop = colorRampPalette(c("firebrick4","orangered",
                             "floralwhite",
                             "deepskyblue","dodgerblue4"))
    div = colorRampPalette(c("firebrick4","orangered",
                             "floralwhite",
                             "deepskyblue","dodgerblue4"))
    
  } else {
    hop = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))
    div = colorRampPalette(c("firebrick4","orangered",
                             "floralwhite",
                             "deepskyblue","dodgerblue4"))
  }
  
  # PCA
  message(paste("# PCA"))
  pic                  = prcomp(t(mi))
  pic$variancefraction = pic$sdev^2/sum(pic$sdev^2)
  
  # Correlation, distance & clustering
  message(paste("# Cols",cols_dist_method,clus_method))
  if (cols_dist_method == "spearman" | cols_dist_method == "pearson" | cols_dist_method == "kendall") {
    mic.cor  = cor(mi,method = cols_dist_method)
    mic.dist = as.dist(1-mic.cor)
    mic.clus = hclust(mic.dist,method=clus_method)
  } else {
    mic.cor  = data.frame()
    mic.dist = dist(t(mi),method=cols_dist_method)
    mic.clus = hclust(mic.dist,method=clus_method)
  }
  
  message(paste("# Rows",rows_dist_method,clus_method))
  if (rows_dist_method == "spearman" | rows_dist_method == "pearson" | rows_dist_method == "kendall") {
    mir.cor  = cor(t(mi),method=rows_dist_method)
    mir.dist = as.dist(1-mir.cor)
    mir.clus = hclust(mir.dist,method=clus_method)
  } else {
    mir.cor  = data.frame()
    mir.dist = dist(mi,method=rows_dist_method)
    mir.clus = hclust(mir.dist,method=clus_method)
  }
  
  # PCoA
  message(paste("# PCoA"))
  pio      = pcoa(mic.dist)
  pio12    = as.data.frame(pio$vectors[,c(1,2,3)])
  pio12$sp = rownames(pio12)
  
  # Plots PCA/PCOA
  if (printpdfpca) {
    
    message(paste("# Plot PCA & PCoA"))
    pdf(paste(outputname,".reddim.pdf",sep=""),height=9,width=8)
    par(mfrow=c(2,2))
    
    # Plot PCA: 12
    plot(pic$x[,c(1,2)],col="red",main="PCA 1&2",
         xlab=paste("PCA 1",signif(pic$variancefraction[1]*100,digits=3),"% variance"),
         ylab=paste("PCA 2",signif(pic$variancefraction[2]*100,digits=3),"% variance"))
    text(pic$x[,c(1,2)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCA: 13
    plot(pic$x[,c(1,3)],col="red",main="PCA 1&3",
         xlab=paste("PCA 1",signif(pic$variancefraction[1]*100,digits=3),"% variance"),
         ylab=paste("PCA 3",signif(pic$variancefraction[3]*100,digits=3),"% variance"))
    text(pic$x[,c(1,3)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCA: 23
    plot(pic$x[,c(2,3)],col="red",main="PCA 2&3",
         xlab=paste("PCA 2",signif(pic$variancefraction[2]*100,digits=3),"% variance"),
         ylab=paste("PCA 3",signif(pic$variancefraction[3]*100,digits=3),"% variance"))
    text(pic$x[,c(2,3)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCA: eigenvalues
    plot(pic$variancefraction,type="b",col="red",main="PCA variance")
    text(pic$variancefraction,labels=signif(pic$variancefrac,digits=3),pos=3,col="black")
    
    # Plot PCoA: 12
    plot(x=pio12$Axis.1,y=pio12$Axis.2,col="blue",main="PCoA 1&2",
         xlab=paste("Axis 1",signif(pio$values$Relative_eig[1]*100,digits=3),"% variance"),
         ylab=paste("Axis 2",signif(pio$values$Relative_eig[2]*100,digits=3),"% variance"))
    text(pio$vectors[,c(1,2)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCoA: 13
    plot(x=pio12$Axis.1,y=pio12$Axis.3,col="blue",main="PCoA 1&3",
         xlab=paste("Axis 1",signif(pio$values$Relative_eig[1]*100,digits=3),"% variance"),
         ylab=paste("Axis 3",signif(pio$values$Relative_eig[3]*100,digits=3),"% variance"))
    text(pio$vectors[,c(1,3)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCoA: 23
    plot(x=pio12$Axis.2,y=pio12$Axis.3,col="blue",main="PCoA 2&3",
         xlab=paste("Axis 2",signif(pio$values$Relative_eig[2]*100,digits=3),"% variance"),
         ylab=paste("Axis 3",signif(pio$values$Relative_eig[3]*100,digits=3),"% variance"))
    text(pio$vectors[,c(2,3)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCoA: eigenvalues
    plot(pio$values$Relative_eig,col="blue",main="PCoA relative eigenvalues",type="b",xlab="Axis rank",ylab="Fraction variance explained")
    text(pio$values$Relative_eig,labels=signif(pio$values$Relative_eig,digits=3),pos=3,col="black")
    lines(pio$values$Broken_stick[pio$values$Broken_stick>.01],
          col="gray",main="PCoA variance broken stick",type="b")
    
    # Plot clustering of columns
    plot(mic.clus, main = paste(cols_are," (",ncol(mi),")",sep=""), sub = paste("clustering:",clus_method,"distances:",cols_dist_method))
    
    dev.off()
  }
  
  # plot heatmaps
  message(paste("# Plot heatmaps"))
  if (printpdfheatmaps){
    
    pdf(paste(outputname,".heatmaps.pdf",sep=""),height=hm_height,width=hm_width)
    par(mfrow=c(1,1))
    
    labrowbool = if(nrow(mi) > 500) { F } else { rownames(mi) }
    labcolbool = if(ncol(mi) > 500) { F } else { colnames(mi) }
    
    # heatmap samples vs variables
    heatmap.2(mi,
              col =hop(31),dendrogram = "both",symkey=F,trace="none",scale="none",
              useRaster=F,keysize = 1,main = paste(cols_are," (",ncol(mi),")"," ~ ",
                                                   rows_are," (",nrow(mi),")",sep=""),symm=F,
              labRow = labrowbool,
              labCol = labcolbool,
              Colv=as.dendrogram(mic.clus),
              Rowv=as.dendrogram(mir.clus),
              margins=c(10,10),na.color = "grey90",key.title = varname)
    
    # heatmap variables correlation
    if ((rows_dist_method == "spearman" | rows_dist_method == "pearson" | rows_dist_method == "kendall") & nrow(mir.cor) < 1000 ) {
      
      labrowbool = if(nrow(mir.cor) > 500) { F } else { rownames(mir.cor) }
      labcolbool = if(ncol(mir.cor) > 500) { F } else { colnames(mir.cor) }
      
      heatmap.2(mir.cor,
                col =div(31),dendrogram = "both",symbreaks=T,trace="none",scale="none",
                useRaster=F,keysize = 1,main = paste("Pairwise",rows_are),symm=F,
                labRow = labrowbool,
                labCol = labcolbool,
                Colv=as.dendrogram(mir.clus),
                Rowv=as.dendrogram(mir.clus),
                margins=c(10,10),na.color = "grey90",key.title = "Corr. coef.")
    }
    
    # heatmap samples correlation
    if ((cols_dist_method == "spearman" | cols_dist_method == "pearson" | cols_dist_method == "kendall") & nrow(mic.cor) < 1000 ) {
      
      labrowbool = if(nrow(mic.cor) > 500) { F } else { rownames(mic.cor) }
      labcolbool = if(ncol(mic.cor) > 500) { F } else { colnames(mic.cor) }
      
      heatmap.2(mic.cor,
                col =div(31),
                dendrogram = "both",symbreaks=T,trace="none",scale="none",
                useRaster=F,keysize = 1,main = paste("Pairwise",cols_are),symm=F,
                Colv=as.dendrogram(mic.clus),
                labRow = labrowbool,
                labCol = labcolbool,
                Rowv=as.dendrogram(mic.clus),
                margins=c(10,10),na.color = "grey90",key.title = "Corr. coef.")
    }
    
    dev.off()
  }
  
  # output: clusterings etc.
  return.list = list(
    "pio"      = pio,
    "pic"      = pic,
    "mic.clus" = if ((cols_dist_method == "spearman" | cols_dist_method == "pearson" | cols_dist_method == "kendall") & nrow(mic.cor) < 1000 ) { mic.clus } else { NA },
    "mic.dist" = if ((cols_dist_method == "spearman" | cols_dist_method == "pearson" | cols_dist_method == "kendall") & nrow(mic.cor) < 1000 ) { mic.dist } else { NA },
    "mir.clus" = if ((rows_dist_method == "spearman" | rows_dist_method == "pearson" | rows_dist_method == "kendall") & nrow(mir.cor) < 1000 ) { mir.clus } else { NA },
    "mir.dist" = if ((rows_dist_method == "spearman" | rows_dist_method == "pearson" | rows_dist_method == "kendall") & nrow(mir.cor) < 1000 ) { mir.dist } else { NA }
  )
  
  return(return.list)
  
}


kmeansfun = function(
  matrix, outputprefix, k_list=c(3,4,5), ylab="counts", 
  iter.max=10, nstart=1, algorithm="Hartigan", returnfile = T,
  alpha=0.2,ylims=c(0,1)) {
  
  library(cluster)
  library(gridExtra)
  library(pheatmap)
  
  mat = matrix
  mat = mat[apply(mat,1,sd) != 0,]
  mat = as.matrix(mat)
  
  # k-means clustering
  mat.kmeans = list() # kmeans clustering data
  mat.ksilho = list() # mean silhouette score for each value of k
  mat.clusters.ann = data.frame(row.names = rownames(mat))
  for (k in k_list) {
    kc=paste("k",k,sep="")
    message("# k-means clust ", kc)
    mat.kmeans[k] = list(kmeans(mat,centers = k,iter.max=iter.max, nstart = nstart, algorithm = algorithm))
    mat.clusters.ann[paste("k",k,sep="")] = as.factor(mat.kmeans[[k]]$cluster)
    mat.ksilho[k] = mean(cluster::silhouette(as.vector(mat.kmeans[[k]]$cluster),dist = dist(mat))[,3])
  }
  
  # plots & table (optional)
  if (returnfile) {
    # plot all profiles according to cluster
    pdf(file=paste(outputprefix,"kmeans.pdf",sep="."),height=2.2*max(k_list),width=8)
    for (k in k_list) {
      par(mfrow=c(max(k_list),2))
      kc=paste("k",k,sep="")
      message("# k-means plots ", kc)
      for (clui in levels(mat.clusters.ann[,kc])) {
        mat.clui = mat[mat.clusters.ann[,kc] == clui,]
        # plot lines with profile
        par(mar=c(5.1,4.1,4.1,2.1))
        plot(x=1:ncol(mat.clui),y=rep(1,ncol(mat.clui)),xlim=c(1,ncol(mat.clui)),
             main=paste("cluster ",kc," k=",clui,", N=",nrow(mat.clui),", silhouette=",signif(mat.ksilho[[k]],3),sep=""),
             pch=NA,xaxt="n",xlab=NA,ylab=ylab,las=1,ylim=ylims)
        axis(1,at=1:ncol(mat.clui),labels=colnames(mat.clui),las=2)
        for (evi in 1:nrow(mat.clui)) {
          lines(x=1:ncol(mat.clui),y=mat.clui[evi,],col=alpha("blue",alpha))
        }
        # plot boxplots with profile
        par(mar=c(5.1,4.1,4.1,2.1))
        boxplot(mat.clui,las=2,col="slategray3",
                pch=NA,ylab=ylab,las=1,ylim=ylims)
      }
    }
    
    # plot clustering groupings at different k
    pheatmap::pheatmap(
      t(mat),cellheight = 0,treeheight_col=0,
      cluster_rows=F,cluster_cols = T,show_colnames = F,show_rownames = F,
      annotation_col = mat.clusters.ann,legend=F)
    
    # plot clustering grouping similarity: Rand index
    plot(0,pch=NA,xlab=NA,ylab=NA,axes=F)
    mat.clusters.rand = data.frame()
    for (pi in 1:ncol(mat.clusters.ann)) {
      for (pj in 1:ncol(mat.clusters.ann)) {
        mat.clusters.rand.i = data.frame(colnames(mat.clusters.ann)[pi],
                                         colnames(mat.clusters.ann)[pj],
                                         mclust::adjustedRandIndex(mat.clusters.ann[,pi],mat.clusters.ann[,pj]))
        mat.clusters.rand = rbind(mat.clusters.rand, mat.clusters.rand.i)
      }
    }
    colnames(mat.clusters.rand) = c("ki","kj","Rand index")
    grid.table(mat.clusters.rand)
    
    dev.off()
    
    # save table with kmeans clustering
    mat.clusters.table = cbind(mat.clusters.ann,mat)
    write.table(mat.clusters.table,
                file=paste(outputprefix,"kmeans.csv",sep="."),
                row.names = T,col.names = T,quote=F,sep = "\t")
    
  }
  
  # return
  output = list("table" = mat.clusters.table, "clusterings" = mat.kmeans)
  return(output)
  
  
}
