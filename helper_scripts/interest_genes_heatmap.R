interest_genes_heatmap = function(interestlist, interestname, source_matrix) {
  
  original_len = length(interestlist)
  interestmatr = subset(source_matrix, rownames(source_matrix) %in% interestlist)
  interestlist = rownames(interestmatr)
  
  interestmatr.ann = data.frame(row.names = interestlist)
  interestmatr.exp = matrix(nrow = 0,ncol = 7)
  
  ann_colors = list()
  
  for (com in complis) {
    tici     = paste(com[5],"_",com[2],"-",com[3],sep="")
    load(paste(outcode,"de_co.",tici,".data.RData",sep=""))
    isposi = (interestlist %in% ex_lis_i$genes_sp) * 1
    isnegi = (interestlist %in% ex_lis_i$genes_sn) * -1
    interestmatr.ann[tici] = as.factor(isposi + isnegi)
    interestmatr.exi       = ex_res_i
    interestmatr.exi$comparison = tici
    interestmatr.exi$gene       = rownames(interestmatr.exi)
    interestmatr.exi$samples_i     = com[2]
    interestmatr.exi$samples_j     = com[3]
    interestmatr.exi               = merge(interestmatr.exi,dign,by.x="gene",by.y="gene_id")
    interestmatr.exi$is_signif_pos = interestmatr.exi$gene %in% ex_lis_i$genes_sp
    interestmatr.exi$is_signif_neg = interestmatr.exi$gene %in% ex_lis_i$genes_sn
    interestmatr.exi$is_signif_any = interestmatr.exi$is_signif_pos | interestmatr.exi$is_signif_neg
    interestmatr.exi       = interestmatr.exi[interestmatr.exi$gene %in% interestlist,]
    interestmatr.exp       = rbind(interestmatr.exp, interestmatr.exi)
    ann_colors[[tici]]     = c("-1" = "magenta3","0"="slategray1","1"="green3")
  }
  
  write.table(interestmatr.exp,file=paste(outcode,"subset.",interestname,".csv",sep=""),
              quote = F,row.names = F, sep="\t")
  
  pdf(file=paste(outcode,"subset.",interestname,".pdf",sep=""),height=10+length(interestlist)/4,width=10)
  rownames(interestmatr) = gsub("Anogam_","",rownames(interestmatr))
  rownames(interestmatr.ann) = gsub("Anogam_","",rownames(interestmatr.ann))
  pheatmap(interestmatr, color = col.fun(21), breaks = seq(0,4,length.out = 20),
           border_color = "white", annotation_row = interestmatr.ann,cellheight = 6, cellwidth = 6, cluster_cols=F, fontsize = 7,
           distance_rows = "correlation",scale="none",annotation_colors = ann_colors,
           gaps_col = seq(min(table(as.numeric(sf$group))),to=length(sf$group), by=min(table(as.numeric(sf$group)))),
           main=paste("DE: ",interestname," (expressed: ",length(interestlist)," out of ",original_len,")",sep=""))
  dev.off()
  
}



interest_genes_heatmap2 = function(interestlist, interestname, source_matrix, gaps_col=0) {
  
  original_len = length(interestlist)
  interestmatr = subset(source_matrix, rownames(source_matrix) %in% interestlist)
  interestlist = rownames(interestmatr)
  
  interestmatr.ann = data.frame(row.names = interestlist)
  interestmatr.exp = matrix(nrow = 0,ncol = 7)
  
  write.table(interestmatr.exp,file=paste(outcode,"subset.",interestname,".csv",sep=""),
              quote = F,row.names = F, sep="\t")
  
  pdf(file=paste(outcode,"subset.",interestname,".pdf",sep=""),height=10+length(interestlist)/4,width=10)
  rownames(interestmatr) = gsub("Anogam_","",rownames(interestmatr))
  rownames(interestmatr.ann) = gsub("Anogam_","",rownames(interestmatr.ann))
  pheatmap(interestmatr, color = col.fun(21), breaks = seq(0,4,length.out = 20),
           border_color = "white", 
           cellheight = 6, cellwidth = 6, cluster_cols=F, fontsize = 7,
           distance_rows = "correlation",scale="none",
           gaps_col = gaps_col,
           main=paste("DE: ",interestname," (expressed: ",length(interestlist)," out of ",original_len,")",sep=""))
  dev.off()
  
}
