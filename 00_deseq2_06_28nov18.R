#### Input vars ####

setwd("~/Dropbox/Transvar/deseq_Aedaeg_Angola/")

si      = "Aedaeg"
outcode = "all"
gfffile = "~/dades/Genomes/Aedaeg_long.annot.gff"
genfile = "~/dades/Genomes/Aedaeg_gDNA.fasta"
txgdict = "~/dades/Genomes/Aedaeg_long.cds.csv"
nmgdict = "~/dades/Genomes/Aedaeg_genes.description.csv"
gomapfi = "~/dades/Genomes/anotacions/Aedaeg_long.pep_eggnog_diamond.emapper.annotations.GO"
panfile = "~/dades/Genomes/anotacions/Aedaeg_long.pep_Pfamscan.seqs"
panform = "pfamscan"
saloutf = "~/dades/Transcriptomes/Aedaeg_Angola_19oct18/salmon/"
samgrup = "samples_expanded.class" # can contain >1 grouping (eg location+treatment)
samclas = c("group","phenotype","exposure","location","groupsim")
expdiss = "~group" # global formula DEseq (comparison-specific formulae indicated below)

# define comparisons:
# 1   -> grouping within which comparisons are made
# 2&3 -> groups 1 & 2
# 4   -> extra grouping, to restrict comparisons within it
#        ALL if no restrictions are intended
# 5   -> group within the extra grouping that will be used;
#        ALL if no restrictions are intended
complis = list(
  c("group","RUN","SNO","ALL","ALL","~group"),
  c("group","RUN","SRO","ALL","ALL","~group"),
  c("group","RUN","RFE","ALL","ALL","~group")
)

# pairs of comparisons (use index from complis) that will be
# compared in the Venn diagrams to see overlap in DE results
# BEWARE: MUST USE SAME CONVENTION FOR OVER/UNDEREXPRESSION!
# Same group of samples must be "first" in both comparisons.
# If no comparison is to be done, just sat c(1,1)
complis_pairs = list(
  complis[c(1,2)]
)

# input variables
pval_threshold  = 0.001
alpha_threshold = 0.001
fc_threshold    = 2      # FC threshold used to define over/under expression
                         # FC = 1  -> log2FC = 0
                         # FC = 2  -> log2FC = 1
                         # FC = 10 -> log2FC = 3.27...
logfc_threshold = log2(fc_threshold)


# load libraries
library(ape)
library(gplots)
library(stringr)
library(plyr)
library(tidyr)
library(topGO)
library(tximport)
library(readr)
library(DESeq2)
library(pheatmap)
library(cluster)
library(VennDiagram)

# load functions
source("~/Dropbox/Scripts/bin/geneSetAnalysis.R")
source("~/Dropbox/Scripts/bin/dimred_18jul18.R")
source("~/Dropbox/Scripts/bin/pval_from_CDF_8nov18.R")
source("~/Dropbox/Scripts/bin/slidingfunctions_18nov18.R")
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))
graphics.off()


#### Input files ####

# sample classification
sf           = read.table(samgrup,header = F)
colnames(sf) = c("sample",samclas)
rownames(sf) = sf$sample

# transcript to gene dictionary
di           = read.table(txgdict)
colnames(di) = c("long_transcript_id","gene_id","transcript_length")
di           = di[,1:2]

# gene name dictionary
dign           = read.table(nmgdict,sep = "\t")
colnames(dign) = c("gene_id","gene_name")


# functional mappings: GO, pfam
gomap = readMappings(gomapfi)
panno = read.table(file = panfile)
if (panform == "simple") {
  colnames(panno) = c("transcript","pstart","pend","domain")
} else if (panform == "pfamscan") {
  colnames(panno) = c("transcript","pstart","pend","pfamid","domain","domseq")
}
panno = merge(panno,di, by.x="transcript",by.y="long_transcript_id")
pannu = subset(panno, select=c("transcript","pfamid","domain","gene_id"))
pannu = pannu[!duplicated(pannu), ]


# load gff
message("# Loading input: GFF")
exon_name     = "CDS"
gene_name     = "gene"
mRNA_name     = "mRNA"
gf            = read.gff(gfffile)
gf$attributes = gsub("ID=","",gf$attributes)
gf$attributes = gsub("Parent=","",gf$attributes)
gf$attributes = gsub(";.*","",gf$attributes)
gf            = subset(gf, type %in% c(gene_name))
levels(gf$type)[levels(gf$type)==exon_name] = "exon"
levels(gf$type)[levels(gf$type)==mRNA_name] = "gene"
gf            = gf[order(gf[,"seqid"],gf[,"start"]),]
gf$source     = "R"
gi            = GRanges(
  gf$seqid,
  IRanges(start=gf$start, end=gf$end , names=gf$attributes),
  strand = gf$strand)



#### DESEQ all-to-all ####

# deseq input data
ex_lif        = file.path(saloutf,paste(si,"_",sf$sample,".salmon_long.out",sep=""), "quant.sf")
names(ex_lif) = sf$sample
ex_inp        = tximport(files = ex_lif, type = "salmon", tx2gene = di)
ex_dat        = DESeqDataSetFromTximport(txi=ex_inp,colData=sf,design=formula(expdiss))
ex_dat        = DESeq(ex_dat,test="Wald") # if tximport is used, a normMatrix with avgTxLength is used as normalization factor (which superseeds size factors)

# matrix of normalised counts 
# uses normalization factors: sample size, effective len...
ex_dnc = na.omit(DESeq2::counts(ex_dat,normalized=T))
ex_dnc = log(ex_dnc,10)
ex_dnc = abs(ex_dnc)
ex_dnc[!is.finite(ex_dnc)] = 0
ex_dnc = ex_dnc[apply(ex_dnc,1,sd) != 0,]

# dimensionality reduction
# prepare data: scale along rows
ex_dnc_st = ex_dnc
ex_dnc_st = t(apply(ex_dnc_st, 1, scale))
colnames(ex_dnc_st) = colnames(ex_dnc)

ex_dnc_dimred = dimredfun(
  matriu = ex_dnc_st, outputname = paste(outcode,"de_all",sep="."), varname = "standardized log norm counts", 
  cols_are = "Samples",rows_are = "Genes", isbidi = F,
  cols_dist_method = "pearson", rows_dist_method = "pearson",
  clus_method = "ward.D2")



#### K-means clustering ####

k_list = c(3,4,5) # list of k values to try for kmeans clustering

# kmeans
ex_dnc_kmeans = kmeansfun(
  matrix=ex_dnc_st,
  outputprefix = paste(outcode,"de_all",sep="."), 
  k_list = k_list, ylab = "standardized log norm counts", alpha=0.1, ylims = c(-3,3))



#### DEseq specific comparisons ####

for (com in complis) {
  
  # define columns with groups
  if (com[5] == "ALL") {
    cli = as.character(sf[sf[,com[1]] %in% com[2],]$sample)
    clj = as.character(sf[sf[,com[1]] %in% com[3],]$sample)
  } else {
    cli = as.character(sf[sf[,com[1]] %in% com[2] & sf[,com[4]] %in% com[5],]$sample)
    clj = as.character(sf[sf[,com[1]] %in% com[3] & sf[,com[4]] %in% com[5],]$sample)
  }
  tici   = paste(com[5],".",com[2],"-",com[3],sep="")
  
  # prepare data for DEseq
  message(paste("# DEseq comparison",tici))
  ex_lif_i = ex_lif[names(ex_lif) %in% c(cli,clj)]
  sf_i     = sf[sf$sample %in% c(cli,clj),]
  sf_i     = droplevels(sf_i)
  ex_inp_i = suppressMessages(tximport(files = ex_lif_i, type = "salmon", tx2gene = di))
  ex_dat_i = suppressMessages(DESeqDataSetFromTximport(txi = ex_inp_i, colData = sf_i, design = formula(com[6])))
  ex_dat_i[[com[1]]] = factor(ex_dat_i[[com[1]]], levels=c(com[3],com[2]))
  ex_dat_i = suppressMessages(DESeq(ex_dat_i, test="Wald"))
  
  # obtain effect size (LFCs), sval and pval
  ex_res_i               = results(ex_dat_i,alpha=alpha_threshold,contrast=com[1:3],name=tici) # unshrunken LFC with pval, FDR-adjusted at alpha
  ex_res_s               = lfcShrink(ex_dat_i,coef=2, type="apeglm", quiet=T, svalue = T)      # shrink LFC using apeglm
  ex_res_i               = as.data.frame(ex_res_i)
  ex_res_i$svalue_shrink = ex_res_s$svalue                                                     # s-value: probability of false signs, ie probability that the effect is actually LFC>0 or LFC<0 (as opposed to pval of probability that effect is 0 or not 0)
  ex_res_i$log2FC_shrink = ex_res_s$log2FoldChange                                             # shrunken LFC values
  
  # matrix of normalized counts
  ex_dnc_i = ex_dnc[,c(cli,clj)]
  ex_dnc_i = ex_dnc_i[apply(ex_dnc_i,1,sd) != 0,]
  ex_yes = order(rowVars(ex_dnc_i),decreasing=T)
  
  # volcano plot & list of over/under expressed genes
  message(paste("# Volcano",tici))
  ex_lis_i = volcanoexp(
    table=ex_res_i,plotname=tici,fileprefix=paste(outcode,"de_co",sep="."),
    pthreshold=alpha_threshold,fc_varname="log2FC_shrink",fcp=com[2],fcn=com[3],pval_varname="padj",
    xlims = c(-10,10),ylims = c(0,-100),minfold = logfc_threshold)
   
  # convert lists of genes to lists of transcripts
  ex_lis_i$tx_sp = as.vector(di[di$gene_id %in% ex_lis_i$genes_sp,]$long_transcript_id)
  ex_lis_i$tx_sn = as.vector(di[di$gene_id %in% ex_lis_i$genes_sn,]$long_transcript_id)
  
  # top genes
  message(paste("# Top genes",tici))
  pdf(file=paste(outcode,".de_co.",tici,"_topgenes.pdf",sep=""),height=12,width=8)
  par(mfrow=c(6,4))
  
  ex_res_i_sortp   = ex_res_i[ex_res_i$padj<pval_threshold,]
  ex_res_i_sortp.p = ex_res_i_sortp[order(ex_res_i_sortp$log2FC_shrink,decreasing = T),]
  ex_res_i_sortp.n = ex_res_i_sortp[order(ex_res_i_sortp$log2FC_shrink),]
  for (geni in 1:24) {
    geni_i = rownames(ex_res_i_sortp.p[geni,])
    geni.n = as.character(dign[dign$gene_id == geni_i,]$gene_name)
    DESeq2::plotCounts(
      ex_dat_i, gene=geni_i, intgroup=com[1], col="slategray",xlab="",
      cex_sub=0.7,cex_axis=0.7,cex_main=0.7,cex_lab=0.7,
      sub=paste("FC=",signif(2^(ex_res_i_sortp.p[geni,]$log2FC_shrink),4)
                ,"| p=",signif(ex_res_i_sortp.p[geni,]$padj,4),"\n",geni.n))
  }
  for (geni in 1:24) {
    geni_i = rownames(ex_res_i_sortp.n[geni,])
    geni.n = as.character(dign[dign$gene_id == geni_i,]$gene_name)
    DESeq2::plotCounts(
      ex_dat_i, gene=geni_i, intgroup=com[1], col="slategray",xlab="",
      cex_sub=0.7,cex_axis=0.7,cex_main=0.7,cex_lab=0.7,
      sub=paste("FC=",signif(2^-(ex_res_i_sortp.n[geni,]$log2FC_shrink),4)
                ,"(^-1) | p=",signif(ex_res_i_sortp.n[geni,]$padj,4),"\n",geni.n))
  }
  dev.off()
  
  # save comparison (can get reloaded later)
  save(list=c("ex_dat_i","ex_lis_i","ex_dnc_i","ex_res_i","sf_i"), 
       file=paste(outcode,"de_co",tici,"data.RData",sep="."))
  
}

for (com in complis) {
  
  # define columns with groups
  if (com[5] == "ALL") {
    cli = as.character(sf[sf[,com[1]] %in% com[2],]$sample)
    clj = as.character(sf[sf[,com[1]] %in% com[3],]$sample)
  } else {
    cli = as.character(sf[sf[,com[1]] %in% com[2] & sf[,com[4]] %in% com[5],]$sample)
    clj = as.character(sf[sf[,com[1]] %in% com[3] & sf[,com[4]] %in% com[5],]$sample)
  }
  
  tici   = paste(com[5],".",com[2],"-",com[3],sep="")
  load(file = paste(outcode,"de_co",tici,"data.RData",sep="."))
  
  # functional enrichment
  message(paste("# Functional enrichment domains",tici))
  hygeofun(list_interest=ex_lis_i$tx_sp,annotation=pannu,gene_col="transcript",ano_col="domain",
           outputname=paste(outcode,"de_co",sep="."),
           name_geneset=paste(tici,"difpos",sep="_"),topnum = 50)
  hygeofun(list_interest=ex_lis_i$tx_sn,annotation=pannu,gene_col="transcript",ano_col="domain",
           outputname=paste(outcode,"de_co",sep="."),
           name_geneset=paste(tici,"difneg",sep="_"),topnum = 50)
  
  message(paste("# Functional enrichment topgo pos",tici))
  suppressMessages(topgofun(list_interest=ex_lis_i$tx_sp,gomap=gomap,
           ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
           outputname=paste(outcode,"de_co",sep="."),
           name_geneset=paste(tici,"difpos",sep="_"),topnum=50))
  message(paste("# Functional enrichment topgo neg",tici))
  suppressMessages(topgofun(list_interest=ex_lis_i$tx_sn,gomap=gomap,
           ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
           outputname=paste(outcode,"de_co",sep="."),
           name_geneset=paste(tici,"difneg",sep="_"),topnum=50))
  
}



#### DE gene overlaps ####

ex_lis_joint = list()

pdf(file=paste(outcode,"de_all.gene_overlaps.pdf",sep="."),height=4,width=4)
for (cop in complis_pairs) {

  ex_lis_joini = list("pai","paj")
  # define comparison pair
  pai      = cop[[1]]
  paj      = cop[[2]]
  pai.tici = paste(pai[5],".",pai[2],"-",pai[3],sep="")
  paj.tici = paste(paj[5],".",paj[2],"-",paj[3],sep="")
  paji     = paste(pai.tici,paj.tici,sep=" ")
  
  # load expression data for comparison
  load(paste(outcode,"de_co",pai.tici,"data.RData",sep="."))
  ex_lis_joini[[1]]$genes_sp = ex_lis_i$genes_sp
  ex_lis_joini[[1]]$genes_sn = ex_lis_i$genes_sn
  load(paste(outcode,"de_co",paj.tici,"data.RData",sep="."))
  ex_lis_joini[[2]]$genes_sp = ex_lis_i$genes_sp
  ex_lis_joini[[2]]$genes_sn = ex_lis_i$genes_sn
  
  # set analysis for OVEREXPRESSED genes
  ex_lis_joint[[paji]]$overexpressed = 
    venn.two(
      list1=ex_lis_joini[[1]]$genes_sp,
      list2=ex_lis_joini[[2]]$genes_sp,
      catname1=pai.tici,
      catname2=paj.tici,
      main = paste("overexp: ",pai.tici,"&",paj.tici),
      col1 = "green3",
      col2 = "blue3"
    )
  
  # set analysis for OVEREXPRESSED genes
  ex_lis_joint[[paji]]$underexpressed = 
    venn.two(
      list1=ex_lis_joini[[1]]$genes_sn,
      list2=ex_lis_joini[[2]]$genes_sn,
      catname1=pai.tici,
      catname2=paj.tici,
      main = paste("underexp: ",pai.tici,"&",paj.tici),
      col1 = "green3",
      col2 = "blue3"
    )
  
}

dev.off()



#### Chromosomal expression bias ####

# input variables
k_len       = 100                     # window length (in # genes)
k_step      = 0.1                     # % of window length to use as fractional step (1=non-overlapping)
chromlist   = levels(gf$seqid)               # chromosome list
chromlist   = chromlist[table(gf$seqid)>100] # exclude chromosomes with less than 100 genes (mind k_len!)
log2fc_thr_chrp = 0                   # log2fc threshold to look for over/underexpressed genes along chr
                                      # it can be different from the threshold employed for short lists (less stringent!)
pval_thr_chrp   = pval_threshold

# loop comparisons
for (com in complis) {
  
  # load expression data for comparison
  tici = paste(com[5],".",com[2],"-",com[3],sep="")
  load(paste(outcode,"de_co",tici,"data.RData",sep="."))
  ex_res_i$gene = rownames(ex_res_i)
  
  # add over/underexpression info
  gf_isover_genes  = ex_res_i[ex_res_i$log2FC_shrink >  log2fc_thr_chrp & ex_res_i$pvalue<pval_thr_chrp,]$gene
  gf_isunder_genes = ex_res_i[ex_res_i$log2FC_shrink < -log2fc_thr_chrp & ex_res_i$pvalue<pval_thr_chrp,]$gene
  gf$isover  = gf$attributes %in% gf_isover_genes  * 1
  gf$isunder = gf$attributes %in% gf_isunder_genes * 1
  
  # count biased genes
  num_overex_genes = sum(gf$isover)
  num_underx_genes = sum(gf$isunder)
  num_biased_genes = sum(gf$isover)+sum(gf$isunder)
  num_total_genes  = nrow(gf)
  fra_overex_genes = sum(gf$isover)/nrow(gf)
  fra_underx_genes = sum(gf$isunder)/nrow(gf)
  fra_biased_genes = (sum(gf$isover)+sum(gf$isunder))/nrow(gf)
  
  # dataframe for output
  ex_gff_i = data.frame()
  
  pdf(file=paste(outcode,".de_co.",tici,"_chromosomalbias.pdf",sep=""),height=18,width=6)
  for (chrom in chromlist) {
    
    message("# expression along chr ",chrom,", comp ",tici)
    par(mfrow=c(6,1))
    
    # format 
    gf_i = gf[gf$seqid==chrom,]
    gf_i = merge(gf_i,ex_res_i,by.x="attributes","gene")
    gf_i = gf_i[!is.na(gf_i$log2FC_shrink),]
    gf_i = gf_i[order(gf_i[,"start"]),]
    
    # variables along windows (sliding functions)
    gf_i.roll              = data.frame(slide.min(gf_i$start,window=k_len,step=k_len*k_step))
    colnames(gf_i.roll)    = "start"
    gf_i.roll$end          = slide.max(gf_i$end,window=k_len,step=k_len*k_step)
    gf_i.roll$chrom        = chrom
    gf_i.roll$isover       = slide.sum(gf_i$isover,window=k_len,step=k_len*k_step)/k_len
    gf_i.roll$isunder      = slide.sum(gf_i$isunder,window=k_len,step=k_len*k_step)/k_len
    
    # fisher tests: overepresentation of overexpressed genes
    gf_i.roll$isover_fisher_p  = sapply(gf_i.roll$isover, function(x) fisher.test(matrix(data=c(
      x*k_len,                            # num genes overexpressed in win
      k_len-(x*k_len),                    # num genes not overexpressed in win
      num_overex_genes-x*k_len,           # num genes overexpresse outside 
      num_total_genes-(num_overex_genes-x*k_len)),   # num genes not overexpressed outside
      nrow = 2, ncol = 2))$p.value
    )
    gf_i.roll$isover_fisher_padj = p.adjust(gf_i.roll$isover_fisher_p,method = "BH")
    
    # fisher tests: overepresentation of underexpressed genes
    gf_i.roll$isunder_fisher_p = sapply(gf_i.roll$isunder, function(x) fisher.test(matrix(data=c(
      x*k_len,                            # num genes under in win
      k_len-(x*k_len),                    # num genes not under in win
      num_underx_genes-x*k_len,           # num genes under outside 
      num_total_genes-(num_underx_genes-x*k_len)),   # num genes not under outside
      nrow = 2, ncol = 2))$p.value)
    gf_i.roll$isunder_fisher_padj = p.adjust(gf_i.roll$isunder_fisher_p,method = "BH")
    
    # fisher tests: overepresentation of over+under genes
    gf_i.roll$isbias_fisher_p = sapply(gf_i.roll$isunder+gf_i.roll$isover, function(x) fisher.test(matrix(data=c(
      x*k_len,                            # num genes biased in win
      k_len-(x*k_len),                    # num genes not biased in win
      num_biased_genes-x*k_len,           # num genes biased outside 
      num_total_genes-(num_biased_genes-x*k_len)), # num genes not biased outside
      nrow = 2, ncol = 2))$p.value)
    gf_i.roll$isbias_fisher_padj = p.adjust(gf_i.roll$isbias_fisher_p,method = "BH")
    
    # genes in window
    gf_i.roll$genes = slide.returnstring(gf_i$attributes,sep = ",",window=k_len,step=k_len*k_step)

    
    # Plots!
    # plot fraction of over-expressed genes
    plot(gf_i.roll$start/1e6,  gf_i.roll$isover, type="l",col="green3",
         main=paste("Overexp genes chr",chrom,"| comp",tici),
         sub=paste("total genes=",nrow(gf_i),", win len=",k_len," genes, step=",k_step*100,"%",", log2FC>",log2fc_thr_chrp," pval<",pval_thr_chrp,sep=""),
         xlab="Mb",ylab="Fraction genes",
         ylim=c(0,1),las=1)
    abline(h=fra_overex_genes,lty=2,col="blue")
    legend("topright",legend = c("overexpressed",paste("genome mean =",signif(fra_overex_genes,3))),col=c("green3","blue"),
           lty = c(1,2),bty = "n")
    
    # plot pval Fisher over-expressed genes
    plot(gf_i.roll$start/1e6, -log(gf_i.roll$isover_fisher_padj,10), type="l",col="green3",
         main="Overexp Fisher's test",
         xlab="Mb",ylab="-log(p)",
         las=1,ylim=c(0,10))
    lines(gf_i.roll$start/1e6, -log(gf_i.roll$isover_fisher_p,10), col="slategray3")
    abline(h=-log(0.01,10),lty=2,col="red")
    legend("topright",legend = c("p = 0.01","padj BH","p"),col=c("red","green3","slategray3"),
           lty = c(2,1,1),bty = "n")
    
    # plot fraction of under-expressed genes
    plot(gf_i.roll$start/1e6,  gf_i.roll$isunder, type="l",col="magenta3",
         main=paste("Underexp genes chr",chrom,"| comp",tici),
         sub=paste("total genes=",nrow(gf_i),", win len=",k_len," genes, step=",k_step*100,"%",", log2FC<-",log2fc_thr_chrp," pval<",pval_thr_chrp,sep=""),
         xlab="Mb",ylab="Fraction genes",
         ylim=c(0,1),las=1)
    abline(h=fra_overex_genes,lty=2,col="blue")
    legend("topright",legend = c("underexpressed",paste("genome mean =",signif(fra_underx_genes,3))),col=c("magenta3","blue"),
           lty = c(1,2),bty = "n")
    
    # plot pval Fisher under-expressed genes
    plot(gf_i.roll$start/1e6, -log(gf_i.roll$isunder_fisher_padj,10), type="l",col="magenta3",
         main="Underexp Fisher's test",
         xlab="Mb",ylab="-log(p)",
         las=1,ylim=c(0,10))
    lines(gf_i.roll$start/1e6, -log(gf_i.roll$isunder_fisher_p,10), col="slategray3")
    abline(h=-log(0.01,10),lty=2,col="red")
    legend("topright",legend = c("p = 0.01","padj BH","p"),col=c("red","magenta3","slategray3"),
           lty = c(2,1,1),bty = "n")
    
    # plot fraction of over+under genes
    plot(gf_i.roll$start/1e6,  gf_i.roll$isover+gf_i.roll$isunder, type="l",col="slategray",
         main="DE genes biased (ov+un)",
         sub=paste("total genes=",nrow(gf_i),", win len=",k_len," genes, step=",k_step*100,"%",", |log2FC|>",log2fc_thr_chrp," pval<",pval_thr_chrp,sep=""),
         xlab="Mb",ylab="Fraction genes",
         ylim=c(0,1),las=1)
    abline(h=fra_biased_genes,lty=2,col="blue")
    legend("topright",legend = c("over + under",paste("genome mean =",signif(fra_biased_genes,3))),col=c("slategray","blue"),
           lty = c(1,2),bty = "n")
    
    # plot pval Fisher over+under genes
    plot(gf_i.roll$start/1e6, -log(gf_i.roll$isbias_fisher_padj,10), type="l",col="slategray",
         main="Biased Fisher's test",
         xlab="Mb",ylab="-log(p)",
         las=1,ylim=c(0,10))
    lines(gf_i.roll$start/1e6, -log(gf_i.roll$isbias_fisher_p,10), col="slategray3")
    abline(h=-log(0.01,10),lty=2,col="red")
    legend("topright",legend = c("p = 0.01","padj BH","p"),col=c("red","slategray","slategray3"),
           lty = c(2,1,1),bty = "n")
    
    ex_gff_i = rbind(ex_gff_i, gf_i.roll)
    
  }
  
  dev.off()
  
  # save table
  write.table(ex_gff_i,file=paste(outcode,".de_co.",tici,"_chromosomalbias.csv",sep=""),
              quote = F,row.names = F)
  
}



#### Save ####

message("# Save")

ex_res_tot = data.frame()
for (com in complis) {
  
  # define columns with groups
  if (com[5] == "ALL") {
    cli = as.character(sf[sf[,com[1]] %in% com[2],]$sample)
    clj = as.character(sf[sf[,com[1]] %in% com[3],]$sample)
  } else {
    cli = as.character(sf[sf[,com[1]] %in% com[2] & sf[,com[4]] %in% com[5],]$sample)
    clj = as.character(sf[sf[,com[1]] %in% com[3] & sf[,com[4]] %in% com[5],]$sample)
  }
  
  # load expression data for comparison
  tici = paste(com[5],".",com[2],"-",com[3],sep="")
  load(paste(outcode,"de_co",tici,"data.RData",sep="."))
  
  ex_dnc_f = na.omit(DESeq2::counts(ex_dat,normalized=T))
  
  # add columns
  ex_res_i$gene          = rownames(ex_res_i)
  ex_res_i$comparison    = tici
  ex_res_i$samples_i     = com[2]
  ex_res_i$normcounts_i  = rowMeans(ex_dnc_f[,cli],na.rm = T)
  ex_res_i$samples_j     = com[3]
  ex_res_i$normcounts_j  = rowMeans(ex_dnc_f[,clj],na.rm = T)
  ex_res_i               = merge(ex_res_i,dign,by.x="gene",by.y="gene_id")
  ex_res_i$is_signif_pos = ex_res_i$gene %in% ex_lis_i$genes_sp
  ex_res_i$is_signif_neg = ex_res_i$gene %in% ex_lis_i$genes_sn
  ex_res_i$is_signif_any = ex_res_i$is_signif_pos | ex_res_i$is_signif_neg
  
  write.table(ex_res_i,file=paste(outcode,"session.deseq_difexp",tici,"csv",sep="."),
              quote = T,row.names = F,sep = "\t")
  
  ex_res_tot             = rbind(ex_res_tot,ex_res_i)
  
}

write.table(ex_res_tot,file=paste(outcode,"session.deseq_difexp.csv",sep="."),
            quote = T,row.names = F,sep = "\t")

save.image(paste(outcode,"session.deseq_all.RData",sep="."))

message("\n\n### FI! ###\n\n")


stop("Ara!")



# 

load(paste(outcode,"de_co","ALL.RUN-SNO","data.RData",sep="."))
ex_res_i_SNO = ex_res_i
load(paste(outcode,"de_co","ALL.RUN-SRO","data.RData",sep="."))
ex_res_i_SRO = ex_res_i


interestlist_SNO = rownames(na.omit(ex_res_i_SNO[ex_res_i_SNO$padj<0.001,]))
interestlist_SRO = rownames(na.omit(ex_res_i_SRO[ex_res_i_SRO$padj<0.001,]))
interestlist_ALL = base::intersect(interestlist_SNO, interestlist_SRO)


interest_genes_heatmap(interestlist = interestlist_ALL, interestname = "ALLsignif",source_matrix = ex_dnc_st)


# function for plotting expression levels across samples with significance of DE
interest_genes_heatmap = function(interestlist, interestname, source_matrix) {
  
  original_len = length(interestlist)
  interestmatr = subset(source_matrix, rownames(source_matrix) %in% interestlist)
  interestlist = rownames(interestmatr)
  
  interestmatr.ann = data.frame(row.names = interestlist)
  interestmatr.exp = matrix(nrow = 0,ncol = 7)
  
  ann_colors = list()
  
  for (com in complis) {
    tici     = paste(com[5],".",com[2],"-",com[3],sep="")
    load(paste(outcode,"de_co",tici,"data.RData",sep="."))
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
  
  write.table(interestmatr.exp,file=paste("genesubset.de_adhoc",interestname,"csv",sep="."),
              quote = F,row.names = F, sep="\t")

  pdf(file=paste("genesubset.de_adhoc",interestname,"pdf",sep="."),height=10+length(interestlist)/4,width=10)
  pheatmap(interestmatr, color = col.fun(21),
           border_color = "white", annotation_row = interestmatr.ann,cellheight = 8,cellwidth = 8,cluster_cols=F,
           distance_rows = "correlation",scale="none",annotation_colors = ann_colors,
           gaps_col = seq(min(table(as.numeric(sf$group))),to=length(sf$group), by=min(table(as.numeric(sf$group)))),
           main=paste("DE: ",interestname," (expressed: ",length(interestlist)," out of ",original_len,")",sep=""))
  dev.off()
  
}



# detoxifying
interestlist = unique(pannu[ 
  pannu$domain == "GST_C" | pannu$domain == "GST_C_2" | pannu$domain == "GST_C_3" | pannu$domain == "GST_C_6" 
  | pannu$domain == "GST_N" | pannu$domain == "GST_N_3" | pannu$domain == "GST_N_4"  
  ,]$gene_id)
interest_genes_heatmap(interestlist = interestlist, interestname = "GST",source_matrix = ex_dnc_st)


interestlist = pannu[ pannu$domain == "p450" ,]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "CYP",source_matrix = ex_dnc_st)

interestlist = pannu[ pannu$domain == "COesterase" ,]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "COesterase",source_matrix = ex_dnc_st)

interestlist = pannu[ pannu$domain == "An_peroxidase" ,]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "An_peroxidase",source_matrix = ex_dnc_st)

interestlist = pannu[ pannu$domain == "ABC2_membrane" | pannu$domain == "ABC_tran" ,]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "ABC",source_matrix = ex_dnc_st)

interestlist = pannu[ pannu$domain == "FAD_binding_5" ,]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "dehydrogenase_FADbin5",source_matrix = ex_dnc_st)



# neur_chan_lbd
interestlist = pannu[ pannu$domain == "Neur_chan_LBD" ,]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "Neur_chan_LBD",source_matrix = ex_dnc_st)

# juvenile hormone binding
interestlist = pannu[ pannu$domain == "JHBP" ,]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "JHBP",source_matrix = ex_dnc_st)

# ion_trans
interestlist = pannu[ pannu$domain == "Ion_trans" ,]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "ion_trans",source_matrix = ex_dnc_st)

# chitin binding
interestlist = pannu[ pannu$domain == "Chitin_bind_4" ,]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "chitinbind4",source_matrix = ex_dnc_st)
interestlist = pannu[ pannu$domain == "CBM_14" ,]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "cbm14",source_matrix = ex_dnc_st)

# cuticle proteins
interestlist = pannu[ pannu$domain == "Cuticle_4" | pannu$domain == "Cuticle_3" | pannu$domain == "CPCFC",]$gene_id
interest_genes_heatmap(interestlist = interestlist, interestname = "Cuticle_3-4-CPCFC",source_matrix = ex_dnc_st)



# adhoc analyses: functions intersection
interestlist = ex_lis_joint[[1]]$overexpressed$list_intersect
interestlist = di[di$gene_id %in% interestlist, "long_transcript_id"]
hygeofun(list_interest=interestlist,annotation=pannu,gene_col="transcript",ano_col="domain",
         outputname=paste(outcode,"de_adhoc",sep="."),
         name_geneset=paste("RUN-SNOSRO_difpos",sep="_"),topnum = 50)
suppressMessages(topgofun(list_interest=interestlist,
                          gomap=gomap,
                          ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
                          outputname=paste(outcode,"de_adhoc",sep="."),
                          name_geneset=paste("RUN-SNOSRO_difpos",sep="_"),topnum=50))
write.table(interestlist,file=paste(outcode,"de_adhoc.RUN-SNOSRO_difpos.txt",sep="."),
            quote = F,row.names = F, col.names = F)

interestlist = ex_lis_joint[[1]]$underexpressed$list_intersect
interestlist = di[di$gene_id %in% interestlist, "long_transcript_id"]
hygeofun(list_interest=interestlist,annotation=pannu,gene_col="transcript",ano_col="domain",
         outputname=paste(outcode,"de_adhoc",sep="."),
         name_geneset=paste("RUN-SNOSRO_difneg",sep="_"),topnum = 50)
suppressMessages(topgofun(list_interest=interestlist,
                          gomap=gomap,
                          ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
                          outputname=paste(outcode,"de_adhoc",sep="."),
                          name_geneset=paste("RUN-SNOSRO_difneg",sep="_"),topnum=50))
write.table(interestlist,file=paste(outcode,"de_adhoc.RUN-SNOSRO_difneg.txt",sep="."),
            quote = F,row.names = F, col.names = F)



