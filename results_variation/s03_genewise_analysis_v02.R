#### Input ####

# genewise differentiation
gen_fn = "differentiaton_aeg.csv"

# genome annotations
gomapfi = "../data_genome/Aedaeg_eggnog_dipNOG.emapper.annotations.GO.csv"
panfile = "../data_genome/Aedaeg_long.pep_Pfamscan.seqs"

# where to store output?
outcode = "differentiaton_gam_out" # (folder + prefix)

# load libraries
library(topGO)
library(fdrtool)
source("../helper_scripts/geneSetAnalysis.R")

# chromosomes to include
chromlist = c("2","3","1")

# genewise differentiation statistics
gen = read.table(gen_fn, sep="\t", header = T)
gen = gen[gen$seqid %in% chromlist,]
# gen = gen[gen$nvar > 10,]

# functional mappings: GO, pfam
gomap = readMappings(gomapfi)
panno = read.table(file = panfile)
colnames(panno) = c("transcript","pstart","pend","pfamid","domain","domseq")
pannu = subset(panno, select=c("transcript","pfamid","domain"))
pannu = pannu[!duplicated(pannu), ]


#### top PBS ####

# remove genes that are invariant in at least one of the
# populations comparisons (1 to 2, 1 to 3, 2 to 3); i.e., PBS = NA
gen_variant = gen[!is.na(gen$pbs_av),]

# adjusted pvalue from standardised PBS
gen_variant$PBS_p_adj = fdrtool(gen_variant$pbs_zscore, statistic = 'normal', plot = F)$qval


pdf(sprintf("%s.dist.pdf", outcode), width = 5, height = 4)
# plot standardised PBS distribution
hist(gen_variant$pbs_zscore, 
     col="blue", border="white",
     breaks = 40, xlab="standardised PBS",
     main=sprintf("Standardised PBS\nn=%i/%i genes",nrow(gen_variant), nrow(gen)))

# plot standardised PBS along chromosome (per gene)
for (chrom in chromlist) {
  plot(x=gen_variant[gen_variant$seqid==chrom,"start"]/1e6, y=gen_variant[gen_variant$seqid==chrom,"pbs_zscore"], col="blue",
       xlim=c(0,500), ylim=c(-10,10), cex=0.5,
       xlab="Mb",ylab="standardised PBS",
       main=sprintf("%s",chrom),las=1)
  abline(h=0,lty=2)
}
dev.off()


#### Functional enrichment ####

# genes with the most extreme highest pbs values
list_genes = as.vector(gen_variant[gen_variant$PBS_p_adj < 0.001, "ID"])
list_genes = list_genes[!is.na(list_genes)]
# functional enrichments
hygeofun(list_interest=list_genes, 
         annotation=pannu, gene_col="transcript", ano_col="domain",
         outputname=outcode,
         name_geneset="PBStop",topnum = 30, padjbool = F)
suppressMessages(topgofun(list_interest=list_genes, 
                          gomap=gomap,
                          ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
                          outputname=outcode,
                          name_geneset="PBStop",topnum=30))

message("Fi!")
