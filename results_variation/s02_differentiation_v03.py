# define working dir
import os
os.chdir("/home/xavi/Documents/fenitrothion-aegypti-Angola/results_variation")

# import libraries
import allel
import h5py
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import itertools
import logging
import scipy.stats

# input files
species_id = "aeg"
listsam_fn = "../data_metadata/samples.list"
metasam_fn = "../data_metadata/samples_classification.csv"
callset_fn = "variants/out_recode.vcf.gz" # to retrieve genotypes with proper ploidy, we need VCF
callhd5_fn = "variants/out_recode.h5"     # to retrieve variant info, HDF5
anngene_fn = "../data_genome/Aedaeg_long.annot.gff"

# define population/group comparison
popl = ["RUN","RFE","SUS"]
popc = "resexp"

# general settings
sns.set(style="ticks",
        font_scale=1.3,
		  rc={"lines.linewidth": 1},
        font="Arial")

logging.basicConfig(
	level=logging.DEBUG, 
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s"
	# handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w"), logging.StreamHandler() ]
	)

## LOAD METADATA

# list of samples
listsam = pd.read_csv(listsam_fn, sep='\t', header=None)
listsam = listsam[0].values

# population metadata
samples = pd.read_csv(metasam_fn, sep='\t')
samples.columns = ["sample","group","resistance","exposure","colony","resexp"]

# load gene annotations: in pandas format, and in pyranges format
# pandas
anngene_pd = allel.gff3_to_dataframe(path=anngene_fn, attributes=["ID","Parent"])

## POPULATION DICTIONARIES

# population dictionary
popdict = dict()
popdict["all"] = []
for popi in popl: 
    popdict[popi]  = samples[samples[popc] == popi].index.tolist()
    popdict["all"] = popdict["all"] + popdict[popi]

## LOAD DATA

# dictionary of defaults to load genotypes and variants
# with scikit allel, with modified values to be able to 
# load polyploid samples
nondefault_numbers = {
    'variants/CHROM': 1, 'variants/POS': 1, 'variants/ID': 1, 'variants/REF': 1, 'variants/ALT': 'A',
    'variants/QUAL': 1,  'variants/DP': 1,  'variants/AN': 1, 'variants/AC': 'A', 'variants/AF': 'A',
    'variants/MQ': 1,    'variants/ANN': 1, 'calldata/DP': 1,
    'calldata/GT': 20, ## default is 2, but we need to use 10 beacuse we have ploidy = 10 (diploid * 5 samples)
    'calldata/GQ': 1, 'calldata/HQ': 2, 'calldata/AD': 'R', 'calldata/MQ0': 1, 'calldata/MQ': 1
}

def load_genotypes(callset_fn, callhd5_fn, popdict, genotype_tag="GT", segregating_pop="all"):
	# load callset
	callhd5 = h5py.File(callhd5_fn , mode="r")
	callset = allel.read_vcf(input=callset_fn, numbers=nondefault_numbers)
	genvars = allel.VariantChunkedTable(callhd5["variants"], names=["CHROM","POS","REF","ALT"])
	genotyp = allel.GenotypeChunkedArray(callset["calldata/%s" % genotype_tag])
	# calculate allele counts
	genalco = genotyp.count_alleles_subpops(subpops=popdict)
	# filters: segregating, no singletons
	is_seg    = genalco[segregating_pop].is_segregating()[:]   # segregating in population of interest? (default: any pop aka "all")
	is_nosing = genalco[segregating_pop][:,:2].min(axis=1) > 1 # no singletons
	is_retain = np.logical_and(is_seg, is_nosing)
	# apply filters
	genotyp_seg = genotyp.subset(sel0=is_retain)
	genalco_seg = genalco.compress(is_retain)
	genvars_seg = genvars[:].compress(is_retain)
	# logging
	logging.info("seg variants: %i / %i ( %.1f )" % ( genotyp_seg.shape[0], genotyp.shape[0], 100 * genotyp_seg.shape[0] / genotyp.shape[0]) )
	logging.info("samples: %i" % ( genotyp_seg.shape[1]) )
	# return
	return genotyp_seg, genalco_seg, genvars_seg

# load phased variants
genotyp_seg, genalco_seg, genvars_seg = load_genotypes(
	callset_fn=callset_fn, 
	callhd5_fn=callhd5_fn, 
	popdict=popdict, 
	genotype_tag="GT")


# function to calculate jack-knifed estiamte of PBS
def average_pbs(ac1, ac2, ac3, blen, normed):
	# calculate per-variant values
	vb = allel.pbs(ac1=ac1, ac2=ac2, ac3=ac3, normed=normed, window_size=blen)
	av = np.nanmean(vb)
	# estimate standard error
	_, se, vj = allel.stats.misc.jackknife( vb, statistic=lambda n: np.mean(n) )
	# output
	return av, se, vb, vj



# prepare transcript-wise output
gen = anngene_pd[anngene_pd["type"] == "mRNA"]
gen = gen.reset_index(drop=True)

# empty arrays for output
gen_num = np.full(len(gen), np.nan)
gen_pbs_av = np.full(len(gen), np.nan)
gen_pbs_se = np.full(len(gen), np.nan)
gen_fst_12_av = np.full(len(gen), np.nan)
gen_fst_13_av = np.full(len(gen), np.nan)
gen_fst_23_av = np.full(len(gen), np.nan)
gen_fst_12_se = np.full(len(gen), np.nan)
gen_fst_13_se = np.full(len(gen), np.nan)
gen_fst_23_se = np.full(len(gen), np.nan)
gen_tad_1 = np.full(len(gen), np.nan)
gen_tad_2 = np.full(len(gen), np.nan)
gen_tad_3 = np.full(len(gen), np.nan)
gen_pid_1 = np.full(len(gen), np.nan)
gen_pid_2 = np.full(len(gen), np.nan)
gen_pid_3 = np.full(len(gen), np.nan)

# loop through genes and populate output fields
for i,gei in gen.iterrows():

	# select variants in gene
	gei_bool = genvars_seg["CHROM"] == gei["ID"]
	
	# num variants in gene
	gen_num[i] = np.sum(gei_bool)

	if np.sum(gei_bool) > 0:

		# average PBS per gene + SE
		gen_pbs_av[i], gen_pbs_se[i], _ , _ = average_pbs(
			ac1=genalco_seg["RFE"].compress(gei_bool), 
			ac2=genalco_seg["RUN"].compress(gei_bool), 
			ac3=genalco_seg["SUS"].compress(gei_bool), 
			blen=1, normed=True)

		# average Fst per gene + SE
		gen_fst_12_av[i], gen_fst_12_se[i], _ , _ = allel.average_hudson_fst(
			ac1=genalco_seg["RFE"].compress(gei_bool), 
			ac2=genalco_seg["RUN"].compress(gei_bool), 
			blen=1)
		gen_fst_13_av[i], gen_fst_13_se[i], _ , _ = allel.average_hudson_fst(
			ac1=genalco_seg["RFE"].compress(gei_bool), 
			ac2=genalco_seg["SUS"].compress(gei_bool),
			blen=1)
		gen_fst_23_av[i], gen_fst_23_se[i], _ , _ = allel.average_hudson_fst(
			ac1=genalco_seg["RUN"].compress(gei_bool), 
			ac2=genalco_seg["SUS"].compress(gei_bool),
			blen=1)

		# average tajima D per gene
		gen_tad_1[i] = allel.tajima_d(
			ac=genalco_seg["RFE"].compress(gei_bool),
			pos=genvars_seg["POS"].compress(gei_bool))
		gen_tad_2[i] = allel.tajima_d(
			ac=genalco_seg["RUN"].compress(gei_bool),
			pos=genvars_seg["POS"].compress(gei_bool))
		gen_tad_3[i] = allel.tajima_d(
			ac=genalco_seg["SUS"].compress(gei_bool),
			pos=genvars_seg["POS"].compress(gei_bool))

		# average genetic diversity
		gen_pid_1[i] = allel.sequence_diversity(
			ac=genalco_seg["RFE"].compress(gei_bool),
			pos=genvars_seg["POS"].compress(gei_bool))
		gen_pid_2[i] = allel.tajima_d(
			ac=genalco_seg["RUN"].compress(gei_bool),
			pos=genvars_seg["POS"].compress(gei_bool))
		gen_pid_3[i] = allel.tajima_d(
			ac=genalco_seg["SUS"].compress(gei_bool),
			pos=genvars_seg["POS"].compress(gei_bool))

		

	# report progress
	if int(i+1) % int(len(gen)/100) == 0 or int(i+1) == len(gen): print(i+1,"/",len(gen))

# output
gen["nvar"] = gen_num
gen["pbs_av"] = gen_pbs_av
gen["pbs_zscore"] = allel.standardize(gen_pbs_av)
gen["pbs_se"] = gen_pbs_se
gen["fst12_av"] = gen_fst_12_av
gen["fst12_se"] = gen_fst_12_se
gen["fst13_av"] = gen_fst_13_av
gen["fst13_se"] = gen_fst_13_se
gen["fst23_av"] = gen_fst_23_av
gen["fst23_se"] = gen_fst_23_se
gen["tad1"] = gen_tad_1
gen["tad2"] = gen_tad_2
gen["tad3"] = gen_tad_3
gen["pid1"] = gen_pid_1
gen["pid2"] = gen_pid_2
gen["pid3"] = gen_pid_3

# write output table
gen.to_csv("differentiaton_%s.csv" % species_id, sep="\t", index=False)

