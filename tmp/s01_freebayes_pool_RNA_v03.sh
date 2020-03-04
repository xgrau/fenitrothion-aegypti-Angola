#!/bin/bash

# define sample list
sample_list=$(cat $1)
echo $sample_list
if [ -z "$sample_list" ] ; then
	echo "NO SAMPLE LIST PROVIDED!"
	exit
fi

# num threads
nt=$2
if [ -z "$nt" ] ; then
	nt=8
fi

# output prefix
out_prefix=${3}
if [ -z "${out_prefix}" ] ; then
	out_prefix="out"
fi


# directories
dir="/home/xavi/dades/Transcriptomes/Aedaeg_Angola_19oct18/"
ref="${dir}/reference_transcripts/Aedaeg_long.cds.fasta"
generate_regions_py="scripts/fasta_generate_regions.py"
vcf_fix_header_py="scripts/vcffirstheader.py"

mkdir -p "${dir}/variants_gvcf_transcripts/"
mkdir -p "${dir}/variants_combined_transcripts/"

# list of input bam
bam_input_list=$( for s in $sample_list ; do 
	echo ${dir}/alignments_transcripts/${s}.bam
done )

# list of input raw vcf
vcf_input_list=$( for s in $sample_list ; do 
	echo ${dir}/variants_gvcf_transcripts/${s}.raw.vcf.gz
done )





## FUNCTIONS

# function: sample-spcific genotyping with freebayes
vcf_freebayes () {
	# generate regions to analyse in parallel, using GNU-parallel
	python3 $generate_regions_py ${ref}.fai 1000000 | grep -v "^UNKN" > ${dir}/variants_gvcf_transcripts/regions.list
   for s in $sample_list ; do
		if [ ! -f "${dir}/alignments_transcripts/${s}.bam.bai" ] ; then
			samtools index -@ ${nt} ${dir}/alignments_transcripts/${s}.bam
		fi
		if [ ! -f "${dir}/variants_gvcf_transcripts/${s}.raw.vcf.gz" ] ; then
			echo "$(date '+%Y-%m-%d %H:%M:%S') CREATE RAW VCF: ${dir}/variants_gvcf_transcripts/${s}.raw.vcf.gz"
			# create raw calls for all samples and all chroms
			parallel -a ${dir}/variants_gvcf_transcripts/regions.list -j ${nt} -k \
			freebayes --region {} -f ${ref} ${dir}/alignments_transcripts/${s}.bam --ploidy 10 --pooled-discrete --pooled-continuous \
			| python2 $vcf_fix_header_py  \
			| vcfstreamsort -w 1000 \
			| vcfsnps \
			| vcfuniq \
			> ${dir}/variants_gvcf_transcripts/${s}.raw.vcf.tmp1 \
			&& mv ${dir}/variants_gvcf_transcripts/${s}.raw.vcf.tmp1 ${dir}/variants_gvcf_transcripts/${s}.raw.vcf \
			&& bgzip ${dir}/variants_gvcf_transcripts/${s}.raw.vcf \
			&& tabix -p vcf ${dir}/variants_gvcf_transcripts/${s}.raw.vcf.gz
		else
			echo "$(date '+%Y-%m-%d %H:%M:%S') FOUND RAW VCF: ${dir}/variants_gvcf_transcripts/${s}.raw.vcf.gz"
		fi
	done
}


# combine all samples into a single VCF
vcf_combine () {
	split ${dir}/variants_gvcf_transcripts/regions.list ${dir}/variants_gvcf_transcripts/regions.list_split
	if [ ! -f "${dir}/variants_combined_transcripts/${out_prefix}_raw.vcf.gz" ] ; then
		echo "$(date '+%Y-%m-%d %H:%M:%S') COMBINE RAW VCFs: ${dir}/variants_combined_transcripts/${out_prefix}_raw.vcf.gz"
		for list in ${dir}/variants_gvcf_transcripts/regions.list_split* ; do
			sed "s/:\([0-9]*\)-\([0-9]*\)$/\t\1\t\2/" ${list} | awk '{ print $1"\t"$2+1"\t"$3+1 }' \
			> ${dir}/variants_gvcf_transcripts/regions.tmplist
			bcftools merge --threads ${nt} --regions-file ${dir}/variants_gvcf_transcripts/regions.tmplist \
			${vcf_input_list} --output-type z \
			--info-rules DP:sum,RO:sum,QR:sum,AO:sum,QA:sum
		done \
		> ${dir}/variants_combined_transcripts/${out_prefix}_raw.vcf.tmp.gz \
		&& mv ${dir}/variants_combined_transcripts/${out_prefix}_raw.vcf.tmp.gz ${dir}/variants_combined_transcripts/${out_prefix}_raw.vcf.gz \
		&& tabix -p vcf ${dir}/variants_combined_transcripts/${out_prefix}_raw.vcf.gz \
		&& rm ${dir}/variants_gvcf_transcripts/regions.list_split* ${dir}/variants_gvcf_transcripts/regions.tmplist
	else
		echo "$(date '+%Y-%m-%d %H:%M:%S') FOUND RAW VCFs: ${dir}/variants_combined_transcripts/${out_prefix}_raw.vcf.gz"
	fi
}

# apply quality filters to the joint VCF
vcf_filter () {
	if [ ! -f "${dir}/variants_combined_transcripts/${out_prefix}_recode.vcf.gz" ] ; then
		echo "$(date '+%Y-%m-%d %H:%M:%S') FILTER VCF: ${dir}/variants_combined_transcripts/${out_prefix}_recode.vcf.gz"
		bcftools filter \
		${dir}/variants_combined_transcripts/${out_prefix}_raw.vcf.gz \
		--include 'QUAL>30 & N_MISSING < 2 & TYPE="snp"' \
		--output-type z \
		--output ${dir}/variants_combined_transcripts/${out_prefix}_recode.vcf.tmp1.gz \
		&& mv ${dir}/variants_combined_transcripts/${out_prefix}_recode.vcf.tmp1.gz ${dir}/variants_combined_transcripts/${out_prefix}_recode.vcf.gz \
		&& tabix -p vcf ${dir}/variants_combined_transcripts/${out_prefix}_recode.vcf.gz

	else
		echo "$(date '+%Y-%m-%d %H:%M:%S') FOUND FILTERED VCF: ${dir}/variants_combined_transcripts/${out_prefix}_recode.vcf.gz"
	fi

}

## MAIN

vcf_freebayes \
&& vcf_combine \
&& vcf_filter


echo "Done"