# Genetic variants from pooled RNA-seq samples

## Mapping & calling

### Requirements

Tools:

* Mapping: `bwa` 0.7.17-r1188 (`bwa mem` aligner)
* Read processing: `biobambam2` (in `conda`), for marking duplicate reads and sorting.
* VCF calling: [`freebayes`](https://github.com/ekg/freebayes) v1.3.1-dirty, combined with GNU `parallel`.
* VCF processing: `vcflib` and (in `conda`).
* Phasing: `phaseit` (including `extractPIRs`) // NOT USEFUL FOR RNA-seq.
* Python 3.

### Steps

1. Create folder structure, use longest isoform per transcript from [Vectorbase](https://www.vectorbase.org/downloads?field_organism_taxonomy_tid%5B%5D=367&field_download_file_format_tid=All&field_status_value=Current) and store it in `${dir}/reference` folder:

   ```bash
   # home directory where you saved your raw data. Should contain an `input` folder with raw reads (fastq)
   dir="/home/xavi/dades/Transcriptomes/Aedaeg_Angola_19oct18/"
   # other folders
   mkdir -p ${dir}/reference_transcripts # save reference fasta here
   ref="${dir}/reference_transcripts/Aedaeg_long.cds.fasta" # this is the reference
   mkdir -p ${dir}/alignments_transcripts # where BAM alignments will go
   nt=10 # number of threads (bwa and freebayes)
   ```

2. Index reference genome:

   ```bash
   bwa index ${ref}
   samtools faidx ${ref}
   /home/xavi/Programes/gatk-4.1.2.0/gatk CreateSequenceDictionary -R ${ref} -O ${ref%%.fa}.dict
   ```

3. Map each sample with `bwa` (run separately for each sample in `../data_metadata/samples.list`), mark duplicates and merge runs with `biobambam`. For `freebayes` genotyping, no base recalibration or indel realignment steps are needed ([here why](https://github.com/ekg/freebayes#calling-variants-from-fastq-to-vcf)).

   ```bash
   s={SAMPLE}
   echo "### $s ###"
   bwa mem -M -t ${nt} ${ref} ${dir}/input/Aedaeg_${s}_1.fastq.gz ${dir}/input/Aedaeg_${s}_2.fastq.gz -R "@RG\tID:${s}\tSM:${s}\tPU:nye\tPL:nye\tLB:${s}" | bamsormadup inputformat=sam threads=${nt} tmpfile=${dir}/alignments_transcripts/tmp_$(date +%s) SO=coordinate > ${dir}/alignments_transcripts/${s}.tmp1.bam && mv ${dir}/alignments_transcripts/${s}.tmp1.bam ${dir}/alignments_transcripts/${s}.bam && rm ${dir}/alignments_transcripts/tmp_*

   ```

4. VCF calling with `freebayes`, filters with `bcftools`.

   ```bash
   bash s01_freebayes_pool_RNA_v02.sh ../data_metadata/samples.list ${nt}
   ```

5. Convert VCF to ZARR and/or HDF5:

   ```python
   import allel
   allel.vcf_to_hdf5('out_recode.vcf.gz', 'out_recode.h5', exclude_fields=['variants/numalt'], overwrite=True)
   ```

6. Calculate gene-wise genetic differentiation statistics ($F_{ST}$ and $PBS$) for each species using `scikit-allel`, which will be stored in the `differentiation_sps.csv` files:

   ```bash
   python s02_differentiation_v03.py
   ```

7. Find top PBS genes in the resistant population, perform functional enrichment analyses:

   ```bash
   Rscript s03_genewise_analysis_v02.R
   ```