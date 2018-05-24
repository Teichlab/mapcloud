# Mapping/genotyping Teichlab computational cloud

This bit of text is going to detail the second available cloud snapshot, mapcloud, which builds upon [basecloud](https://github.com/Teichlab/basecloud) to include tools for analysing 10X/SS2/spatial transcriptomics data, both in terms of regular mapping/quantification and genotyping cells potentially from multiple donors. If you aren't familiar with basecloud, visit its repository first - it features a borderline excessive tutorial on using OpenStack, and everything in there carries over to here. The stuff exclusive to mapcloud is:

* **Cell Ranger** (v2.0.2 for internal consistency, v2.1.1 for VDJ) plus references: 10X GRCh38, mm10 and ERCC 1.2.0 releases, GRCh38 2.0.0 VDJ reference
	- the ERCC reference can be used to generate SS2-minded GRCh38/mm10+ERCC references:
	```
	cd ~/cellranger
	cat GRCh38/genes/genes.gtf ERCC/genes/genes.gtf > genes.gtf
	cat GRCh38/fasta/genome.fa ERCC/fasta/genome.fa > genome.fa
	cellranger-2.0.2/cellranger mkref --genome=GRCh38-ERCC --fasta=genome.fa --genes=genes.gtf --nthreads=`grep -c ^processor /proc/cpuinfo`
	rm genome.fa && rm genes.gtf
	
	cd ~/cellranger
	cat mm10/genes/genes.gtf ERCC/genes/genes.gtf > genes.gtf
	cat mm10/fasta/genome.fa ERCC/fasta/genome.fa > genome.fa
	cellranger-2.0.2/cellranger mkref --genome=mm10-ERCC --fasta=genome.fa --genes=genes.gtf --nthreads=`grep -c ^processor /proc/cpuinfo`
	rm genome.fa && rm genes.gtf
	```
* **salmon, kallisto** with GRCh38+ERCC and GRCm38+ERCC references
* **HTSeq** for STAR+HTSeq SS2 analysis; the pipeline uses the exact STAR version shipped with Cell Ranger and the references' SA indices for analysis consistency
* **st_pipeline** for STAR+HTSeq analysis of spatial transcriptomics data generated using SciLifeLab technology
* **VCF with SNPs for all protein coding genes** for use with 10X/SS2 GRCh38 mapping only
* **automated pipeline code base** for 10X/SS2/spatial transcriptomics analysis, plus per-cell donor genotype calling

All of this lives in the home directory (`~`), with self-explanatory folder names.

## Getting the data from iRODS

If you have a particular `$ID` of choice (a sample identifier for 10X or 10X-VDJ, or a RUN_LANE combination for SS2), the data lives in `/archive/HCA/$TECH/$ID`, where `$TECH` is 10X, 10X-VDJ or SS2 according to the technology that generated the data. In the case of 10X, the filtered count matrix can be downloaded as `iget -Kr /archive/HCA/10X/$ID/outs/filtered_gene_bc_matrices`. The main output for SS2 can be acquired via `iget -Kr /archive/HCA/SS2/$ID/outs`. The whole `$ID` folder can be downloaded via `iget -Kr` if the complete output is desired (akin to what was available in the LUSTRE folder before the iRODS migration), or looked at via `ils` to assess what exactly is available.

The main output for a space-delimited collection of 10X samples can be downloaded as follows:

	#!/bin/bash
	set -e

	#your space-delimited samples go here
	for SAMPLE in 
	do
		iget -Kr /archive/HCA/10X/$SAMPLE/outs/filtered_gene_bc_matrices $SAMPLE
		iget -K /archive/HCA/10X/$SAMPLE/outs/web_summary.html $SAMPLE
	done

### Using the pipelines

The 10X/SS2/spatial transcriptomics pipelines all exist as templates within the GitHub repository. All templates go from downloading data from iRODS to a complete mapping/quantification output. In the case of 10X/SS2, the complete output is automatically uploaded back to iRODS (`/archive/HCA` specifically) for storage. The spatial transcriptomics pipeline output is very minuscule and is not automatically uploaded to iRODS at this time. The genotyping pipeline can be used as an addition for 10X/SS2 runs, but only if the GRCh38 reference is used. It is a reimplementation of Davis/Raghd/Angela's picard/GATK approach with samtools/bcftools to improve run time.

To use the pipeline, copy over your desired template to a new directory on your mounted volume and edit the variables right at the start of the script. Select between a GRCh38 and mm10 reference in the `REFERENCE=` line. If running 10X data which has FASTQ files stored on iRODS, you can switch the file type away from CRAM if so desired in the `FILETYPE=` line. Paste your list of sample IDs (10X/spatial) or RUN_LANE combinations (SS2), separated by spaces, into the for loop and you're ready to go! Engaging the genotyping pipeline involves commenting out a single line from the middle of the templates, everything is clearly marked. Running EmptyDrops for 10X data and TPMs for SS2 data are enabled in a similar manner.

### Extra setup

You create a mapcloud just like you would a basecloud, but you select the mapcloud snapshot instead of basecloud in Zeta. While mapcloud was built on `m1.large`, the 10X/SS2 pipelines were always ran on larger flavours, with the Zeta equivalent being `m1.2xlarge`. The latter of these is recommended if you intend to make use of the pipelines, to the point where it was used to propose the increased allocation request as ten allotments of this size's resources.

You still have to go through all the same motions as with basecloud when setting one of these up, plus a quick `git clone` if you wish to use the pipelines contained within this repository.

	cd /mnt && git clone https://github.com/Teichlab/mapcloud

### 10X pipeline
* **Input:** Sanger study IDs, or any other uniquely identifying field in the iRODS metadata; GRCh38/mm10 choice of reference; cram/fastq file format to download; optional second pair of iRODS metadata fields (such as `run_id`) to help narrow down the download scope if needed; VCF file location/top expressed gene count to filter SNPs for genotyping if relevant
* Downloads CRAM/FASTQ from iRODS, converts CRAM to FASTQ if needed
* Runs Cell Ranger (version 2.0.2)
* Optionally runs EmptyDrops, creating a `final-count-matrix` folder with the Cell Ranger, EmptyDrops and cell set union count matrices in `.RDS` format, storing it in `outs`
* Optionally runs the genotyping pipeline, creating a `.vcf` file named after the sample in `outs`; only works with a GRCh38 reference
* Uploads complete output to `/archive/HCA/10X`, creating individual folders for each of your samples
* VDJ variant runs Cell Ranger 2.1.1 with a GRCh38 reference only and uploads the complete output to `/archive/HCA/10X-VDJ`, creating individual folders for each of your samples

### SS2 pipeline
* **Input:** RUN_LANE combinations identifying all your plates (e.g. 24013_1); GRCh38/mm10 choice of reference; VCF file location/top expressed gene count to filter SNPs for genotyping if relevant
* Downloads CRAM from iRODS, converts CRAM to FASTQ
* Downloads all of the iRODS metadata, extracts the `sample_supplier_name` field and creates a map of CRAM file to the resulting values in `outs/sampleInfo.txt`
* Aligns the data with STAR (version 2.5.1b, shipped with Cell Ranger 2.0.2)
* Quantifies the data with HTSeq
* Collects the counts and STAR mapping information for all the cells and creates:
	-	`outs/countMatrixNames.txt` tags the cells by the `sample_supplier_name` metadata
	-	`outs/countMatrix.txt` tags the cells by the CRAM file name
	-	`outs/uniqueMappedPercent.txt` captures the unique STAR mapping percentage
* Optionally computes TPMs, saving them as `tpm.RDS` and `tpmNames.RDS` in `outs`; only works with GRCh38 and mm10 references
* Optionally runs the genotyping pipeline, creating `snpCalls.vcf` and `snpCallsNames.vcf` file in `outs`; only works with a GRCh38 reference
* Uploads complete output to `/archive/HCA/SS2`, creating individual folders for each of your RUN_LANE combinations

### Spatial transcriptomics pipeline
* **Input:** Sanger study IDs, or any other uniquely identifying field in the iRODS metadata; GRCh38/mm10 choice of reference
* Downloads CRAM from iRODS, converts CRAM to FASTQ
* Runs st_pipeline to map the data with STAR (requiring a newer version than what's shipped with Cell Ranger 2.0.2) and quantify it with HTSeq
* Creates a second count matrix file, converting ENSEMBL IDs to gene names (tagged as `genenames`) and runs quality control on it
