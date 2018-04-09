# Mapping/genotyping Teichlab computational cloud

This bit of text is going to detail the second available cloud snapshot, mapcloud, which builds upon [basecloud](https://github.com/Teichlab/basecloud) to include tools for analysing 10X/SS2/spatial transcriptomics data, both in terms of regular mapping/quantification and genotyping cells potentially from multiple donors. If you aren't familiar with basecloud, visit its repository first - it features a borderline excessive tutorial on using OpenStack, and everything in there carries over to here. The stuff exclusive to mapcloud is:

* **Cell Ranger** (v2.0.2 for internal consistency) plus references: 10X GRCh38 1.2.0 release and a custom-built GRCh37.75 consistent with the above in terms of GTF `gene_biotype` filtering
* **salmon, kallisto** with GRCh38+ERCC and GRCm38+ERCC references
* **HTSeq** for STAR+HTSeq SS2 analysis; the pipeline uses the exact STAR version shipped with Cell Ranger and the references' SA indices for analysis consistency
* **st_pipeline** for STAR+HTSeq analysis of spatial transcriptomics data generated using SciLifeLab technology
* **VCF with SNPs for all protein coding genes** for use with 10X/SS2 GRCh37 mapping only
* **automated pipeline code base** for 10X/SS2/spatial transcriptomics analysis, plus per-cell donor genotype calling

### Using the pipelines

The 10X/SS2/spatial transcriptomics pipelines all exist as templates within the GitHub repository. All templates go from downloading data from iRODS to a complete mapping/quantification output. There are two modes for 10X and SS2: local and server. Local mode stores a selection of the output (the count matrix and some quality control report) for each sample in the directory you run it in, while server mode copies the complete output, including alignment BAMs and things to that nature, over to your specified remote location. There's only one spatial transcriptomics pipeline as it generates output matching what the local mode would want by default. The genotyping pipeline can be used as an addition for 10X/SS2 runs, but only if the GRCh37 reference is used. It is a reimplementation of Davis/Raghd/Angela's picard/GATK approach with samtools/bcftools to improve run time.

To use the pipeline, copy over your desired template to a new directory on your mounted volume and edit the variables right at the start of the script. Select between a GRCh37 and GRCh38 reference. If using the server mode, specify the remote path to copy to. If running 10X data which has FASTQ files stored on iRODS, you can switch the file type away from CRAM if so desired. Paste your list of sample IDs (10X/spatial) or RUN_LANE combinations (SS2), separated by spaces, into the for loop and you're ready to go! Engaging the genotyping pipeline involves commenting out a single line from the middle of the templates, everything is clearly marked. Running EmptyDrops for 10X data is enabled in a similar manner.

### Extra setup

You create a mapcloud just like you would a basecloud, but you select the mapcloud snapshot instead of basecloud in Delta. You can only use two flavours: `k1.2xlarge` (28 cores) or `s1.4xlarge` (54 cores). The former of these is recommended, to the point where it was used to propose the increased allocation request as ten allotments of this size's resources.

You still have to go through all the same motions as with basecloud when setting one of these up, plus a couple extra steps if you wish to use the pipelines contained within this repository (such as actually getting the code). The pipelines make use of the `$SSHNAME` system variable for minimum fuss automated result transfer to the farm, so if you foresee yourself using the mode wherein the complete output is automatically copied over to Lustre for safekeeping then call the code snippet below, subbing in your Sanger user ID, and store your farm login password in `~/.sshpass`. You're the only person on here so it's not going anywhere, don't worry.

	printf 'export SSHNAME=<user-id>\n' >> ~/.bashrc && exec bash
	cd /mnt && git clone https://github.com/Teichlab/mapcloud

### 10X pipeline
* **Input:** Sanger study IDs, or any other uniquely identifying field in the iRODS metadata; GRCh37/38 choice of reference; cram/fastq file format to download; optional second pair of iRODS metadata fields (such as `run_id`) to help narrow down the download scope if needed; VCF file location/top expressed gene count to filter SNPs for genotyping if relevant
* Downloads CRAM/FASTQ from iRODS, converts CRAM to FASTQ if needed
* Runs Cell Ranger (version 2.0.2)
* Optionally runs EmptyDrops, creating a `final-count-matrix` folder with the Cell Ranger, EmptyDrops and cell set union count matrices in `.RDS` format, storing it in `outs`
* Optionally runs the genotyping pipeline, creating a `.vcf` file named after the sample in `outs`; only works with a GRCh37 reference
* Local output stores the `filtered_gene_bc_matrices` folder as `cellranger` in a folder named after the sample, if EmptyDrops/genotyping were ran the output is copied into that folder too

### SS2 pipeline
* **Input:** RUN_LANE combinations identifying all your plates (e.g. 24013_1); GRCh37/38 choice of reference; VCF file location/top expressed gene count to filter SNPs for genotyping if relevant
* Downloads CRAM from iRODS, converts CRAM to FASTQ
* Downloads all of the iRODS metadata, extracts the `sample_supplier_name` field and creates a map of CRAM file to the resulting values in `outs/sampleInfo.txt`
* Aligns the data with STAR (version 2.5.1b, shipped with Cell Ranger 2.0.2)
* Quantifies the data with HTSeq
* Collects the counts and STAR mapping information for all the cells and creates:
	-	`outs/countMatrixNames.txt` tags the cells by the `sample_supplier_name` metadata
	-	`outs/countMatrix.txt` tags the cells by the CRAM file name
	-	`outs/uniqueMappedPercent.txt` captures the unique STAR mapping percentage
* Optionally computes TPMs, saving them as `tpm.RDS` and `tpmNames.RDS` in `outs`; only works with `GRCh38` and `mm10`
* Optionally runs the genotyping pipeline, creating `snpCalls.vcf` and `snpCallsNames.vcf` file in `outs`; only works with a GRCh37 reference
* Local output stores the `outs` folder and drops the STAR alignments, complete HTSeq output and iRODS metadata dump

### Spatial transcriptomics pipeline
* **Input:** Sanger study IDs, or any other uniquely identifying field in the iRODS metadata; GRCh37/38 choice of reference
* Downloads CRAM from iRODS, converts CRAM to FASTQ
* Runs st_pipeline to map the data with STAR (requiring a newer version than what's shipped with Cell Ranger 2.0.2) and quantify it with HTSeq
* Creates a second count matrix file, converting ENSEMBL IDs to gene names (tagged as `genenames`) and runs quality control on it
