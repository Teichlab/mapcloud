# 10x/Smart-seq2 mapping cloud

This is a collection of scripts detailing the setup and use of a cloud for analysing both 10x and Smart-seq2 data in as compatible a manner as possible. This is accomplished by the Smart-seq2 pipeline making active use of the version of STAR shipped with Cell Ranger, along with the corresponding indices stored within the reference folders. Two references are currently available - Cell Ranger's stock GRCh38 (version 1.2.0) and a faithful recreation thereof, in terms of `gene_biotype` filtering of the GTF, from GRCh37.75. The cloud also comes with Seurat 2.1.0 for downstream analyses if desired.

Using the scripts is quite straightforward. Start off by following the setup instructions from the two scripts in the `setup` directory, which will create a functional snapshot of a worker cloud with all the software and references on it. You can then use this snapshot to spin up as many worker clouds as you have resources to support, with minimal setup required for each one. Once in possession of a functional worker cloud, pick an appropriate run script template (10x/ss2, select local output or complete farm output storage) and specify your samples/plates to run. In the case of server-copying templates, provide the correct path to the `rsync` that drops the files onto the farm. It's quite easy to spot.

**10X Pipeline:**
*	**Input:** Sanger study IDs, or any other uniquely identifying field in the iRODS metadata; GRCh37/38 choice of reference; cram/fastq iRODS file download; optional second pair of iRODS metadata fields (such as `run_id`) to help narrow down the download scope if needed
*	Downloads CRAM/FASTQ from iRODS, converts CRAM to FASTQ if needed
*	Runs Cell Ranger (version 2.0.2)
*	Optionally runs EmptyDrops, creating a `final-count-matrix` folder with the Cell Ranger, EmptyDrops and cell set union count matrices in `.RDS` forma, storing it in `outs`
*	Local output stores the `filtered_gene_bc_matrices` folder as `cellranger` in a folder named after the sample, if EmptyDrops was ran `final-count-matrix` is copied into that folder

**Smart-seq2 Pipeline:**
*	**Input:** RUN_LANE combinations identifying all your plates (e.g. 24013_1); GRCh37/38 choice of reference
*	Downloads CRAM from iRODS, converts CRAM to FASTQ
*	Downloads all of the iRODS metadata, extracts the `sample_supplier_name` field and creates a map of CRAM file to the resulting values in `outs/sampleInfo.txt`
*	Aligns the data with STAR (version 2.5.1b, shipped with Cell Ranger 2.0.2)
*	Quantifies the data with HTSeq
*	Collects the counts and STAR mapping information for all the cells and creates:
	-	`outs/countMatrixNames.txt` tags the cells by the `sample_supplier_name` metadata
	-	`outs/countMatrix.txt` tags the cells by the CRAM file name
	-	`outs/uniqueMappedPercent.txt` captures the unique STAR mapping percentage
*	Local output stores the `outs` folder and drops the STAR alignments, complete HTSeq output and iRODS metadata dump
