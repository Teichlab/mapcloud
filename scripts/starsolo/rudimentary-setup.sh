#set up as follows
#(needed to regenerate star index as new star versions that support solo
#have some sort of slightly different index thing going on)
~/STAR-2.7.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --runThreadN 26 --genomeDir ~/STAR-2.7.3a/GRCh38 --genomeFastaFiles ~/cellranger/GRCh38/fasta/genome.fa --sjdbGTFfile ~/cellranger/GRCh38/genes/genes.gtf
mkdir ~/STAR-2.7.3a/whitelists && cd ~/STAR-2.7.3a/whitelists
cp ~/cellranger/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt.gz v3.txt.gz
gunzip v3.txt.gz
cp ~/cellranger/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/737K-august-2016.txt v2.txt