#preparing a worker cloud for picture day!
#(as in, the stuff that can be done pre-snapshot, so no disc mounting etc)

#log onto delta.internal.sanger.ac.uk, instances, launch instance, name it something informative
#instance flavour: k1.2xlarge (28 cores, loads of drive space for references etc)
#boot from snapshot, set snapshot to basecloud
#access & security tab, tick default + ssh + icmp
#networks tab, make sure cloudforms is dragged in
#once spawned, press the little arrow on the far right of the instance's row and associate a floating IP

#cellranger prep! ho! version 2.0.2 with a GRCh38 1.2.0 reference
mkdir ~/cellranger && cd ~/cellranger
#originally this would be downloaded from 10X themselves, but they don't offer archival downloads
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/software/solexa/pkg/cellranger/2.0.2/dist/cellranger-2.0.2.tar.gz .
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0.tar.gz.download GRCh38.tar.gz
tar -xzvf cellranger-2.0.2.tar.gz && tar -xzvf GRCh38.tar.gz && sudo rm *.tar.gz
#move the GRCh38 reference to a folder called GRCh38 for ease of script access later
mv refdata-cellranger-GRCh38-1.2.0/ GRCh38/

#custom-made GRCh37.75 reference. needs GTF and genome FASTA
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/lustre/scratch117/cellgen/team205/FetalCellAtlas/remapV2/genome.fa .
#filter the GTF to begin with, then construct the reference
~/cellranger/cellranger-2.0.2/cellranger mkgtf Homo_sapiens.GRCh37.75.gtf Homo_sapiens.GRCh37.75.filtered.gtf \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lincRNA \
                   --attribute=gene_biotype:antisense \
                   --attribute=gene_biotype:IG_C_gene \
                   --attribute=gene_biotype:IG_C_pseudogene \
                   --attribute=gene_biotype:IG_D_gene \
                   --attribute=gene_biotype:IG_J_gene \
                   --attribute=gene_biotype:IG_J_pseudogene \
                   --attribute=gene_biotype:IG_V_gene \
                   --attribute=gene_biotype:IG_V_pseudogene \
                   --attribute=gene_biotype:TR_V_gene \
                   --attribute=gene_biotype:TR_V_pseudogene \
                   --attribute=gene_biotype:TR_J_gene \
                   --attribute=gene_biotype:TR_J_pseudogene \
                   --attribute=gene_biotype:TR_D_gene \
                   --attribute=gene_biotype:TR_C_gene
~/cellranger/cellranger-2.0.2/cellranger mkref --genome=GRCh37 --fasta=genome.fa --genes=Homo_sapiens.GRCh37.75.filtered.gtf --nthreads=28
rm genome.fa && rm *.gtf

#bamcollate2 is slow and not particularly good at using resources
#so let's run a bunch of those in parallel to speed the proceedings using this thing up a notch. bam!
sudo apt-get -y install parallel

#set stuff up for smartseq2 - STAR and HTSeq
#STAR lives in ~/cellranger/cellranger-2.0.2/STAR/2.5.1b/STAR
sudo pip3 install HTSeq

#download GRCh38/GRCm38 transcriptomes
cd ~ && rsync -Pr <user-id>@farm3-login.internal.sanger.ac.uk:/lustre/scratch117/cellgen/team269/kp9/cdna+genemaps .

#salmon!
mkdir ~/salmon && cd ~/salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz
tar -xzvf Salmon-0.9.1_linux_x86_64.tar.gz && rm rm Salmon-0.9.1_linux_x86_64.tar.gz
mv Salmon-latest_linux_x86_64/ salmon-0.9.1
salmon-0.9.1/bin/salmon index -t ~/cdna+genemaps/GRCh38-ERCC/GRCh38.p10_E89_cDNA.fasta -i GRCh38-ERCC
salmon-0.9.1/bin/salmon index -t ~/cdna+genemaps/GRCm38-ERCC/GRCm38.p5_E89_cDNA.fasta -i GRCm38-ERCC

#kallisto!
mkdir ~/kallisto && cd ~/kallisto
wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
tar -xzvf kallisto_linux-v0.44.0.tar.gz && rm kallisto_linux-v0.44.0.tar.gz
mv kallisto_linux-v0.44.0 kallisto-0.44.0
kallisto-0.44.0/kallisto index -i GRCh38-ERCC.idx ~/cdna+genemaps/GRCh38-ERCC/GRCh38.p10_E89_cDNA.fasta
kallisto-0.44.0/kallisto index -i GRCm38-ERCC.idx ~/cdna+genemaps/GRCm38-ERCC/GRCm38.p5_E89_cDNA.fasta

#spatial transcriptomics! along with matplotlibby/tk'y dependencies it doesn't bother to tell you about!
sudo apt-get -y install python-dev
cd ~ && wget https://bootstrap.pypa.io/get-pip.py
sudo python get-pip.py && rm get-pip.py
#numpy still needs attention before everything else
sudo pip install numpy
sudo pip install stpipeline matplotlib
sudo apt-get -y install python-tk
echo 'backend : Agg' | sudo tee -a ~/.config/matplotlib/matplotlibrc
wget https://github.com/alexdobin/STAR/archive/2.5.4b.tar.gz
tar -xzvf 2.5.4b.tar.gz && rm 2.5.4b.tar.gz
#for stpipeline purposes, STAR is in ~/STAR-2.5.4b/bin/Linux_x86_64/STAR
#grab the repository for design barcode files
cd ~ && git clone https://github.com/SpatialTranscriptomicsResearch/st_pipeline

#grab protein coding SNPs off farm for genomic pipeline
mkdir ~/vcf-protein-coding && cd ~/vcf-protein-coding
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/lustre/scratch117/cellgen/team205/singlecell_fetal/data/references/dbsnp_138.hg19.protein_coding.recode.dedupped.vcf.gz .
gunzip dbsnp_138.hg19.protein_coding.recode.dedupped.vcf.gz
#kick out all the chr to make it be consistent with the cellranger reference
sed 's/chr//g' -i dbsnp_138.hg19.protein_coding.recode.dedupped.vcf

#with that, the cloud is almost ready for picture day!
#just need to comment out a thing that happens when a new cloud gets spun up in /etc/fstab
#as if it isn't undone, making a cloud from the snapshot will try to find a mount, fail and die
sudo sed 's/\/dev\/vdb/#\/dev\/vdb/g' -i /etc/fstab

#go back to delta, instances, create snapshot of the instance, name it something useful (mapcloud comes to mind)
#once it saves, ssh back into the instance you snapshotted and undo the commenting out you just did
#if the instance becomes inaccessible, go back to delta and hard reboot it
sudo sed 's/#\/dev\/vdb/\/dev\/vdb/g' -i /etc/fstab
