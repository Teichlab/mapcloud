#preparing a worker cloud for picture day!
#(as in, the stuff that can be done pre-snapshot, so no disc mounting etc)

#log onto eta.internal.sanger.ac.uk, instances, launch instance, name it something informative
#source: basecloud (snapshot)
#flavor: m1.large (the smaller the better, the more sizes can use this later)
#networks: ensure cloudforms_network is dragged in
#security groups: ssh and icmp (on top of default)
#once spawned, press the little arrow on the far right of the instance's row and associate a floating IP

#cellranger prep! ho! version 2.0.2 with a GRCh38 1.2.0 reference, plus also mm10 and ERCC
#(and also set up 2.1.1 with the official VDJ download for VDJ purposes)
mkdir ~/cellranger && cd ~/cellranger
#originally this would be downloaded from 10X themselves, but they don't offer archival downloads
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/software/solexa/pkg/cellranger/2.0.2/dist/cellranger-2.0.2.tar.gz .
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/software/solexa/pkg/cellranger/2.1.1/dist/cellranger-2.1.1.tar.gz .
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/software/solexa/pkg/cellranger/3.0.2/dist/cellranger-3.0.2.tar.gz .
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0.tar.gz.download GRCh38.tar.gz
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-mm10-1.2.0.tar.gz.download mm10.tar.gz
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-ercc92-1.2.0.tar.gz.download ercc92.tar.gz
tar -xzvf cellranger-2.0.2.tar.gz && tar -xzvf cellranger-2.1.1.tar.gz && tar -xzvf cellranger-3.0.2.tar.gz
tar -xzvf GRCh38.tar.gz && tar -xzvf mm10.tar.gz && tar -xzvf ercc92.tar.gz && sudo rm *.tar.gz
curl -O http://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0.tar.gz
tar -xzvf refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0.tar.gz && rm refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0.tar.gz
#move the references to folders with shorter names for ease of script access later
mv refdata-cellranger-GRCh38-1.2.0/ GRCh38/
mv refdata-cellranger-mm10-1.2.0/ mm10/
mv refdata-cellranger-ercc92-1.2.0/ ERCC/
mv refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0/ GRCh38-VDJ/

#set stuff up for smartseq2 - STAR and HTSeq
#STAR lives in ~/cellranger/cellranger-2.0.2/STAR/2.5.1b/STAR
sudo pip3 install HTSeq

#create transcriptome files from 10X references for maximum consistency
#(along with a gene map and a translation of gene IDs to gene names)
sudo apt-get -y install cufflinks
cd ~ && mkdir transcriptomes+genemaps && cd transcriptomes+genemaps
for REF in GRCh38 mm10
do
	mkdir $REF && cd $REF
	gffread -w transcripts.fa -g ~/cellranger/$REF/fasta/genome.fa ~/cellranger/$REF/genes/genes.gtf
	grep -P '\ttranscript\t' ~/cellranger/$REF/genes/genes.gtf | sed -r 's/.*gene_id "([^"]+)".*transcript_id "([^"]+)".*/\2\t\1/' > genemap.tsv
	grep -P '\tgene\t' ~/cellranger/$REF/genes/genes.gtf | sed -r 's/.*gene_id "([^"]+)".*gene_name "([^"]+)".*/\1\t\2/' > genenames.tsv
	cd ..
done
#ERCCs don't get to have a complicated gtf, so need to build the genemap/genenames a little different
mkdir ERCC && cd ERCC
gffread -w transcripts.fa -g ~/cellranger/ERCC/fasta/genome.fa ~/cellranger/ERCC/genes/genes.gtf
grep -P '\texon\t' ~/cellranger/$REF/genes/genes.gtf | sed -r 's/.*gene_id "([^"]+)".*transcript_id "([^"]+)".*/\2\t\1/' > genemap.tsv
cp genemap.tsv genenames.tsv

#salmon!
mkdir ~/salmon && cd ~/salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.13.1/salmon-0.13.1_linux_x86_64.tar.gz
tar -xzvf salmon-0.13.1_linux_x86_64.tar.gz && rm salmon-0.13.1_linux_x86_64.tar.gz
mv salmon-latest_linux_x86_64/ salmon-0.13.1
salmon-0.13.1/bin/salmon index -t ~/transcriptomes+genemaps/GRCh38/transcripts.fa -i GRCh38
salmon-0.13.1/bin/salmon index -t ~/transcriptomes+genemaps/mm10/transcripts.fa -i mm10

#kallisto!
mkdir ~/kallisto && cd ~/kallisto
wget https://github.com/pachterlab/kallisto/releases/download/v0.45.0/kallisto_linux-v0.45.0.tar.gz
tar -xzvf kallisto_linux-v0.45.0.tar.gz && rm kallisto_linux-v0.45.0.tar.gz
mv kallisto_linux-v0.45.0 kallisto-0.45.0
kallisto-0.45.0/kallisto index -i GRCh38.idx ~/transcriptomes+genemaps/GRCh38/transcripts.fa
kallisto-0.45.0/kallisto index -i mm10.idx ~/transcriptomes+genemaps/mm10/transcripts.fa

#homer!
mkdir ~/homer && cd ~/homer && wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install mm10 hg38

#spatial transcriptomics! along with matplotlibby/tk'y dependencies it doesn't bother to tell you about!
sudo pip install stpipeline matplotlib
sudo apt-get -y install python-tk
#make st-qa.py use matplotlib in Agg mode
sudo sed '0,/^[^#]/{s/\(^[^#]\)/import matplotlib\nmatplotlib.use("Agg")\n\1/}' -i /usr/local/bin/st_qa.py
#this needs a new'ish STAR
cd ~ && wget https://github.com/alexdobin/STAR/archive/2.5.4b.tar.gz
tar -xzvf 2.5.4b.tar.gz && rm 2.5.4b.tar.gz
#for stpipeline purposes, STAR is in ~/STAR-2.5.4b/bin/Linux_x86_64/STAR

#grab gnomAD exome SNPs for genotyping pipeline
mkdir ~/gnomad-vcf && cd ~/gnomad-vcf
wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz
gunzip gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz
#this makes the non-major-chromosome names in GRCh38 compatible with the cellranger reference
sed -r 's/^[0-9]+\_(.*)v([0-9])\_[A-z]+/\1.\2/' -i gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf
#rename for ease of use, as always
mv gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf GRCh38.vcf

#increase file cap so demultiplexing 10X works in the genomics
echo '## Example hard limit for max opened files' | sudo tee -a /etc/security/limits.conf
echo 'ubuntu        hard nofile 131072' | sudo tee -a /etc/security/limits.conf
echo '## Example soft limit for max opened files' | sudo tee -a /etc/security/limits.conf
echo 'ubuntu        soft nofile 131072' | sudo tee -a /etc/security/limits.conf

#cellranger atac, not mirrored on farm so grab from source
#need to do this slow as this thing's filling up
mkdir ~/cellranger-atac && cd ~/cellranger-atac
curl -o cellranger-atac-1.1.0.tar.gz "http://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-1.1.0.tar.gz?Expires=1560210728&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWF0YWMvY2VsbHJhbmdlci1hdGFjLTEuMS4wLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU2MDIxMDcyOH19fV19&Signature=miadbrKyMfxpmoLisns4t~X0hzlmb5qLCT715pggdHeS3b6ablsBvpG7CFDg5WIzq6TzMA6-gVhkSZeaKeyF0XSr5wuMAY9T-0Ac5USTWb9Pi8OHfkXSpMky48Oh9niGAeDr4s72RejMmNMJS9zAp5zCsBtWZK0EA9ZtfkU9MeAX-SahQJP1I2NDil5lxImte0qtqZqZ9Q2Jc8JWZixp4vM1QOH5sgOSuUgso~NhtNw~LKErL0IpYGwYL7ZwPN8xC8RWwkoxtob0Yr-e~3PivRfzUqwONjC5L5~xicXhQ6jTgVzW-fehMZh8CsISRFjqc3S3d2Jqiaog53APoGQiSA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -xzvf cellranger-atac-1.1.0.tar.gz && rm cellranger-atac-1.1.0.tar.gz
curl -O http://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-1.1.0.tar.gz
tar -xzvf refdata-cellranger-atac-GRCh38-1.1.0.tar.gz && rm refdata-cellranger-atac-GRCh38-1.1.0.tar.gz
#I can't fit the mouse. sorry.
mv refdata-cellranger-atac-GRCh38-1.1.0 GRCh38

#the cloud is now ready for picture day!
#go back to eta, instances, create snapshot of the instance, name it something useful (mapcloud comes to mind)