#preparing a worker cloud for picture day!
#(as in, the stuff that can be done pre-snapshot, so no disc mounting etc)

#log onto zeta.internal.sanger.ac.uk, instances, launch instance, name it something informative
#source: basecloud (snapshot)
#flavor: m1.xlarge (the smaller the better, the more sizes can use this later)
#networks: ensure cloudforms_network is dragged in
#security groups: ssh and icmp (on top of default)
#once spawned, press the little arrow on the far right of the instance's row and associate a floating IP

#cellranger prep! ho! version 2.0.2 with a GRCh38 1.2.0 reference, plus also mm10 and ERCC
#(and also set up 2.1.1 with the official VDJ download for VDJ purposes)
mkdir ~/cellranger && cd ~/cellranger
#originally this would be downloaded from 10X themselves, but they don't offer archival downloads
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/software/solexa/pkg/cellranger/2.0.2/dist/cellranger-2.0.2.tar.gz .
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/software/solexa/pkg/cellranger/2.1.1/dist/cellranger-2.1.1.tar.gz .
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0.tar.gz.download GRCh38.tar.gz
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-mm10-1.2.0.tar.gz.download mm10.tar.gz
rsync -P <user-id>@farm3-login.internal.sanger.ac.uk:/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-ercc92-1.2.0.tar.gz.download ercc92.tar.gz
tar -xzvf cellranger-2.0.2.tar.gz && tar -xzvf cellranger-2.1.1.tar.gz
tar -xzvf GRCh38.tar.gz && tar -xzvf mm10.tar.gz && tar -xzvf ercc92.tar.gz && sudo rm *.tar.gz
curl -O http://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0.tar.gz
tar -xzvf refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0.tar.gz && rm refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0.tar.gz
#move the references to folders with shorter names for ease of script access later
mv refdata-cellranger-GRCh38-1.2.0/ GRCh38/
mv refdata-cellranger-mm10-1.2.0/ mm10/
mv refdata-cellranger-ercc92-1.2.0/ ERCC/
mv refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0/ GRCh38-VDJ/

#custom-made GRCh37.75 reference. needs GTF and genome FASTA (stuck with toplevel as PECAM1 lives on a patch)
#however, it is unuseable as STAR requires stupid amounts of RAM at toplevel
#as such, leaving the code chunk below as a historical curiosity more than anything. don't run it!
#engage mass comment mode just in case (the whole : <<'END' [...] END thing):
: <<'END'
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
wget http://ftp.ensemblorg.ebi.ac.uk/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz && gunzip Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
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
~/cellranger/cellranger-2.0.2/cellranger mkref --genome=GRCh37 --fasta=Homo_sapiens.GRCh37.75.dna.toplevel.fa \
                   --genes=Homo_sapiens.GRCh37.75.filtered.gtf --nthreads=`grep -c ^processor /proc/cpuinfo` --memgb=37
rm Homo_sapiens.GRCh37.75.dna.toplevel.fa && rm *.gtf
END

#set stuff up for smartseq2 - STAR and HTSeq
#STAR lives in ~/cellranger/cellranger-2.0.2/STAR/2.5.1b/STAR
sudo pip3 install HTSeq

#download GRCh38/GRCm38 transcriptomes from gencode, and tidy them up with a helper script
#(and also get ERCCs from github while we're doing that)
cd ~ && git clone https://github.com/Teichlab/mapcloud
mkdir transcriptomes+genemaps && cp -r mapcloud/setup/ERCC92 transcriptomes+genemaps
cd transcriptomes+genemaps && mkdir GRCh38 && mkdir GRCm38 && mkdir GRCh38-ERCC && mkdir GRCm38-ERCC
cd GRCh38 && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.transcripts.fa.gz
gunzip gencode.v28.transcripts.fa.gz && python3 ~/mapcloud/setup/gencode_parse.py gencode.v28.transcripts.fa
cd ../GRCm38 && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.transcripts.fa.gz
gunzip gencode.vM17.transcripts.fa.gz && python3 ~/mapcloud/setup/gencode_parse.py gencode.vM17.transcripts.fa
cd ..
cat GRCh38/transcripts.fasta ERCC92/transcripts.fasta > GRCh38-ERCC/transcripts.fasta
cat GRCm38/transcripts.fasta ERCC92/transcripts.fasta > GRCm38-ERCC/transcripts.fasta
cat GRCh38/genemap.tsv ERCC92/genemap.tsv > GRCh38-ERCC/genemap.tsv
cat GRCm38/genemap.tsv ERCC92/genemap.tsv > GRCm38-ERCC/genemap.tsv
#that's all we needed from the helper scripts
sudo rm -r ~/mapcloud

#salmon!
mkdir ~/salmon && cd ~/salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz
tar -xzvf Salmon-0.9.1_linux_x86_64.tar.gz && rm rm Salmon-0.9.1_linux_x86_64.tar.gz
mv Salmon-latest_linux_x86_64/ salmon-0.9.1
salmon-0.9.1/bin/salmon index -t ~/transcriptomes+genemaps/GRCh38-ERCC/transcripts.fasta -i GRCh38-ERCC
salmon-0.9.1/bin/salmon index -t ~/transcriptomes+genemaps/GRCm38-ERCC/transcripts.fasta -i GRCm38-ERCC

#kallisto!
mkdir ~/kallisto && cd ~/kallisto
wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
tar -xzvf kallisto_linux-v0.44.0.tar.gz && rm kallisto_linux-v0.44.0.tar.gz
mv kallisto_linux-v0.44.0 kallisto-0.44.0
kallisto-0.44.0/kallisto index -i GRCh38-ERCC.idx ~/transcriptomes+genemaps/GRCh38-ERCC/transcripts.fasta
kallisto-0.44.0/kallisto index -i GRCm38-ERCC.idx ~/transcriptomes+genemaps/GRCm38-ERCC/transcripts.fasta

#spatial transcriptomics! along with matplotlibby/tk'y dependencies it doesn't bother to tell you about!
sudo apt-get -y install python-dev
cd ~ && wget https://bootstrap.pypa.io/get-pip.py
sudo python get-pip.py && rm get-pip.py
#numpy still needs attention before everything else
sudo pip install numpy
sudo pip install stpipeline matplotlib
sudo apt-get -y install python-tk
#make st-qa.py use matplotlib in Agg mode
sudo sed '0,/^[^#]/{s/\(^[^#]\)/import matplotlib\nmatplotlib.use("Agg")\n\1/}' -i /usr/local/bin/st_qa.py
#this needs a new'ish STAR
wget https://github.com/alexdobin/STAR/archive/2.5.4b.tar.gz
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

#the cloud is now ready for picture day!
#go back to zeta, instances, create snapshot of the instance, name it something useful (mapcloud comes to mind)
