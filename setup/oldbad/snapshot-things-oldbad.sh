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

#old gencode-based transcriptomes
#figured out a way to make transcriptomes out of 10x references, so this got obsoleted.

#download GRCh38/GRCm38 transcriptomes from gencode, and tidy them up with a helper script
#(and also get ERCCs from github while we're doing that)
cd ~ && git clone https://github.com/Teichlab/mapcloud
mkdir transcriptomes+genemaps && cp -r mapcloud/setup/oldbad/ERCC92 transcriptomes+genemaps
cd transcriptomes+genemaps && mkdir GRCh38 && mkdir GRCm38 && mkdir GRCh38-ERCC && mkdir GRCm38-ERCC
cd GRCh38 && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.transcripts.fa.gz
gunzip gencode.v28.transcripts.fa.gz && python3 ~/mapcloud/setup/oldbad/gencode_parse.py gencode.v28.transcripts.fa
cd ../GRCm38 && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.transcripts.fa.gz
gunzip gencode.vM17.transcripts.fa.gz && python3 ~/mapcloud/setup/oldbad/gencode_parse.py gencode.vM17.transcripts.fa
cd ..
cat GRCh38/transcripts.fasta ERCC92/transcripts.fasta > GRCh38-ERCC/transcripts.fasta
cat GRCm38/transcripts.fasta ERCC92/transcripts.fasta > GRCm38-ERCC/transcripts.fasta
cat GRCh38/genemap.tsv ERCC92/genemap.tsv > GRCh38-ERCC/genemap.tsv
cat GRCm38/genemap.tsv ERCC92/genemap.tsv > GRCm38-ERCC/genemap.tsv
#that's all we needed from the helper scripts
sudo rm -r ~/mapcloud

#OLD CELLRANGER. now just roll with 4.0.0
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

#ERCC GENEMAP STUFF.
#ERCCs don't get to have a complicated gtf, so need to build the genemap/genenames a little different
mkdir ERCC && cd ERCC
gffread -w transcripts.fa -g ~/cellranger/ERCC/fasta/genome.fa ~/cellranger/ERCC/genes/genes.gtf
grep -P '\texon\t' ~/cellranger/$REF/genes/genes.gtf | sed -r 's/.*gene_id "([^"]+)".*transcript_id "([^"]+)".*/\2\t\1/' > genemap.tsv
cp genemap.tsv genenames.tsv

#OLD SPATIAL. now just roll with visium as this got antiquated.
#spatial transcriptomics! along with matplotlibby/tk'y dependencies it doesn't bother to tell you about!
sudo pip install stpipeline matplotlib
sudo apt-get -y install python-tk
#make st-qa.py use matplotlib in Agg mode
sudo sed '0,/^[^#]/{s/\(^[^#]\)/import matplotlib\nmatplotlib.use("Agg")\n\1/}' -i /usr/local/bin/st_qa.py
#this needs a new'ish STAR
cd ~ && wget https://github.com/alexdobin/STAR/archive/2.5.4b.tar.gz
tar -xzvf 2.5.4b.tar.gz && rm 2.5.4b.tar.gz
#for stpipeline purposes, STAR is in ~/STAR-2.5.4b/bin/Linux_x86_64/STAR