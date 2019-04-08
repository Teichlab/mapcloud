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