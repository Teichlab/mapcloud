#preparing a worker cloud for picture day!
#(as in, the stuff that can be done pre-snapshot, so no disc mounting etc)

#log onto theta.internal.sanger.ac.uk, instances, launch instance, name it something informative
#source: basecloud-theta-v1 (snapshot)
#flavor: m1.xlarge (the smaller the better, the more sizes can use this later)
#networks: ensure cloudforms_network is dragged in
#security groups: ssh and icmp (on top of default)
#once spawned, press the little arrow on the far right of the instance's row and associate a floating IP

#cellranger prep! ho! version 4.0.0, plus 3.0.0-based references and a 4.0.0 VDJ reference
mkdir ~/cellranger && cd ~/cellranger
#grab cellranger and spaceranger from source
curl -o cellranger-4.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-4.0.0.tar.gz?Expires=1603340087&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci00LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MDMzNDAwODd9fX1dfQ__&Signature=BuZqRlQ43lVaPSN-wBsOk5sRKEs9AyheT-bVLsyYZEXE2j-g5e~uLBu0laMFDi1VjZ1d~YkwQtzYGvphzx9vwQlHgZi4E~8HX-eyFQutaviTv5QuzeeB-60wgZ3BacJbwnK4co~J9mXX0S7zaPVlivpsskYyrvZUoXBnhvdrGhqX~JEPx8GPYOxR5Vps3J3otj0ctjezE-yuRz1pXR0Qyhy4zyZY318x42Fz3y8WiLG6jBTpoq9XKyeD7QIpxO4-EcaXwONBqP3yiG60b9U-tSllrHIxq7LC7OP7ChBHqEPuZQyeWekmqQdWbjUMbmOt9Vrau8akJL-Pa67TChXiWw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
curl -o spaceranger-1.1.0.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-1.1.0.tar.gz?Expires=1603344262&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvc3BhdGlhbC1leHAvc3BhY2VyYW5nZXItMS4xLjAudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjAzMzQ0MjYyfX19XX0_&Signature=Ma5Ab34hLcfVBWYyP2spPkleSc0jPVWdBaGQTxsg5rj~D0iKAuwcdiuCb~C8yQlsBdcjGbbp1Hl~jSxnM9CYYdjy6gK8yrDK4Ry~vptBsgpwDaJbxquwcnUeyirKJ8bGih~249Pra5U9IoZ42sDQxSAo1ahbhsUJrPH3psGQkDYgVL7e34X~lgg9A~agOTU8~d9KD98WMONPeKz9whSBL2-w2uvG2STcEiAtgB2Na84~jZw81fPQZODaZd6j812Z8CRvihy0FCIIjdgcowwIRdIGfESykd2jlMwXsA1N8oV3mwNStmtpFuoGmZyBYB22Te64Q6f9BKvpobKVPECY1A__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
#grab references, partly from official backup, partly from lustre folder to avoid duping up efforts
rsync -P <user-id>@farm5-login:/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-3.0.0.tar.gz.download GRCh38.tar.gz
rsync -Pr <user-id>@farm5-login:/lustre/scratch117/cellgen/team205/kp9/cellranger-refs/GRCh38-3.0.0-viral .
curl -O https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-4.0.0.tar.gz
tar -xzvf cellranger-4.0.0.tar.gz && rm cellranger-4.0.0.tar.gz
tar -xzvf spaceranger-1.1.0.tar.gz && rm spaceranger-1.1.0.tar.gz
tar -xzvf GRCh38.tar.gz && rm GRCh38.tar.gz
tar -xzvf refdata-cellranger-vdj-GRCh38-alts-ensembl-4.0.0.tar.gz && rm refdata-cellranger-vdj-GRCh38-alts-ensembl-4.0.0.tar.gz
mv refdata-cellranger-GRCh38-3.0.0 GRCh38-3.0.0
mv refdata-cellranger-vdj-GRCh38-alts-ensembl-4.0.0 GRCh38-VDJ-4.0.0

#STAR for starsoloing
cd ~ && wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
tar -xzvf 2.7.3a.tar.gz && rm 2.7.3a.tar.gz
mkdir ~/STAR-2.7.3a/whitelists && cd ~/STAR-2.7.3a/whitelists
cp ~/cellranger/cellranger-4.0.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz v3.txt.gz
gunzip v3.txt.gz
cp ~/cellranger/cellranger-4.0.0/lib/python/cellranger/barcodes/737K-august-2016.txt v2.txt
#references were created previously based on 10x's 3.0.0 with a stock star index call
cd ~/STAR-2.7.3a
rsync -Pr <user-id>@farm5-login:/lustre/scratch117/cellgen/team205/kp9/star-refs/GRCh38-3.0.0 .
rsync -Pr <user-id>@farm5-login:/lustre/scratch117/cellgen/team205/kp9/star-refs/GRCh38-3.0.0-viral .

#cellranger atac from source
mkdir ~/cellranger-atac && cd ~/cellranger-atac
curl -o cellranger-atac-1.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-1.2.0.tar.gz?Expires=1603344369&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hdGFjL2NlbGxyYW5nZXItYXRhYy0xLjIuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MDMzNDQzNjl9fX1dfQ__&Signature=iIqzedxT4T7eLqSKAyAJTA1c9hFCZimCWhaOBxRND1vtaNCQj4WT9pE-qhvNk7sfVSbCA9F43ASf0P8H2pplHGAnlzzOTEjIslIbCEdRVrXCQtpozWa05SMrMuVib4CxD6XaoChsz6b8SroO3MEmM2PRPQPdBin2g26CGWtrCQCUeKQV93nQlI5wJUINLLPDEqaQqWpkdnzdfVn~X3noNuuK~ohCh7CzacbcY3-1V1-jBAEWmvsbgxHThCOZhMm93FDn0JxG81cVv4uF6Rv11O7ch6O5GJie6p69zjgOj0z44YPXwCYd2KeyUOAf4xi-C4785cZ4CHZWmP8-tmWuyQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -xzvf cellranger-atac-1.2.0.tar.gz && rm cellranger-atac-1.2.0.tar.gz
curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-1.2.0.tar.gz
tar -xzvf refdata-cellranger-atac-GRCh38-1.2.0.tar.gz && rm refdata-cellranger-atac-GRCh38-1.2.0.tar.gz
mv refdata-cellranger-atac-GRCh38-1.2.0 GRCh38-1.2.0

#set up HTSeq for ss2
sudo pip3 install HTSeq

#create transcriptome files from 10X references for maximum consistency
#(along with a gene map and a translation of gene IDs to gene names)
sudo apt-get -y install cufflinks
cd ~ && mkdir transcriptomes+genemaps && cd transcriptomes+genemaps
for REF in GRCh38-3.0.0
do
	mkdir $REF && cd $REF
	gffread -w transcripts.fa -g ~/cellranger/$REF/fasta/genome.fa ~/cellranger/$REF/genes/genes.gtf
	grep -P '\ttranscript\t' ~/cellranger/$REF/genes/genes.gtf | sed -r 's/.*gene_id "([^"]+)".*transcript_id "([^"]+)".*/\2\t\1/' > genemap.tsv
	grep -P '\tgene\t' ~/cellranger/$REF/genes/genes.gtf | sed -r 's/.*gene_id "([^"]+)".*gene_name "([^"]+)".*/\1\t\2/' > genenames.tsv
	cd ..
done

#salmon!
mkdir ~/salmon && cd ~/salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.3.0/salmon-1.3.0_linux_x86_64.tar.gz
tar -xzvf salmon-1.3.0_linux_x86_64.tar.gz && rm salmon-1.3.0_linux_x86_64.tar.gz
mv salmon-latest_linux_x86_64/ salmon-1.3.0
salmon-1.3.0/bin/salmon index -t ~/transcriptomes+genemaps/GRCh38-3.0.0/transcripts.fa -i GRCh38-3.0.0

#kallisto!
mkdir ~/kallisto && cd ~/kallisto
wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz
tar -xzvf kallisto_linux-v0.46.1.tar.gz && rm kallisto_linux-v0.46.1.tar.gz
mv kallisto kallisto-0.46.1
kallisto-0.46.1/kallisto index -i GRCh38-3.0.0.idx ~/transcriptomes+genemaps/GRCh38-3.0.0/transcripts.fa

#homer!
mkdir ~/homer && cd ~/homer && wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install hg38

#souporcell!
mkdir ~/souporcell && cd ~/souporcell
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1tj2j8QZuGz8sylHgWbnejWyUn8n6m0Y8' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1tj2j8QZuGz8sylHgWbnejWyUn8n6m0Y8" -O souporcell.sif && rm -rf /tmp/cookies.txt
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=13aebUpEKrtjliyT9rYzRijtkNJVUk5F_' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=13aebUpEKrtjliyT9rYzRijtkNJVUk5F_" -O common_variants_grch38.vcf && rm -rf /tmp/cookies.txt

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
#go back to theta, instances, create snapshot of the instance, name it something useful (mapcloud comes to mind)