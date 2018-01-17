#preparing a worker cloud for picture day!
#(as in, the stuff that can be done pre-snapshot, so no disc mounting etc)

#log onto delta.internal.sanger.ac.uk, instances, launch instance, name it something informative
#instance flavour: k1.2xlarge (28 cores) - best use of 116 available cores and 4 allocatable floating IPs
#boot from image, set image to trusty-isg-docker
#access & security tab, tick default + ssh + icmp
#networks tab, make sure cloudforms is dragged in
#once spawned, press the little arrow on the far right of the instance's row and associate a floating IP

#update to begin proceedings
sudo apt-get update

#set up irods
wget ftp://ftp.renci.org/pub/irods/releases/4.1.10/ubuntu14/irods-icommands-4.1.10-ubuntu14-x86_64.deb
wget ftp://ftp.renci.org/pub/irods/releases/4.1.10/ubuntu14/irods-runtime-4.1.10-ubuntu14-x86_64.deb
wget ftp://ftp.renci.org/pub/irods/releases/4.1.10/ubuntu14/irods-dev-4.1.10-ubuntu14-x86_64.deb
#this will cry about missing dependencies, just ignore it and call the apt-get -y -f install below
sudo dpkg -i irods-icommands-4.1.10-ubuntu14-x86_64.deb irods-runtime-4.1.10-ubuntu14-x86_64.deb irods-dev-4.1.10-ubuntu14-x86_64.deb
sudo apt-get -y -f install && rm *.deb
#nick configuration from the farm
rsync -Pr kp9@farm3-login.internal.sanger.ac.uk:/nfs/users/nfs_k/kp9/.irods ~
#delete the last line of the configuration file, which is some farm weirdness and doesn't apply here
sed ':a;N;$!ba;s/,\n    "irods_plugins_home" : "\/opt\/renci\/icommands\/plugins\/"//g' -i ~/.irods/irods_environment.json
#add internal.sanger.ac.uk to where this file looks for things
sudo sed 's/search openstacklocal/search openstacklocal internal.sanger.ac.uk/g' -i /etc/resolv.conf
#now just feed in your password. boom. done
iinit

#cellranger prep! ho! version 2.0.2 with a GRCh38 1.2.0 reference
mkdir ~/cellranger && cd ~/cellranger
#originally this would be downloaded from 10X themselves, but they don't offer archival downloads
rsync -P kp9@farm3-login.internal.sanger.ac.uk:/software/solexa/pkg/cellranger/2.0.2/dist/cellranger-2.0.2.tar.gz .
rsync -P kp9@farm3-login.internal.sanger.ac.uk:/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0.tar.gz.download GRCh38.tar.gz
tar -xzvf cellranger-2.0.2.tar.gz && tar -xzvf GRCh38.tar.gz && sudo rm *.tar.gz
#move the GRCh38 reference to a folder called GRCh38 for ease of script access later
mv refdata-cellranger-GRCh38-1.2.0/ GRCh38/

#custom-made GRCh37.75 reference. needs GTF and genome FASTA
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
rsync -P kp9@farm3-login.internal.sanger.ac.uk:/lustre/scratch117/cellgen/team205/FetalCellAtlas/remapV2/genome.fa .
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

#set up bamcollate2 and samtools for processing CRAM files
#bamcollate2 and libmaus2 are not super nice when compiled from source, but can be apt-getted apparently
sudo add-apt-repository -y ppa:gt1/staden-io-lib-trunk-tischler
sudo add-apt-repository -y ppa:gt1/libmaus2
sudo add-apt-repository -y ppa:gt1/biobambam2
sudo apt-get update && sudo apt-get -y install libmaus2-dev biobambam2

#bamcollate2 is slow and not particularly good at using resources
#so let's run a bunch of those in parallel to speed the proceedings using this thing up a notch. bam!
sudo apt-get -y install parallel

#we also need samtools, fresh ones with CRAM support
#and in turn, samtools wants htslib present
cd ~ && git clone https://github.com/samtools/htslib && cd htslib
sudo apt-get -y install libbz2-dev liblzma-dev libcurl4-openssl-dev
autoheader
autoconf
./configure
make
sudo make install
#DO NOT DELETE ~/htslib!! OR SAMTOOLS WILL CRY!!

#so now we can get samtools proper going
cd ~ && git clone https://github.com/samtools/samtools && cd samtools
sudo apt-get -y install libncurses5-dev
#ignore the warning autoheader spits out
autoheader
autoconf -Wno-syntax
./configure
make
sudo make install
#leave this just in case, you never know

#pre-download the CRAM cache as this is always going to be relevant
cd ~ && iget -K /seq/24013/24013_1#1.cram && samtools collate 24013_1#1.cram test && rm *am

#set stuff up for smartseq2 - STAR and HTSeq
#STAR lives in ~/cellranger/cellranger-2.0.2/STAR/2.5.1b/STAR
sudo easy_install pip && sudo pip install HTSeq

#toss in a version-controlled setup of R and Seurat just in case
#recent R setups require some non-standard apt-get repositories, so let's set them up
echo "deb https://cran.ma.imperial.ac.uk/bin/linux/ubuntu trusty/" | sudo tee -a /etc/apt/sources.list
#then add the corresponding apt key using this call
#sometimes fails because reasons, just call it again in that case
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
#good enough, let's install R but screw version controlling it as it refuses to install 3.4.2
sudo apt-get update && sudo apt-get -y install r-base

#deal with Seurat's mountain of dependencies, some of which are external and quiet about it
sudo apt-get -y install libxml2-dev openjdk-7-jdk
#and now setting up internal R dependencies!
sudo R -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("edgeR","BiocParallel"));  install.packages("devtools", repos="https://cran.ma.imperial.ac.uk/"); devtools::install_version("Seurat", version="2.1.0", repos="https://cran.ma.imperial.ac.uk/")'

#set up EmptyDrops
cd ~ && git clone https://github.com/TimothyTickle/hca-jamboree-cell-identification
cd hca-jamboree-cell-identification/src/poisson_model && sudo R CMD INSTALL package
cd ~ && sudo rm -r hca-jamboree-cell-identification

#configure git in case gitting needs to happen
git config --global user.name ktpolanski
git config --global user.email krzysztof_polanski@o2.pl

#park ssh password here
nano ~/.sshpass

#with that, the cloud is almost ready for picture day!
#just need to comment out a thing that happens when a new cloud gets spun up in /etc/fstab
#as if it isn't undone, making a cloud from the snapshot will try to find a mount, fail and die
sudo sed 's/\/dev\/vdb/#\/dev\/vdb/g' -i /etc/fstab

#go back to delta, instances, create snapshot of the instance, name it something useful (mapcloud here)
#once it saves, ssh back into the instance you snapshotted and undo the commenting out you just did
#if the instance becomes inaccessible, go back to delta and hard reboot it
sudo sed 's/#\/dev\/vdb/\/dev\/vdb/g' -i /etc/fstab
