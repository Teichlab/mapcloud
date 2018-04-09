#!/usr/bin/Rscript

library(methods)

#function for getting annotations off biomaRt
getAnnotations = function(genes,gene_ids='ensembl_gene_id',attributes=c("ensembl_gene_id","chromosome_name","gene_biotype", "start_position", "end_position","description"),dataset='hsapiens_gene_ensembl',biomart='ENSEMBL_MART_ENSEMBL',host='www.ensembl.org'){
  if(is.null(host)){
    bmart = biomaRt::useMart(biomart=biomart,dataset=dataset)
  }else{
    bmart = biomaRt::useMart(biomart=biomart,dataset=dataset,host=host)
  }
  dat = biomaRt::getBM(attributes=attributes,filters=gene_ids,values=genes,mart=bmart)
  return(dat)
}

getTPM = function(runlane, dataset)
{
	#read count matrix and soak up the matching annotations
	tpm = read.csv(paste(runlane,'/outs/countMatrix.txt',sep=''),row.names=1,header=TRUE,sep='\t',check.names=FALSE)
	geneDat = getAnnotations(rownames(tpm),dataset=dataset)
	#sometimes stuff gets duplicated in the geneDat, so just keep the first occurrence of each gene ID
	geneDat = geneDat[!duplicated(geneDat[,1]),]
	#ensure that we've got the same gene order in both the count matrix and the annotation
	#need this as a two-step as order(match()) doesn't like unmapped values
	tpm = tpm[rownames(tpm) %in% geneDat[,1],]
	tpm = tpm[order(match(rownames(tpm),geneDat[,1])),]
	#easier TPM computation, as our indices always match
	#start with RPK and then convert to TPM
	#http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
	len = geneDat[,5] - geneDat[,4] + 1
	tpm = tpm/len
	tpm = t(t(tpm)/(colSums(tpm)/1000000))
	#save the TPM, and then convert the column names to the ones in countMatrixNames.txt
	saveRDS(tpm,paste(runlane,'/outs/tpm.RDS',sep=''))
	names = read.csv(paste(runlane,'/outs/countMatrixNames.txt',sep=''),row.names=1,header=TRUE,sep='\t',check.names=FALSE)
	colnames(tpm) = colnames(names)
	saveRDS(tpm,paste(runlane,'/outs/tpmNames.RDS',sep=''))
}

args <- commandArgs(TRUE)
runlane <- args[1]
reference <- args[2]

#set up dataset of choice
#remove -ERCC from reference if it shows up
reference = gsub('-ERCC','',reference)
if (reference=='GRCh38') {
	getTPM(runlane, 'hsapiens_gene_ensembl')
} else if (reference=='mm10') {
	getTPM(runlane, 'mmusculus_gene_ensembl')
} else {
	cat(paste('Skipping TPM computation, unrecognised reference: ',reference,'\n',sep=''))
}
