#!/usr/bin/Rscript

library(methods)
library(Seurat)
library(dplyr)
library(Matrix)

args = commandArgs(TRUE)
samplename = args[1]
topgenes = args[2]

exprs = Read10X(data.dir=paste(samplename,'/outs/filtered_gene_bc_matrices/GRCh38',sep=''))
exprs = CreateSeuratObject(raw.data = exprs, min.cells = 3, min.genes = 200)
exprsum = Matrix::rowSums(exprs@data)
exprsum = sort(exprsum, decreasing=TRUE)
genes = names(exprsum)[1:topgenes]
#make the names GTF-compatible
for (i in 1:length(genes))
{
	genes[i] = paste('"',genes[i],'";',sep='')
}
write.table(genes,'topgenes.txt',sep='\n',quote=FALSE,row.names=FALSE,col.names=FALSE)
