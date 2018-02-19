#!/usr/bin/Rscript

library(methods)
library(Seurat)
library(dplyr)
library(Matrix)

args = commandArgs(TRUE)
topgenes = args[1]

exprs = read.csv('outs/countMatrixNames.txt',sep='\t',row.names=1)
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

