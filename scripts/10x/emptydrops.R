#!/usr/bin/Rscript

library(methods)
library(DropletUtils)
library(BiocParallel)
library(Seurat)
library(Matrix)

args <- commandArgs(TRUE)
samplename <- args[1]

rawdata = Read10X(data.dir=list.dirs(paste(samplename,'/outs/raw_gene_bc_matrices',sep=''))[2])
cellranger.cm = Read10X(data.dir=list.dirs(paste(samplename,'/outs/filtered_gene_bc_matrices',sep=''))[2])

#do.while() emulation - repeat this until it "converges"
#i.e. there are no FALSES in sig that are also TRUES in out$LIMITED
#sig is a vector with TRUE for "real" cells and false otherwise
#while out$Limited is a vector specifying if extra simulation iterations are required for a cell
niters = 10000
#folder creation for output in case emptydrops errors out
dir.create(paste(samplename,'/outs/final-count-matrix',sep=''))
repeat
{
	#if emptydrops errors out, copy the cellranger count matrix, barf out the error to a text file and quit
	tryCatch({out = emptyDrops(rawdata, niters=niters, BPPARAM=MulticoreParam(detectCores()))},
		error = function(err) {write(as.character(err),paste(samplename,'/outs/final-count-matrix/EmptyDropsError.txt',sep=''))
			file.copy(list.dirs(paste(samplename,'/outs/raw_gene_bc_matrices',sep=''))[2],paste(samplename,'/outs/final-count-matrix',sep=''),recursive=TRUE)
			file.rename(list.dirs(paste(samplename,'/outs/final-count-matrix',sep=''))[2],paste(samplename,'/outs/final-count-matrix/cellranger',sep=''))
			quit()})
	sig = out$FDR <= 0.01 & !is.na(out$FDR)
	mask1 = !sig
	mask2 = out$Limited & !is.na(out$Limited)
	if (sum(mask1&mask2) == 0)
	{
		#we have converged. we can leave
		break
	}
	#we haven't converged. next iteration please
	niters = niters + 10000
}

#let's make the plot emptydrops makes too, just in case
png(paste(samplename,'/outs/final-count-matrix/emptydrops-plot.png',sep=''), width=12, height=12, units="in", pointsize=12, res=120)
plot(out$Total, out$LR, log="x", col=ifelse(sig, "red", "black"))
o <- order(out$Total)
lines(out$Total[o], out$Expected[o], col="dodgerblue")
legend("topleft", sprintf("%.2f", sum(sig)/length(sig) * 100), bty="n")
dev.off()

#okay, we now have what EmptyDrops believes to be the real cells, in sig
#let's get the explicit cell barcodes that get picked by each
ed.cells = rownames(out)[sig]
cr.cells = colnames(cellranger.cm)

#so we can create the union!
final.cells = union(ed.cells, cr.cells)
final.cm = rawdata[,colnames(rawdata)%in%final.cells]

#create file with info on which source each barcode is from
holder = rep('Both',dim(final.cm)[2])
emptydrops_mask = final.cells %in% ed.cells
cellranger_mask = final.cells %in% cr.cells
holder[cellranger_mask & !emptydrops_mask] = 'CellRanger'
holder[emptydrops_mask & !cellranger_mask] = 'EmptyDrops'
write.csv(data.frame(Barcode=final.cells,CellSource=holder),paste(samplename,'/outs/final-count-matrix/CellSource.csv',sep=''),quote=FALSE,row.names=FALSE)

#export the count matrices
dir.create(paste(samplename,'/outs/final-count-matrix/emptydrops',sep=''))
dir.create(paste(samplename,'/outs/final-count-matrix/cellranger',sep=''))
dir.create(paste(samplename,'/outs/final-count-matrix/union',sep=''))
emptydrops.cm = rawdata[,colnames(rawdata)%in%ed.cells]
writeMM(emptydrops.cm, paste(samplename,'/outs/final-count-matrix/emptydrops/matrix.mtx',sep=''))
writeMM(cellranger.cm, paste(samplename,'/outs/final-count-matrix/cellranger/matrix.mtx',sep=''))
writeMM(final.cm, paste(samplename,'/outs/final-count-matrix/union/matrix.mtx',sep=''))
write(paste(ed.cells,'-1',sep=''),paste(samplename,'/outs/final-count-matrix/emptydrops/barcodes.tsv',sep=''))
write(paste(cr.cells,'-1',sep=''),paste(samplename,'/outs/final-count-matrix/cellranger/barcodes.tsv',sep=''))
write(paste(final.cells,'-1',sep=''),paste(samplename,'/outs/final-count-matrix/union/barcodes.tsv',sep=''))
file.copy(paste(list.dirs(paste(samplename,'/outs/raw_gene_bc_matrices',sep=''))[2],'/genes.tsv',sep=''),paste(samplename,'/outs/final-count-matrix/emptydrops/genes.tsv',sep=''))
file.copy(paste(list.dirs(paste(samplename,'/outs/raw_gene_bc_matrices',sep=''))[2],'/genes.tsv',sep=''),paste(samplename,'/outs/final-count-matrix/cellranger/genes.tsv',sep=''))
file.copy(paste(list.dirs(paste(samplename,'/outs/raw_gene_bc_matrices',sep=''))[2],'/genes.tsv',sep=''),paste(samplename,'/outs/final-count-matrix/union/genes.tsv',sep=''))
