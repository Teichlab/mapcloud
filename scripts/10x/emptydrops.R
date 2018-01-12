#!/usr/bin/Rscript

library(methods)
library(Seurat)
library(EmptyDrops)

args <- commandArgs(TRUE)
samplename <- args[1]
reference <- args[2]

rawdata = Read10X(data.dir=paste(samplename,'/outs/raw_gene_bc_matrices/',reference,sep=''))
cellranger.cm = Read10X(data.dir=paste(samplename,'/outs/filtered_gene_bc_matrices/',reference,sep=''))

#do.while() emulation - repeat this until it "converges"
#i.e. there are no FALSES in sig that are also TRUES in out$LIMITED
#sig is a vector with TRUE for "real" cells and false otherwise
#while out$Limited is a vector specifying if extra simulation iterations are required for a cell
npts = 20000
repeat
{
	out = detectCells(rawdata, npts=npts, BPPARAM=MulticoreParam(28))
	sig = out$FDR <= 0.01 & !is.na(out$FDR)
	mask1 = !sig
	mask2 = out$Limited & !is.na(out$Limited)
	if (sum(mask1&mask2) == 0)
	{
		#we have converged. we can leave
		break
	}
	#we haven't converged. next iteration please
	npts = npts + 10000
}

#let's make the plot emptydrops makes too, just in case
dir.create(paste(samplename,'/outs/final-count-matrix',sep=''))
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
saveRDS(final.cm, paste(samplename,'/outs/final-count-matrix/final-union-count-matrix.RDS',sep=''))

#but since we have the EmptyDrops and CellRanger cells already, may as well export those too
emptydrops.cm = rawdata[,colnames(rawdata)%in%ed.cells]
saveRDS(emptydrops.cm, paste(samplename,'/outs/final-count-matrix/emptydrops-count-matrix.RDS',sep=''))
saveRDS(cellranger.cm, paste(samplename,'/outs/final-count-matrix/cellranger-count-matrix.RDS',sep=''))
