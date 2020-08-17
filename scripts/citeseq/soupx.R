#!/usr/bin/Rscript

library(SoupX)

#run with one positional argument - the citeseq combined sample identifier
#(that cellranger uses to create its output folder)
args = commandArgs(TRUE)
dataDir = file.path(args[1],'outs')
sc = load10X(dataDir, includeFeatures=c('Gene Expression', 'Antibody Capture'))
#this might sometimes actually fail. so try() it just in case
message = try(sc <- autoEstCont(sc))
if (class(message) != 'try-error')
{
	out = adjustCounts(sc)
	#create soupx subfolder within filtered_feature_bc_matrix
	#and write out a cellranger3-formatted mtx export there
	featurePath = file.path(dataDir,'filtered_feature_bc_matrix')
	soupxPath = file.path(featurePath,'soupx')
	dir.create(soupxPath)
	Matrix::writeMM(out, file.path(soupxPath,'matrix.mtx'))
	system(paste('gzip', file.path(soupxPath,'matrix.mtx')))
	for (fName in c('barcodes.tsv.gz', 'features.tsv.gz'))
	{
		system(paste('cp', file.path(featurePath,fName), file.path(soupxPath,fName)))
	}
}