#!/usr/bin/Rscript

library(SoupX)

#run with one positional argument - the citeseq combined sample identifier
#(that cellranger uses to create its output folder)
args = commandArgs(TRUE)
dataDir = file.path(args[1],'outs')
#load just the CITE part
sc = load10X(dataDir, includeFeatures=c('Antibody Capture'))
#this might sometimes actually fail. so try() it just in case
message = try(sc <- autoEstCont(sc, soupQuantile=0.25, tfidfMin=0.2, forceAccept=TRUE))
if (class(message) != 'try-error')
{
	out = adjustCounts(sc)
	#create soupx-cite subfolder within filtered_feature_bc_matrix
	#and write out a cellranger3-formatted mtx export there
	featurePath = file.path(dataDir,'filtered_feature_bc_matrix')
	soupxPath = file.path(featurePath,'soupx-cite')
	dir.create(soupxPath)
	Matrix::writeMM(out, file.path(soupxPath,'matrix.mtx'))
	system(paste('gzip', file.path(soupxPath,'matrix.mtx')))
	system(paste('cp', file.path(featurePath,'barcodes.tsv.gz'), file.path(soupxPath,'barcodes.tsv.gz')))
	#the genes need a little postprocessing nudge - filter it to just the AB_ ones, then compress
	system(paste('zcat', file.path(featurePath,'features.tsv.gz'), '| grep "AB_" >', file.path(soupxPath,'features.tsv.gz')))
	system(paste('rm', file.path(soupxPath,'features.tsv.gz')))
	system(paste('gzip', file.path(soupxPath,'features.tsv')))
}