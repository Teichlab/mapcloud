library(SoupX)

for (feature in c('Gene','GeneFull'))
{
	#the raw data is the same for each feature choice
	raw = Seurat::Read10X(paste0('counts/',feature,'/raw'))
	for (filter in c('cr3','filtered'))
	{
		#skip if there's no scrublet file - that means doublet detection failed
		#(and something is off with the sample)
		if (!file.exists(paste0('counts/',feature,'/',filter,'/scrublet.csv')))
			next
		
		#read in the filtered cells and create soupx object
		filt = Seurat::Read10X(paste0('counts/',feature,'/',filter))
		sc = SoupChannel(raw,filt)

		#soupx needs clusters, so to get those we load up the scrublet results and get the primary cluster for each cell
		#and in turn, to get those, we split the cluster on the comma and take the first part
		scrublet = read.csv(paste0('counts/',feature,'/',filter,'/scrublet.csv'), row.names=1, stringsAsFactors=FALSE)
		get_cluster = function(x) strsplit(x,',')[[1]][1]
		scrublet[,'scrublet_leiden'] = sapply(scrublet[,'scrublet_leiden'], get_cluster)
		
		#now that we have clusters, we can resume the soupx'ing
		sc = setClusters(sc,scrublet[,'scrublet_leiden'])
		#this might sometimes actually fail. so try() it just in case
		message = try(sc <- autoEstCont(sc))
		if (class(message) == 'try-error')
			next
		out = adjustCounts(sc)
		#create a new folder and write the results there
		#not just the count matrix - also the soup profile and rho estimates
		dir.create(paste0('counts/',feature,'/',filter,'/soupx'))
		file.copy(paste0('counts/',feature,'/',filter,'/barcodes.tsv'), paste0('counts/',feature,'/',filter,'/soupx'))
		file.copy(paste0('counts/',feature,'/',filter,'/genes.tsv'), paste0('counts/',feature,'/',filter,'/soupx'))
		Matrix::writeMM(out, paste0('counts/',feature,'/',filter,'/soupx/matrix.mtx'))
		write.csv(sc$soupProfile, paste0('counts/',feature,'/',filter,'/soupx/soupProfile.csv'), quote=FALSE)
		write.csv(data.frame(sc$metaData$rho, row.names=colnames(out)), paste0('counts/',feature,'/',filter,'/soupx/rho.csv'), quote=FALSE)
	}
}