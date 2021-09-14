import argparse
import os

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--samples', dest='samples', type=argparse.FileType('r'), required=True, help='File with sample IDs to process, one ID (or pair of IDs separated with "+" for CITE) per line.')
	parser.add_argument('--command', dest='command', type=str, required=True, help='Cellranger command (count/vdj) to run.')
	parser.add_argument('--reference', dest='reference', type=str, required=True, help='Reference to use. Must be present as a folder in ~/cellranger.')
	parser.add_argument('--version', dest='version', type=str, default='4.0.0', help='Cellranger version to use. Must be present as a folder in ~/cellranger. Default: 4.0.0')
	parser.add_argument('--library_type', dest='library_type', type=str, help='Library type to specify during data download. Some sample IDs have multiple different libraries of different types associated with them.')
	parser.add_argument('--feature_ref', dest='feature_ref', type=str, default='features.csv', help='CITE only. Feature reference file to use. Must be present in folder. Default: features.csv')
	parser.add_argument('--chain', dest='chain', type=str, help='VDJ only. Chain to force in Cellranger. GD triggers dandelion post-processing.')
	parser.add_argument('--primers', dest='primers', type=str, help='VDJ only. File with inner enrichment primers. Must be present in folder.')
	parser.add_argument('--no_upload', dest='no_upload', action='store_true', help='Flag. If provided, will not upload to iRODS and just keep the results on the drive.')
	parser.add_argument('--dry', dest='dry', action='store_true', help='Flag. If provided, will just print the commands that will be called rather than running them.')
	args = parser.parse_args()
	#TODO: sanity check input
	#have VDJ chain set based on library_type, if applicable and necessary
	if args.library_type is not None:
		if args.chain is None:
			if args.library_type == 'TCR':
				args.chain = 'TR'
			elif args.library_type == 'BCR':
				args.chain = 'IG'
	return args

def runcommand(command, dry):
	if dry:
		print(command)
	else:
		os.system(command)

def make_fastqs(sample, args, dest='fastq'):
	command = 'bash /mnt/mapcloud/scripts/10x/utils/sample_fastq.sh '+sample
	#account for library_type if it's there
	if args.library_type is not None:
		command = command + ' "'+args.library_type+'"' 
	runcommand(command, args.dry)
	#move to destination folder
	runcommand('mkdir '+dest+' && mv *.fastq.gz '+dest, args.dry)

def main():
	args = parse_args()
	#loop over the samples
	for sample in args.samples:
		sample = sample.rstrip()
		#start preparing the command
		command = '/home/ubuntu/cellranger/cellranger-'+args.version+'/cellranger '+args.command
		#download the fastqs
		#is this a CITE sample, so two IDs joined by "+"?
		is_cite = False
		if "+" in sample:
			is_cite = True
			#need two sets of fastq files downloaded
			gex = sample.split("+")[0]
			cite = sample.split("+")[1]
			make_fastqs(gex, args, dest='fastq_gex')
			make_fastqs(cite, args, dest='fastq_cite')
			#since we're here, may as well do the CITE-specific prep
			#set up libraries.csv
			runcommand('echo "fastqs,sample,library_type" > libraries.csv', args.dry)
			runcommand('echo $(pwd)"/fastq_gex,'+gex+',Gene Expression" >> libraries.csv', args.dry)
			runcommand('echo $(pwd)"/fastq_cite,'+cite+',Antibody Capture" >> libraries.csv', args.dry)
			#cellranger cannot handle pluses. create internal ID by joining gex and cite with a dash
			cr_cite_id = gex+'-'+cite
			#add to the command
			command = command + ' --id='+cr_cite_id+' --libraries=libraries.csv --feature-ref='+args.feature_ref
		else:
			#just a simple regular single sample
			make_fastqs(sample, args)
			#which passes fastqs like this, count or vdj
			command = command + ' --id='+sample+' --fastqs=fastq'
		#count/vdj pass their arguments a little differently
		if args.command == 'count':
			command = command + ' --transcriptome=/home/ubuntu/cellranger/'+args.reference
		elif args.command == 'vdj':
			command = command + ' --reference=/home/ubuntu/cellranger/'+args.reference
			#is there a chain we need to pass/account for?
			if args.chain is not None:
				#if it's GD, pass TR
				if args.chain == 'GD':
					command = command + ' --chain=TR'
				else:
					command = command + ' --chain='+args.chain
			#are there primers?
			if args.primers is not None:
				command = command + ' --inner-enrichment-primers='+args.primers
		#that's the whole command set up
		runcommand(command, args.dry)
		#CITE requires repositioning of folder to plus-joined sample names
		#and soupx post-processing
		if is_cite:
			runcommand('mv '+cr_cite_id+' '+sample, args.dry)
			runcommand('Rscript /mnt/mapcloud/scripts/citeseq/soupx.R '+sample, args.dry)
		#remove temporary cellranger files
		runcommand("ls -d "+sample+"/* | grep -v 'outs' | xargs rm -r", args.dry)
		#VDJ requires folder renaming, and possible GD postprocessing
		if args.command == 'vdj':
			if args.chain is not None:
				if args.chain == 'GD':
					#GD postprocessing, renames folder to dandelion
					runcommand('bash /mnt/mapcloud/scripts/10x/dandelion.sh '+sample, args.dry)
				else:
					#chain is not empty, so we forced a chain. move folder
					runcommand('mv '+sample+'/outs '+sample+'/'+args.chain.lower(), args.dry)
			else:
				#TODO: infer chain from output and move output
				1+1
		#clean up fastqs
		runcommand('rm -rf fastq && rm -rf fastq_gex && rm -rf fastq_cite && rm -f libraries.csv', args.dry)
		#upload time... if we're to do so
		if not args.no_upload:
			if args.command == 'count':
				runcommand('bash /mnt/mapcloud/scripts/irods-upload.sh 10X '+sample, args.dry)
			elif args.command == 'vdj':
				runcommand('bash /mnt/mapcloud/scripts/irods-upload.sh 10X-VDJ '+sample, args.dry)
			#clean up mapping
			runcommand('rm -rf '+sample, args.dry)

if __name__ == "__main__":
	main()