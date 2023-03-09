import argparse
import sys
import os

from braceexpand import braceexpand

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--samples', dest='samples', type=str, required=True, help='Sample IDs to process. Can be either space-delimited IDs, with brace expansion supported, or a  file with one sample ID per line. For CITE samples, please provide the GEX and CITE sample IDs joined with a "+"; brace expansion can be resolved on the IDs in parallel.')
	parser.add_argument('--cellranger', dest='cellranger', type=str, required=True, help='Path to cellranger release to use.')
	parser.add_argument('--command', dest='command', type=str, required=True, help='Cellranger command (count/vdj) to run.')
	parser.add_argument('--reference', dest='reference', type=str, required=True, help='Path to cellranger reference to use.')
	parser.add_argument('--library-type', dest='library_type', type=str, help='Optional. Library type to specify during data download. Some sample IDs have multiple libraries of different types associated with them. For CITE, provide two library types separated by "+", potentially skipping one by leaving the side empty.')
	parser.add_argument('--chemistry', dest='chemistry', type=str, help='Optional. 10X chemistry argument to pass to Cellranger.')
	parser.add_argument('--feature-ref', dest='feature_ref', type=str, help='CITE only. Path to feature reference file to use.')
	parser.add_argument('--chain', dest='chain', type=str, help='VDJ only. Optional. Chain to force in Cellranger. GD triggers dandelion post-processing.')
	parser.add_argument('--primers', dest='primers', type=str, help='VDJ only. Optional. Path to file with inner enrichment primers.')
	parser.add_argument('--dandelion', dest='dandelion', type=str, help='VDJ only. Path to dandelion singularity container if chain is set to GD.')
	parser.add_argument('--no-upload', dest='no_upload', action='store_true', help='Flag. If provided, will not upload to iRODS and just keep the results on the drive.')
	parser.add_argument('--no-bams', dest='no_bams', action='store_true', help='Flag. If provided, will delete the BAM files from the mapping output.')
	parser.add_argument('--dry', dest='dry', action='store_true', help='Flag. If provided, will just print the commands that will be called rather than running them.')
	args = parser.parse_args()
	#TODO: sanity check input
	#process sample input - is it a file?
	if os.path.isfile(args.samples):
		#it is a file. read its contents and trim off newlines
		with open(args.samples, 'r') as fid:
			lines = fid.readlines()
		args.samples = [i.rstrip() for i in lines]
	else:
		#it is not a file, it's a space-delimited list of IDs with brace expansion
		#identify all the lumps that need to be fed into brace expansion
		lumps = args.samples.split(' ')
		#start a fresh list in args.samples and append it with per-lump output
		args.samples = []
		for lump in lumps:
			#is this a CITE lump? judge by presence of "+"
			if "+" in lump:
				#brace expansion on GEX and CITE separately
				for gex, cite in zip(braceexpand(lump.split("+")[0]), braceexpand(lump.split("+")[1])):
					args.samples.append(gex+"+"+cite)
			else:
				#not CITE, one regular brace expansion will do
				for sample in braceexpand(lump):
					args.samples.append(sample)
	#ensure that the cellranger/reference paths lack any "~" surprises
	#(some cellrangers dislike that)
	args.cellranger = os.path.expanduser(args.cellranger)
	args.reference = os.path.expanduser(args.reference)
	#args.cellranger needs to point all the way to the binary
	#add the extra layer if it's just the folder
	if args.cellranger.split('/')[-1] != 'cellranger':
		args.cellranger = args.cellranger+'/cellranger'
	#if the feature reference is unset, make sure none of the samples are CITE
	if args.feature_ref is None:
		#CITE samples will have a "+" in their IDs
		if any("+" in sample for sample in args.samples):
			sys.stderr.write('CITE samples present, --feature-ref needs to be set.\n')
			sys.exit(1)
	#have VDJ chain set based on library_type, if applicable and necessary
	if args.library_type is not None:
		if args.chain is None:
			if args.library_type == 'TCR':
				args.chain = 'TR'
			elif args.library_type == 'BCR':
				args.chain = 'IG'
	#do we have the dandelion container, if we need to use it?
	if args.chain == 'GD':
		if args.dandelion is None:
			sys.stderr.write('GD samples present, --dandelion needs to be set.\n')
			sys.exit(1)
	#set up path to mapcloud/scripts, i.e. where this is
	#get the realpath to this file and then strip out the scripts.py at the end
	args.location = '/'.join(os.path.realpath(__file__).split('/')[:-1])
	return args

def runcommand(command, dry):
	if dry:
		print(command)
	else:
		#add pipefail just in case something goes wrong
		#(for this to work, need to specify bash like so)
		code = os.system('/bin/bash -c "set -eo pipefail"; '+command)
		#check that the command ran fine
		if code != 0:
			#we encountered an error
			sys.stderr.write('Error encountered while running: '+command+'\n')
			sys.exit(1)

def make_fastqs(sample, args, dest='fastq'):
	command = 'bash '+args.location+'/10x/utils/sample_fastq.sh '+sample
	#account for library_type if it's there
	if args.library_type is not None:
		#is this CITE? judge by dest value
		if dest == 'fastq_gex':
			#GEX. take first of pair of library types
			ltype = args.library_type.split('+')[0]
			#if present - this might be empty. if so, do nothing
			if ltype != '':
				command = command + ' "'+ltype+'"'
		elif dest == 'fastq_cite':
			#CITE. take second of pair of library types
			ltype = args.library_type.split('+')[1]
			#if present - this might be empty. if so, do nothing
			if ltype != '':
				command = command + ' "'+ltype+'"'
		else:
			#no need to do any splitting, just roll with the thing
			command = command + ' "'+args.library_type+'"' 
	runcommand(command, args.dry)
	#move to destination folder
	runcommand('mkdir '+dest+' && mv *.fastq.gz '+dest, args.dry)

def main():
	args = parse_args()
	#loop over the samples
	for sample in args.samples:
		#start preparing the command
		command = args.cellranger+' '+args.command
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
			command = command + ' --transcriptome='+args.reference
		elif args.command == 'vdj':
			command = command + ' --reference='+args.reference
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
		#set chemistry if provided
		if args.chemistry is not None:
			command = command + ' --chemistry='+args.chemistry
		#that's the whole command set up
		runcommand(command, args.dry)
		#CITE requires repositioning of folder to plus-joined sample names
		#and soupx post-processing
		if is_cite:
			runcommand('mv '+cr_cite_id+' '+sample, args.dry)
			runcommand('Rscript '+args.location+'/citeseq/soupx.R '+sample, args.dry)
		#remove temporary cellranger files
		runcommand("ls -d "+sample+"/* | grep -v 'outs' | xargs rm -r", args.dry)
		#potentially kick bams
		if args.no_bams:
			runcommand('rm '+sample+'/outs/*.bam', args.dry)
		#VDJ requires folder renaming, and possible GD postprocessing
		if args.command == 'vdj':
			if args.chain is not None:
				if args.chain == 'GD':
					#GD postprocessing, renames folder to dandelion
					runcommand('bash '+args.location+'/10x/dandelion.sh '+sample+' '+args.dandelion, args.dry)
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
				runcommand('bash '+args.location+'/irods-upload.sh 10X '+sample, args.dry)
			elif args.command == 'vdj':
				runcommand('bash '+args.location+'/irods-upload.sh 10X-VDJ '+sample, args.dry)
			#clean up mapping
			runcommand('rm -rf '+sample, args.dry)

if __name__ == "__main__":
	main()