#sponge up the translations
with open('outs/sampleInfo.txt','r') as fid:
	lines = fid.readlines()

#and now cook 'em up into a dictionary
trans = {}
for line in lines:
	line = line.strip().split('\t')
	#lose the weird thing up to the /
	line[0] = line[0].split('/')[1]
	#and now we can add to the dict
	trans[line[0]+'.cram'] = line[1]

with open('outs/snpCalls.vcf','r') as fid1:
	with open('outs/snpCallsNames.vcf','w') as fid2:
		for line in fid1:
			if line.startswith('#CHROM\tPOS'):
				#we have our header line. swap the .cram for the actual cell name
				line = line.rstrip().split('\t')
				#0:8 are pre-cell parts of the header
				for j in range(9,len(line)):
					line[j] = line[j].split('.')[0]+'-'+trans[line[j]]
				line = '\t'.join(line)+'\n'
			#write the line out to file, be it the reformatted header or any old content
			fid2.write(line)
