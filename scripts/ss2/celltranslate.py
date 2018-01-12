#sponge up the untranslated cells and make a list of them
with open('outs/cellhold.txt','r') as fid:
	cells = fid.readline()
cells = cells.split('\t')

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
	trans[line[0]] = line[1]

#translation o clock!
#note the shift by 1 - it starts with an ENSG that we want to leave in place
for i in range(1,len(cells)):
	cells[i] = trans[cells[i]]

#okay, that's a wrap. export!
with open('outs/cellnamehold.txt','w') as fid:
	fid.write('\t'.join(cells))
