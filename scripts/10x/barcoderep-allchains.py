#!/usr/bin/env python3 

import sys
import argparse

parser = argparse.ArgumentParser(description="Create new entries for the secondary chains in the barocde_report file. This script Will rename the barcode.")
parser.add_argument('-b', help="", required=True, dest="barcode_report")
parser.add_argument('--frac', help="the abundance needs to be more the fraction of the primary chain", default=0.1, dest="frac")

args = parser.parse_args()

barcodeReport = args.barcode_report
frac = float(args.frac)

def GetChainType(v, j, c):
	s = ""
	if (c != "*" and c != "."):
		s = c 
	elif (j != "*" and j != "."):
		s = j
	elif (v != "*" and v != "."):
		s = v
	else:
		return 7
	
	if (s[0:3] == "IGH"):
		return 0
	elif (s[0:3] == "IGK"):
		return 1
	elif (s[0:3] == "IGL"):
		return 2
	elif (s[0:3] == "TRA"):
		return 3
	elif (s[0:3] == "TRB"):
		return 4
	elif (s[0:3] == "TRG"):
		return 5
	elif (s[0:3] == "TRD"):
		return 6
	else:
		return 7


def GetCellType(v, j, c, defaultType = "*"):
	chainType = GetChainType(v, j, c)
	if (chainType <= 2):
		return "B"
	elif (chainType <= 4):
		return "abT"
	elif (chainType <= 6):
		return "gdT"
	else:
		return defaultType

#AGAGTGGTCTATCCTA-1	abT	TRBV6-6*01,TRBD2*01,TRBJ2-5*01,*,TGTGCCAGTCTACTTGGGGGGACCCAGTACTTC,CASLLGGTQYF,3.61,AGAGTGGTCTATCCTA-1_48281,100.00,0	TRAV10*01,*,TRAJ34*01,*,TGTGTGGTGAGCGCCCGCACCGACAAGCTCATCTTT,CVVSARTDKLIF,1.00,AGAGTGGTCTATCCTA-1_65563,100.00,0	*	*
fp = open(barcodeReport)
for line in fp:
	#a master counter of all chains for the spot
	k = 0
	if (line[0] == "#"):
		print(line.rstrip())
		continue

	cols = line.rstrip().split() 
	barcode = cols[0]

	# Output the primary chain
	outputCols = cols[:]
	#split chains into separate entries if there are two primaries
	if (outputCols[2] != "*") and (outputCols[3] != "*"):
		for ind in [3,2]:
			outputColsHold = outputCols.copy()
			outputColsHold[ind] = "*"
			outputColsHold[0] = barcode + "_" + str(k)
			print("\t".join(outputColsHold))
			k += 1
	else:
		outputCols = cols[:]
		outputCols[0] = barcode + "_" + str(k)
		print("\t".join(outputCols))
		k += 1
	
	#output all secondary chains
	for chain in [1,2]:
		# Expand the secondary chain
		secondaryEntry = cols[3 + chain]
		if (cols[1 + chain] == "*" or secondaryEntry == "*"):
			continue
		primaryAbund = float(cols[1 + chain].split(',')[6])
		secondaryCols = secondaryEntry.split(";")
		for i in range(2,len(outputCols)):
			outputCols[i] = "*"

		for c in secondaryCols:
			outputCols[0] = barcode + "_" + str(k)
			outputCols[chain + 1] = c
			subCols = c.split(',')
			abund = float(subCols[6])
			outputCols[1] = GetCellType(subCols[0], subCols[2], subCols[3], cols[1])
			if (abund < primaryAbund * frac):
				continue
			print("\t".join(outputCols))
			k += 1
fp.close()