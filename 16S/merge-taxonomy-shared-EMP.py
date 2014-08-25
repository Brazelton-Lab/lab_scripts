# transfers taxonomy from .taxonomy file to each OTU in .shared.transposed file
# to make shared.transposed file:
# open .shared file in excel
# copy and paste>transpose all data 
# save as tab-delimited text
# open in TextWrangler and chance line breaks to Unix (LF)
# then run script:
# usage:
# python merge-taxonomy-shared.py myfile.taxonomy myfile.shared.transposed

import sys
taxfilename = sys.argv[1]
sharedfilename = sys.argv[2]
outfilename = sharedfilename + ".taxonomy"

D = {}
with open(taxfilename) as taxfile:
	for line in taxfile:
		cols = line.split('\t')
		otu = cols[0].split('_')
		#print cols
		otu = otu[-1]
		try: D[otu] = cols[1].replace(';','\t') 	
		except: pass
#print D
with open(sharedfilename) as sharedfile:
	for line in sharedfile:
		if "label" in line:
			with open(outfilename,'a') as outfile: outfile.write(line)
		elif "Group\t" in line: 
			with open(outfilename,'a') as outfile: outfile.write(line)
		elif "numOtus\t" in line: 
			with open(outfilename,'a') as outfile: outfile.write(line)
		else:
			with open(outfilename,'a') as outfile:
				outfile.write(line.strip('\n'))
				otu = line.split('\t')
				try: 
					otunum = otu[0]
					#print otunum
					outfile.write('\t' + D[otunum])			
				except: outfile.write(line)
