import sys
infilename = sys.argv[1]
outfilename = infilename + '.transposed'

with open(infilename) as infile:
	matrix = []
	for line in infile:
		items = line.split('\t')
		L = []
		for item in items:
			item = item.strip('\n')
			L.append(item)
		matrix.append(L)

#print matrix	
transposed = map(None,*matrix)
#print transposed

with open(outfilename, 'a') as outfile:
	for line in transposed:
		for item in line:
			if item == None: pass
			else: 
				outfile.write(item)
				outfile.write('\t')
		outfile.write('\n')	