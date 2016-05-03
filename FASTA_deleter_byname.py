#! /usr/bin/env python
# deletes entries from FASTA alignment that match sequence in given text file

############################
# enter the name of the original .fasta file:
fasta_filename = 'all.filter.fasta.phylipok.trimmed.prunedrefs.nolongs'

# enter the name of the text file containing the sequences to be deleted - one per line
seqs_to_be_del = 'RAxML_classification.CloErysFullPruned.June2013.EPA1.longbranches.txt'
# must have unix linebreaks

# enter the name of the new .fasta file with entries deleted
outfile_name = 'all.filter.fasta.phylipok.trimmed.prunedrefs.nolongs2'
###############################

FASTA_file = open(fasta_filename).read()
FASTAs = FASTA_file.split('>')

del_list = []
textfile = open(seqs_to_be_del)
for line in textfile:
	if len(line) < 2: pass
	else:
		name = line.strip('\n')
		name = name.strip('>')
		del_list.append(name)
#print keep_list

outfile = open(outfile_name, 'a')

count = 0
for FASTA in FASTAs[1:]:
	FASTA_split = FASTA.split('\n')
	header = FASTA_split[0]
	header = header.split(' ')
	header = header[0]
	header = header.split('\t')
	header = header[0]
	#print header
	if header in del_list: pass
	else:
		outfile.write('>')
		outfile.write(FASTA)
		count = count + 1
print len(FASTAs)-1,
print 'sequences in original FASTA'		
print len(del_list),
print 'sequences in delete file'
print count,
print 'sequences in new FASTA'
		
textfile.close()				
outfile.close()
