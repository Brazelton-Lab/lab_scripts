#! /usr/bin/env python
# deletes entries from FASTA alignment that match sequence in given text file

"""
Copyright:

    FASTA_deleter_byname.py Delete FASTA entries match alignment
    Copyright (C) 2016  William Brazelton

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

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
