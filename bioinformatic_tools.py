'''Multiple tools related to bioinformatics.

Contains a variety of tools that perform functions often used in
bioinformatics. Each tool aims to accomplish a single common task
such as returning a list of desired files.

Copyright:

    bioinformatic_tools.py Misc. functions for use in bioinformatics
    Copyright (C) 2016  William Brazelton, Alex Hyer

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
'''

import os
import sys
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML
import datetime

#BCBio is not a common module, therefore importing it is optional.
try:
    from BCBio import GFF
except ImportError:
    pass

__version__ = '0.002.000'

def unique_blast_databases(databaseDirectory):
    '''Returns a list of unique databases in a given directory.

    Returns a list of database names intended for use in BLAST
    searches. Only databases ending in ".nhr" and ".00", if
    a database has multipe partitions, are selected, trimmed
    of their ".nhr" and ".00" extension, and returned as a list.
    '''
    uniqueDatabases = []
    for path, dirs, files in os.walk(databaseDirectory):
        for file in files:
            if file.endswith('.nhr'):
                removedExtension = file.replace('.nhr', '')
                parseRemainder = removedExtension.split('.')
                try:
                    indexTest = int(parseRemainder[-1])
                    isContig = False
                except ValueError:
                    indexTest = None
                    if parseRemainder[-1] == 'contigs':
                        isContig = True
                if indexTest == 0 or indexTest == None:
                    coreName = parseRemainder[0]
                    if isContig == True:
                        coreName += '.contigs'
                    uniqueDatabases.append(coreName)
    return uniqueDatabases

def files_in_folder(folderDirectory, fileTypes = 'fasta'):
    '''Returns a list of all files files in a folder.

    fileTypes = [list_of_desired_file_extensions]

    Returns a list of all files ending in any file extension
    specified when the function is called. If the fileType is
    not specified manually, the function defaults to searching
    for FASTA files ending in ".fasta", ".fa", ".fas", ".ffn",
    ".fna", ".faa", ".frn", and ".mpfa."
    '''
    filesInFolder = []
    if fileTypes == 'fasta':
        folderFileTypes = ['fasta', 'fa', 'fas', 'ffn', 'fna',\
                          'faa', 'frn', 'mpfa']
    else:
        folderFileTypes = fileTypes
    for path, dirs, files in os.walk(folderDirectory):
        for file in files:
            for eachType in folderFileTypes:
                if file.endswith(eachType):
                    filesInFolder.append(file)
        break
    return filesInFolder

def convert_file(in_file, out_file):
    '''Converts between bioinformatic file types.

    Detects file type, including GFF, of in_file and out_file and converts
    in_file -> out_file if possible. If the conversion is not possible a
    detailed error message is given.
    '''
    with open(in_file, 'r') as in_handle:
        with open(out_file, 'w') as out_handle:
            firstFileType = (in_file.rsplit('.', 1))[-1]
            secondFileType = (out_file.rsplit('.', 1))[-1]
            try:
                if secondFileType == 'gff':
                    GFF.write(SeqIO.parse(in_handle, firstFileType), out_handle)
                else:
                    SeqIO.convert(in_file, firstFileType, out_file,\
                                  secondFileType)
            except ValueError:
                acceptableFileTypeList = ['abif', 'ace', 'embl', 'fasta',\
                                          'fastq', 'fastq-sanger',\
                                          'fastq-solexa', 'fastq-illumina',\
                                          'genbank', 'gb', 'ig', 'igmt',\
                                          'nexus', 'phd', 'phlip', 'pir',\
                                          'seqxml', 'sff', 'sff-trim',\
                                          'stockholm', 'swiss', 'tab',\
                                          'qual', 'uniprot-xml']
                print('One of the file\'s format is invalid.')
                print('The following list is a list of all formats accpetable')
                print('for the file you wish to convert (file1). If your file1')
                print('type is in this list then the file you wish to convert')
                print('to (file2) is invalid. If your file2 type is in this')
                print('list but you are seeing this error then the file you')
                print('wish to write to is invalid. This function can also')
                print('write to GFF files.\n')
                for fileType in acceptableFileTypeList:
                    print(fileType)
                sys.exit()

def cigar_from_blast(query, match, sbjct, cigarAge = 'old'):
    '''Converts a BLAST alignment into a CIGAR line.

    cigarAge = 'new' or 'old'

    Reads a given query sequence, alignment sequence (match), and subject
    sequence given in the BLAST format and generates a CIGAR line from them.
    cigarAge 'old' assigns all alignments, matches or mismatches, as 'M' while
    cigarAge 'new' assigns all matching alignments as '=' and all mismatching
    alignments as 'X'.
    '''
    positionList = ['' for base in sbjct]
    #Find deletions
    position = 0
    for base in query:
        if base == '-':
            positionList[position] = 'D'
        position += 1
    #Find matches
    position = 0
    for base in match:
        if base == '+' or base == '|' or base.isalpha():
            if cigarAge == 'old':
                positionList[position] = 'M'
            elif cigarAge == 'new' and base != '+':
                positionList[position] = '='
        position += 1
    #Find insertions
    position = 0
    for base in sbjct:
        if base == '-':
            positionList[position] = 'I'
        position += 1
    position = 0
    #Find mismatches
    position = 0
    for place in positionList:
        if place == '':
            if cigarAge == 'old':
                positionList[position] = 'M'
            elif cigarAge == 'new':
                positionList[position] = 'X'
        position += 1
    #Concatenate list into CIGAR line
    cigarLine = ''
    previous = positionList[0]
    repeats = 0
    position = 0
    for alignment in positionList:
        if alignment == previous:
            repeats += 1
        else:
            cigarLine += str(repeats) + positionList[position-1]
            repeats = 1
            previous = alignment
        if position == len(positionList)-1:
            cigarLine += str(repeats) + positionList[position]
        position += 1
    return cigarLine


def blast_to_sam_and_bam(blastXMLFile, outputFile, cigarAge = 'old'):
    '''Converts BLAST XML and tabular output into sam and bam files.

    cigarAge = 'old' or 'new' (see cigar_from_blast docstring for details)
    outputFile = 'file name without .sam extension'

    CAUTION: This function is currently untested.

    WARNING: The SAM and BAM files output by this function are only
    viable with BLASTN searches.

    Takes new BLAST engine XML output (option -m 7 for BLASTN output)
    and converts them to a sam file, bam file, sorted bam file, and
    indexed bam file. Uses the cigar_from_blast function.
    '''
    with open(blastXMLFile, 'r') as in_file:
        fileHeader = '@HD\tVN:1.0'
        sequenceHeaders = []
        programHeaders = []
        alignmentLines = []
        results = NCBIXML.parse(in_file)
        for result in results:
            sequenceHeader = '@SQ\tSN:' + str(result.database) + '\tLN:' +\
                             str(result.num_letters_in_database)
            sequenceHeaders.append(sequenceHeader)
            programHeader = '@PG\tID:' + str(result.application) + '\tPN:' +\
                            str(result.application) + '\tVN:' +\
                            str(result.version)
            programHeaders.append(programHeader)
            for alignment in result.alignments:
                for hsp in alignment.hsps:
                    qname = str(alignment.title)
                    flag = '4'
                    rname = '*'
                    pos = str(hsp.sbjct_start)
                    mapq = '255'
                    cigar = cigar_from_blast(hsp.query, hsp.match,\
                                             hsp.sbjct, cigarAge)
                    rnext = '*'
                    pnext = '0'
                    tlen = str(hsp.sbjct_end-hsp.sbjct_start)
                    seq = str(hsp.query)
                    qual = '*'
                    alignmentLine = qname + '\t' + flag + '\t' + rname +\
                                  '\t' + pos + '\t' + mapq + '\t' + cigar\
                                  + '\t' + rnext + '\t' + pnext + '\t' + \
                                  tlen + '\t' + seq + '\t' + qual
                    alignmentLines.append(alignmentLine)
        with open(outputFile + '.sam', 'w') as out_file:
            out_file.write(fileHeader)
            for sequenceHeader in sequenceHeaders:
                out_file.write('\n' + sequenceHeader)
            for programHeader in programHeaders:
                out_file.write('\n' + programHeader)
            for alignmentLine in alignmentLines:
                out_file.write('\n' + alignmentLine)
    subprocess.call(['samtools', 'view', '-hSbo', outputFile + '.bam',\
                     outputFile + '.sam'])
    subprocess.call(['samtools', 'sort', outputFile + '.bam',\
                     outputFile + '.sorted'])
    subprocess.call(['samtools', 'index', outputFile + '.sorted.bam'])

def qualityCheck(out_log_name, file_type, *files):
    #Records sequence number and file size of files so user can check quality
    with open(out_log_name + '.log.txt', 'a') as out_handle:
        today = datetime.datetime.today()
        out_handle.write(today.strftime('Date of run: %m/%d/%Y'))
        out_handle.write(today.strftime('\nTime of run: %H:%M:%S'))
        command_line = ''
        for i in sys.argv:
            command_line += i
            command_line += ' '
        out_handle.write('\nCommand: ' + command_line)
        for file in files:
            with open(file, 'r') as in_handle:
                seq_number = 0
                for seq_record in SeqIO.parse(in_handle, file_type):
                    seq_number += 1
                file_size = os.path.getsize(file)
                out_handle.write('\n\nFile Name: ' + file)
                out_handle.write('\nNumber of sequences: ' + str(seq_number))
                out_handle.write('\nFile size: ' + str(file_size))

def fastaOrFastq(in_file, cat = False):
    #Determines if file is a FASTA or FASTQ file, can read concatenated files
    fastaFileTypes = ['fasta', 'fa', 'fas', 'ffn', 'fna',
                      'faa', 'frn', 'mpfa', 'rbp']
    fastqFileTypes = ['fastq', 'fq']
    fileTypesSet = [fastaFileTypes, fastqFileTypes]
    split_name = in_file.split('.')
    for fileTypeSet in fileTypesSet:
        for fileType in fileTypeSet:
            whereToCheck = -1
            if cat == True:
                whereToCheck = -2
            if split_name[whereToCheck] == fileType:
                return fileTypeSet[0]


