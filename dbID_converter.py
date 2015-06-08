#!/usr/bin/env python

from __future__ import print_function

'''convert between database ids

THIS PROGRAM IS UNDER CONSTRUCTION AND CURRENTLY DOESN't DO WHAT IT CLAIMS
To DO.
'''

__authors__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__version__ = '0.0.0.3'

import argparse
from collections import defaultdict
import math
import os
import sys

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: ' + message)
        self.print_help()
        sys.exit(2)

def append_table(conversion_files, table_file, annotation_file):
    '''creates a table to allow for easy conversion between database IDs

    Input:

        conversion_files:
                A list of files from M5nr of the file name "*.md52id2func".
                Each file must have three columns containg the md52 checksum,
                database-specific id, and annotation on the gene respecitvely.

        table_file:
                The conversion table file to append the converesion entries to

        annotation_file:
                The annotation file to append annotations to

    Output:
            Appends conversion entries to conversion table and annotations
            to the annotation file from the conversion files
    '''

    header, tableLines = read_table(table_file, 'md5')
    header = header.strip().split('\t')
    originalColumnNumber = len(header) - 1
    table = {}
    for key in tableLines:
        items = []
        for item in header:
            if not item == 'md5':
                items.append(tableLines[key][item])
        table[key] = items
    annotations = read_annotation(annotation_file)
    for file in conversion_files:
        db = file.split('.')[0]
        try:
            db = db.split(os.sep)[-1]
        except:
            pass
        with open(file, 'rU') as in_handle:
            conversionLines = {}
            for line in conversion_file_iter(in_handle):
                conversionLines[line['md5']] = (line['dbID'],\
                                                line['annotation'])
            header.append(db)
            for key in table:
                if key in annotations:
                    table[key] = [annotations[key]] + table[key]
                else:
                    table[key] = ['None'] + table[key]
                if key in conversionLines:
                    table[key].append(conversionLines[key][0])
                else:
                    table[key].append('None')
            for key in conversionLines:
                if key not in table:
                    table[key] = [conversionLines[key][1]] +\
                                 ['None' for i in range(originalColumnNumber)]\
                                 + [conversionLines[key][0]]
        header = '\t'.join(header) + '\n'
        write_tables(table_file, annotation_file, header, table)

def create_table(conversion_files, out_table, out_annotation):
    '''creates a table to allow for easy conversion between database IDs

    Input:

        conversion_files:
                A list of files from M5nr of the file name "*.md52id2func".
                Each file must have three columns containg the md52 checksum,
                database-specific id, and annotation on the gene respecitvely.

        out_table:
                File to conversion write table to

        out_annotation:
                File to write annotation table to

    Output:
            Two table files: one consisting of all the different database
            identifiers for a given gene and another containing the M5nr
            identifier and annotation for each gene.

    Note: This function can be very RAM intensive.
    '''

    header = ['md5']
    annotation = ''
    numberOfFiles = len(conversion_files)
    identifiers = defaultdict(lambda: [line['annotation']] + \
                              ['None' for i in range(numberOfFiles)])
    for file in enumerate(conversion_files):
        db = file[1].split('.')[0]
        try:
            db = db.split(os.sep)[-1]
        except:
            pass
        header.append(db)
        with open(file[1], 'rU') as in_handle:
            for line in conversion_file_iter(in_handle):
                identifiers[line['md5']][file[0] + 1] = line['dbID']
    headerStr = '\t'.join(header) + '\n'
    write_tables(out_table, out_annotation, headerStr, identifiers)

def conversion_file_iter(handle):
    '''iterates over each line of a conversion file formatted as follows

    M5nr ID	Database ID	Annotations	Database

    Each line is returned as a dictionary consisting fo the first three
    columns with the respective keys: md5, dbID, annotation
    '''

    for line in handle:
        lineData = {}
        # [:-1] is to leave out the database name from M5nr
        cols = line.strip().split('\t')[:-1]
        lineData['md5'] = cols[0]
        lineData['dbID'] = cols[1]
        lineData['annotation'] = ' - '.join(cols[2:])
        yield lineData

def format_create_table_args(args):
    '''Formats args to fit generalized conversion table related functions'''
    
    if args.create:
        create_table(args.create, args.conversionTable, args.annotationTable)
    elif args.append:
       append_table(args.append, args.conversionTable, args.annotationTable)
    
def read_annotation(annotation_file):
    '''reads in the annotation_file as a dictionary'''

    annotations = {}
    with open(annotation_file, 'rU') as in_handle:
        for line in in_handle:
            lineSplit = line.strip().split('\t')
            annotations[lineSplit[0]] = lineSplit[1]
    return annotations

def read_table(table_file, key_db):
    '''reads in the table file as a nested dictionary and returns the dictionary

    Input:

        table_file:
                The table file to be read as a dictionary

        key_db:
                Which database values to use as the key 

    Output:
            A tuple containing the table header and a dictionary.
            The dictionary where the each key is an entry in the database
            specified by key_db. Each value is a dictionary where each
            key is a different database and each value is the equivalent
            database ID to the parent dictionary database ID.

    '''

    table = {}
    with open(table_file, 'rU') as in_handle:
        firstLine = in_handle.readline()
        splitFirstLine = firstLine.strip().split('\t')
        if key_db in splitFirstLine:
            keyColumn = splitFirstLine.index(key_db)
        else:
            print('{0} not in {1}'.format(key_db, table_file))
            print('Databases in {0} follow:'.format(table_file))
            databases = '\n'.join(splitFirstLine)
            print(databases)
            sys.exit(1)
        for line in in_handle:
            temp = {}
	    lineSplit = line.strip().split('\t')
            dbKey = lineSplit[keyColumn]
            if dbKey != 'None':
                for id in enumerate(lineSplit):
                    temp[splitFirstLine[id[0]]] = id[1]
                table[dbKey] = temp
    return (firstLine, table)

def write_tables(out_table, out_annotation, header_str, identifiers):
    '''writes a conversion table and annotation table

    Input:

        out_table:
                The converesion table to write

        out_annotation:
                The annotation table to write

        header_str:
                The header to write to the converion table

        identifiers:
                A dictionary of each M5nr checksum and each equivalent
                database ID

    Output:
            A database ID conversion table and a annotation table
    '''

    with open(out_table, 'w') as conversion_handle:
        conversion_handle.write(header_str)
        with open(out_annotation, 'w') as annotation_handle:
            annHeader = 'md5\tannotation\n'
            annotation_handle.write(annHeader)
            for key in identifiers:
                annotation = identifiers[key][0]
                dbIDs = key + '\t' + '\t'.join(identifiers[key][1:]) + '\n'
                annotationOutput = key + '\t' + annotation + '\n'
                conversion_handle.write(dbIDs)
                annotation_handle.write(annotationOutput)

if __name__ == '__main__':
    parentParser = MyParser(description = __doc__,
                            formatter_class = argparse.\
                            RawDescriptionHelpFormatter)
    parentParser.add_argument('--version',
                              action = 'version',
                              version = '0.0.0.2')
    subParsers = parentParser.add_subparsers(title = 'sub-commands')

    createTableParser = subParsers.add_parser('createTable',
                                               help = 'creates conversion '+\
                                                      'table between databases')
    createTableParser.add_argument('conversionTable',
                                    help = 'conversion table to be created '+\
                                           'or appended to')
    createTableParser.add_argument('annotationTable',
                                    help = 'annotation table for conversion table')
    editMode = createTableParser.add_mutually_exclusive_group(required = True)
    editMode.add_argument('--create',
                          nargs = '+',
                          help = 'create new conversion table')
    editMode.add_argument('--append',
                          nargs = '+',
                          help = 'add database conversion to table')
    createTableParser.set_defaults(func=format_create_table_args)

    convertIdFile = subParsers.add_parser('convertIdFile',
                                      help = 'converts one database ID to '+\
                                             'another')
    convertIdFile.add_argument('in_file',
                           help = 'The character delimited file to convert '+\
                                  'database IDs in')
    convertIdFile.add_argument('out_file',
                           help = 'File to print results to')
    convertIdFile.add_argument('column',
                           help = 'the column of the file containing the '+\
                                  'database IDs to convert with 0-indexing')
    convertIdFile.add_argument('originalDB',
                           help = 'the database type to convert from')
    convertIdFile.add_argument('targetDB',
                           help = 'the database type to convert to')
    convertIdFile.add_argument('--delimiter',
                           default = '\t',
                           help = 'the character delimiting the in_file, '+\
                                  'default is "\\t"')
    convertIdFile.add_argument('--log_file',
                           help = 'log file to redirect output to')
    args = parentParser.parse_args()
    args.func(args)
    sys.exit(0)
