#!/usr/bin/env python

from __future__ import print_function

"""convert between database ids

THIS PROGRAM IS UNDER CONSTRUCTION AND CURRENTLY DOESN't DO WHAT IT CLAIMS
To DO.

Copyright:

    dbID_converter.py Convert IDs for same sequence across databases
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
"""

__authors__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__version__ = '0.0.0.4'

import argparse
from collections import defaultdict
import copy
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

    header, originalTable = index_table(table_file, 'md5')
    annotations = read_annotation(annotation_file)
    header = header.strip().split('\t')
    originalHeader = copy.deepcopy(header)
    originalNumberColumns = len(header)
    tables = [originalTable]
    files = [table_file]
    for file in conversion_files:
        db = file.split('.')[0]
        try:
            db = db.split(os.sep)[-1]
        except:
            pass
        header.append(db)
        firstLine, table = index_table(file, 0, header = False)
        tables.append(table)
        files.append(file)
    conversionHeader = ['md5', 'ID', 'Annotation']
    outTable = table_file + '.tmp'
    with open(outTable, 'w') as out_table:
       out_table.write('\t'.join(header) + '\n')
    outAnn = annotation_file + '.tmp'
    with open(outAnn, 'w') as out_ann:
        out_ann.write('md5\tannotation\n')
    uniqueIDs = unique_ids(tables)
    del tables
    for id in uniqueIDs:
        firstItem = True
        toWrite = {}
        toWrite[id] = []
        for file, index in zip(files, uniqueIDs[id]):
            if firstItem:
                if index == 'None':
                    toWrite[id] = ['None' for i in range(originalNumberColumns)]
                else:
                    with open(file, 'rU') as in_handle:
                        data = read_indexed_table(in_handle, originalHeader,\
                                                  index)
                        for item in originalHeader:
                            toWrite[id].append(data[item])
                    toWrite[id][0] = annotations[id]
                firstItem = False
            else:
                if index == 'None':
                    toWrite[id].append('None')
                else:
                    with open(file, 'rU') as in_handle:
                        data = read_indexed_table(in_handle, conversionHeader,\
                                                  index)
                        toWrite[id].append(data['ID'])
                        if toWrite[id][0] == 'None':
                            # Remove database at end of line
                            newAnn = data['Annotation'].split(' ')[:-1]
                            toWrite[id][0] = ' '.join(newAnn)
        write_tables(outTable, outAnn, '', toWrite)
    os.remove(table_file)
    os.remove(annotation_file)
    os.rename(outTable, table_file)
    os.rename(outAnn, annotation_file)

def create_table(conversion_files, table_file, annotation_file):
    '''creates a table to allow for easy conversion between database IDs

    Input:

        conversion_files:
                A list of files from M5nr of the file name "*.md52id2func".
                Each file must have four columns containg the md52 checksum,
                database-specific id, annotation, and database for a given
                gene.

        out_table:
                File to conversion write table to

        out_annotation:
                File to write annotation table to

    Output:
            Two table files: one consisting of all the different database
            identifiers for a given gene and another containing the M5nr
            identifier and annotation for each gene.
    '''

    header = ['md5']
    db = conversion_files[0].split('.')[0]
    try:
        db = db.split(os.sep)[-1]
    except:
        pass
    header.append(db)
    headerStr = '\t'.join(header) + '\n'
    with open(table_file, 'w') as table_handle:
        table_handle.write(headerStr)
        with open(annotation_file, 'w') as ann_handle:
            ann_handle.write('md5\tannotation\n')
            with open(conversion_files[0], 'rU') as in_handle:
                for line in conversion_file_iter(in_handle):
                    tableOut = '{}\t{}\n'.format(line['md5'], line['dbID'])
                    annOut = '{}\t{}\n'.format(line['md5'], line['annotation'])
                    table_handle.write(tableOut)
                    ann_handle.write(annOut)
    conversionFiles = conversion_files[1:]
    append_table(conversionFiles, table_file, annotation_file)

def create_sub_table(out_table, databases, table_file, annotation_file):
    '''creates a new conversion table from existing table

    Input:

        out_table:
                The new table to write

        databases:
                The list of database names to write to the sub-table

        table_file:
                The table file to produce a sub-table from

        annotation_file:
                The annotation file of the original table file

    Output:
            A new table and annotation file containing only the
            specified databases
    '''

    firstLine, table = index_table(table_file, 'md5')
    annotations = read_annotation(annotation_file)
    tableDatabases = firstLine.strip().split('\t')
    for database in databases:
        if not database in tableDatabases:
            print('{} not in {}'.format(database, table_file))
            print('The databases in {} follow:'.format(table_file))
            print('\n'.join(tableDatabases))
            sys.exit(1)
    with open(out_table, 'w') as table_handle:
        table_handle.write('md5' + '\t' + '\t'.join(databases)+'\n')
    outAnnotation = out_table + '.ann'
    with open(outAnnotation, 'w') as annotation_handle:
        annotation_handle.write('md5\tannotation\n')
    with open(table_file, 'rU') as in_handle:
        for id in table:
            allNone = True
            index = table[id]
            data = read_indexed_table(in_handle, tableDatabases, index)
            for database in databases:
                if data[database] != 'None':
                    allNone = False
                    break
            if not allNone:
                outDict = {}
                outDict[id] = ['None' for i in range(len(databases))]
                for database in enumerate(databases):
                    outDict[id][database[0]] = data[database[1]]
                outDict[id] = [annotations[id]] + outDict[id]
                write_tables(out_table, outAnnotation, '', outDict)

def convert_column_ids(in_file, out_file, column_number, table_file,\
                        original_db, target_db, delimiter = '\t',\
                        log_file = None):
    '''stores a file and changes the IDs of a column within a file

    Input:

        in_file:
                CSV file containing IDs to be converted, assumes a header

        out_file:
                File to write to

        column_number:
                The column containing the IDs to be converted

        original_db:
                The database ID number contained in the in_file

        target_db:
                The database ID to convert to

        delimiter:
                The character delimiting the columns in in_file,
                default is "\\t"

        log_file:
                Log file to write errors to. If None, prints errors
                to stdin
    '''

    firstLine, table = index_table(table_file, original_db)
    with open(table_file, 'rU') as table_handle:
        header = table_handle.readline().strip().split('\t')
        with open(in_file, 'rU') as in_handle:
            with open(out_file, 'w') as out_handle:
                # Write header from in_handle
                out_handle.write(in_handle.readline())
                for line in in_handle:
                    lineSplit = line.strip().split(delimiter)
                    dbID = lineSplit[column_number]
                    try:
                        targetID = 'None'
                        assert dbID in table
                        try:
                            index = table[dbID]
                            data = read_indexed_table(table_handle, header,\
                                                      index)
                            targetID = data[target_db]
                        except KeyError:
                            message = '{} is not in {}'\
                                      .format(target_db, table_file)
                            if log_file:
                                with open(log_file, 'a') as out_handle:
                                    out_handle.write(message)
                            else:
                                print(message)
                            sys.exit(1)
                        if targetID == 'None':
                            raise AssertionError
                        lineSplit[column_number] = targetID
                        toWrite = delimiter.join(lineSplit) + '\n'
                        out_handle.write(toWrite)
                    except AssertionError:
                        message = '{} has no equivalent ID in {}'\
                                  .format(dbID, target_db)
                        if log_file:
                            with open(log_file, 'a') as out_handle:
                                out_handle.write(message)
                        else:
                            print(message)

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

def format_table_tools_args(args):
    '''Formats args to fit generalized conversion table related functions'''

    if args.create:
        create_table(args.create, args.conversionTable, args.annotationTable)
    elif args.append:
       append_table(args.append, args.conversionTable, args.annotationTable)
    elif args.create_sub_table:
        outTable = args.create_sub_table[0]
        databases = args.create_sub_table[1:]
        create_sub_table(outTable, databases, args.conversionTable,\
                         args.annotationTable)

def format_convert_id_file_args(args):
    '''Formats args to fit generalized file id conversion function'''

    column_number = int(args.column) - 1
    convert_column_ids(args.inFile, args.outFile, column_number,\
                       args.tableFile, args.originalDB, args.targetDB,\
                       delimiter = args.delimiter,\
                       log_file = args.log_file)

def index_table(table_file, key_db, header = True):
    '''reads in the table file as a nested dictionary and returns the dictionary

    Input:

        table_file:
                The table file to be indexed as a dictionary

        key_db:
                Which database values to use as the key

        header:
                Whether or not file has a header, default = True

    Output:
            A tuple containing the table header and a dictionary.
            The dictionary contains the database IDs for each ID
            in key_db as keys and the values are the location
            of that ID in the table file.
    '''

    table = {}
    with open(table_file, 'rU') as in_handle:
        if header:
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
        else:
            keyColumn = int(key_db)
            firstLine = ''
        lastPosition = in_handle.tell()
        while True:
            line = in_handle.readline()
            if not line:
                break
	    lineSplit = line.strip().split('\t')
            dbKey = lineSplit[keyColumn]
            if dbKey != 'None' and dbKey not in table:
                table[dbKey] = lastPosition
            lastPosition = in_handle.tell()
    return (firstLine, table)

def read_annotation(annotation_file):
    '''reads in the annotation_file as a dictionary'''

    annotations = {}
    with open(annotation_file, 'rU') as in_handle:
        for line in in_handle:
            lineSplit = line.strip().split('\t')
            annotations[lineSplit[0]] = lineSplit[1]
    return annotations

def read_indexed_table(table_handle, header, index):
    '''reads line of indexed table handle

    Input:

        table_handle:
                File handle of the indexed table file

        header:
                The header of the table file

        index:
                The index of the line to read

    Output:

            A dictionary where each line is a database in the header
            and each value is that databases ID
    '''

    table_handle.seek(index)
    line = table_handle.readline()
    lineSplit = line.strip().split('\t')
    # Combines all elements in the line beyond the length of the header
    # into a single element at the end of the list lineSplit
    if len(lineSplit) > len(header):
        toStr = ' '.join(i for i in lineSplit[len(header) - 1:])
        lineSplit = lineSplit[:len(header) - 1]
        lineSplit.append(toStr)
    lineDict = {}
    for key, value in zip(header, lineSplit):
       lineDict[key] = value
    return lineDict

def unique_ids(tables):
    '''yields unique IDs across tables

    Input:

        tables:
                A list of dictionaries

    Output:
        Returns the a dictionary containing only unique IDs
        across all dictionaries in tables. The keys are the
        unique IDs and the values are lists containing the value
        for that ID in each dictionary. If the dictionary does not
        containg the ID, the value at that position in the list
        is "None".
    '''

    uniqueIDs = {}
    for table in tables:
       for id in table:
           if id not in uniqueIDs:
               uniqueIDs[id] = []
    for id in uniqueIDs:
        for table in tables:
            if id in table:
                uniqueIDs[id].append(table[id])
            else:
                uniqueIDs[id].append('None')
    return uniqueIDs

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

    with open(out_table, 'a') as conversion_handle:
        if header_str:
            conversion_handle.write(header_str)
        with open(out_annotation, 'a') as annotation_handle:
            if header_str:
                annHeader = 'md5\tannotation\n'
                annotation_handle.write(annHeader)
            for id in identifiers:
                annotation = identifiers[id][0]
                dbIDs = id + '\t' + '\t'.join(identifiers[id][1:]) + '\n'
                annotationOutput = id + '\t' + annotation + '\n'
                conversion_handle.write(dbIDs)
                annotation_handle.write(annotationOutput)

if __name__ == '__main__':
    parentParser = MyParser(description = __doc__,
                            formatter_class = argparse.\
                            RawDescriptionHelpFormatter)
    parentParser.add_argument('--version',
                              action = 'version',
                              version = '0.0.0.4')
    subParsers = parentParser.add_subparsers(title = 'sub-commands')

    tableToolsParser = subParsers.add_parser('tableTools',
                                               help = 'creates conversion '+\
                                                      'table between databases')
    tableToolsParser.add_argument('conversionTable',
                                    help = 'conversion table to be created '+\
                                           'or appended to')
    tableToolsParser.add_argument('annotationTable',
                                    help = 'annotation table for conversion table')
    editMode = tableToolsParser.add_mutually_exclusive_group(required = True)
    editMode.add_argument('--create',
                          nargs = '+',
                          help = 'create new conversion table from '+\
                                 'conversion files')
    editMode.add_argument('--create_sub_table',
                          nargs = '+',
                          help = 'create new conversion table from '+\
                                 'databases in existing table, first '+\
                                 'item is table file name')
    editMode.add_argument('--append',
                          nargs = '+',
                          help = 'add database conversion data to table')
    tableToolsParser.set_defaults(func=format_table_tools_args)

    convertIdFile = subParsers.add_parser('convertIdFile',
                                      help = 'converts one database ID to '+\
                                             'another')
    convertIdFile.add_argument('inFile',
                           help = 'The character delimited file to convert '+\
                                  'database IDs in')
    convertIdFile.add_argument('outFile',
                           help = 'File to print results to')
    convertIdFile.add_argument('column',
                           help = 'the column of the file containing the '+\
                                  'database IDs to convert with 1-indexing')
    convertIdFile.add_argument('tableFile',
                           help = 'conversion table')
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
    convertIdFile.set_defaults(func=format_convert_id_file_args)
    args = parentParser.parse_args()
    args.func(args)
    sys.exit(0)
