#!/usr/bin/env Rscript

suppressPackageStartupMessages(library('argparse'))

doc.string = paste('Make a barplot given HumanN2 abudnance files')

parser = ArgumentParser(description=doc.string)
parser$add_argument('--files', '-f',
                    required=TRUE,
                    nargs='+',
                    help='HumanN2 table files')
parser$add_argument('--output', '-o',
                    required=TRUE,
                    help='name of image file to output')
parser$add_argument('--level', '-l',
                    default=1,
                    type='integer',
                    help='Pathway level to analyze [Default: 1]')
parser$add_argument('--base', '-b',
                    default=0,
                    type='integer',
                    help='base to take log of for all data, 0=don\'t take log [Default: 0]')
parser$add_argument('--file_type', '-t',
                    default='png',
                    choices=c(
                      'bmp',
                      'jpeg',
                      'png',
                      'tiff'
                    ),
                    help='image file type to write [Default: png]')
args = parser$parse_args()
args$level = args$level + 1

# Create and populate list of table data frames from files
tables = list()
for (i in 1:length(args$files)) {
  table = read.table(args$files[i])
  split.file.name = strsplit(args$files[i], '/')
  file.parts = split.file.name[[1]]
  attr(table, 'file_name') = tail(file.parts, n=1)
  tables[[length(tables)+1]] = table
}

# Create and seed master table
master.table = data.frame(0)

# Go through each table and create column data in master.table
for (i in 1:length(tables)) {
  current.table = tables[[i]]
  master.table[,ncol(master.table)+1] = 0
  col.num = ncol(master.table)
  colnames(master.table)[col.num] = attr(current.table, 'file_name')
  categories = rownames(master.table)
  
  # Scan through current.table rows for abundance data
  for (row in 1:nrow(current.table)) {
    abundance = current.table[row, 1]
    category = as.vector(current.table[row, args$level])[1]
    category.num = match(category, categories)

    # If category doesn't yet exist in master.table, make it
    if (all(is.na(category.num))) {
      master.table[nrow(master.table)+1,] = 0
      rownames(master.table)[nrow(master.table)] = category
      categories = row.names(master.table)
      category.num = match(category, categories)
    }
    
    # Finally, add appropriate value
    master.table[category.num, col.num] = 
      master.table[category.num, col.num] + abundance
  }
}

# Remove seed from table
master.table = master.table[-1, -1]

# Take log of data if specified
if (!args$base==0) {
  for (r in 1:nrow(master.table)) {
    for (c in 1:ncol(master.table)) {
      if (!master.table[r, c]==0) {
        master.table[r, c] = log(master.table[r, c], base=args$base)
      }
    }
  }
}

# Write raw output table
table.out = paste(args$output, '.tbl', sep='')
write.table(master.table, file=table.out, sep='\t', na='0')

# Create image file chosen by user
image.out = paste(args$output, '.', args$file_type, sep='')
if (args$file_type=='bmp'){
  bmp(filename=image.out, width=1000, height=1000)
} else if (args$file_type=='jpeg') {
  jpeg(filename=image.out, width=1000, height=1000)
} else if (args$file_type=='png') {
  png(filename=image.out, width=1000, height=1000)
} else if (args$file_type=='tiff') {
  tiff(filename=image.out, width=1000, height=1000)
}

# Set parameters used multiple times later
legend = rownames(master.table)
colors = rainbow(length(legend))

# Detect and define needed graph boundraries
bottom.width = length(legend) + 3
par(xpd=TRUE, mar=par()$mar+c(bottom.width, 0, 0, 0))

# Format x-axis label
if (!args$base==0) {
  x.label = paste('Abundance (log', args$base, '(RPK))', sep='')
} else {
  x.label = 'Abudnance (RPK)'
}

# Finally plot the data
barplot(as.matrix(master.table), horiz=TRUE,
        col=colors,
        main='Metagenomic Pathway Profiles',
        xlab=x.label)
legend(-0.25, legend, fill=colors)