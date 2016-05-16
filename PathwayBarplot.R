#!/usr/bin/env Rscript
# Author: Alex Hyer
# Last Updated: 2016-05-11

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
parser$add_argument('--ratio', '-r',
                    action='store_true',
                    help='specifies to graph all abundances as ratios of total abundance')
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
  table = read.table(args$files[i], header=FALSE, fill=TRUE)
  split.file.name = strsplit(args$files[i], '/')
  file.parts = split.file.name[[1]]
  attr(table, 'file_name') = substr(tail(file.parts, n=1), 1, 13)
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
    if (!category=='UNDEFINED') {
      master.table[category.num, col.num] = 
        master.table[category.num, col.num] + abundance
    }
  }
}

# Remove seed from table
master.table = master.table[-1, -1, drop=FALSE]

# Take ratio or log of data if specified
if (args$ratio==TRUE) {
  for (c in 1:ncol(master.table)) {
    total = sum(master.table[,c])
    for (r in 1:nrow(master.table)) {
      master.table[r, c] = master.table[r, c] / total
    }
  }
} else if (!args$base==0) {
  for (r in 1:nrow(master.table)) {
    for (c in 1:ncol(master.table)) {
      if (!master.table[r, c]==0) {
        master.table[r, c] = log(master.table[r, c], base=args$base)
      }
    }
  }
}

# Write raw output table
table.out = paste(args$output, '.tsv', sep='')
write.table(master.table, file=table.out, sep='\t', na='0')

# Calculate image proprotions
scale = nrow(master.table)
image.height = 1000 + 8 * scale


# Create image file chosen by user
image.out = paste(args$output, '.', args$file_type, sep='')
if (args$file_type=='bmp'){
  bmp(filename=image.out, width=1000, height=image.height)
} else if (args$file_type=='jpeg') {
  jpeg(filename=image.out, width=1000, height=image.height)
} else if (args$file_type=='png') {
  png(filename=image.out, width=1000, height=image.height)
} else if (args$file_type=='tiff') {
  tiff(filename=image.out, width=1000, height=image.height)
}

# Set parameters used multiple times later
legend = rownames(master.table)
colors = rainbow(length(legend))

layout(rbind(1,2), heights=c(3,1))

# Format x-axis label
if (args$ratio==TRUE){
  x.label = 'Ratio of Total Abudnance per Sample'
} else if (!args$base==0) {
  x.label = paste('Abundance (log', args$base, '(RPK))', sep='')
} else {
  x.label = 'Abudnance (RPK)'
}

# Finally, plot the data
barplot(as.matrix(master.table), horiz=TRUE,
        col=colors,
        main='Metagenomic Pathway Profiles',
        xlab=x.label)

# Format legend
par(mar=c(0,0,0,0))
plot.new()
legend('center', legend, fill=colors, ncol=2)
