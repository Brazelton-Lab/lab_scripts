#!/usr/bin/env Rscript
# Author: Alex Hyer
# Last Updated: 11-10-2015

suppressPackageStartupMessages(library('argparse'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('phyloseq'))

doc.string = paste('Create a barplot from Mothur taxonomic data')

parser = ArgumentParser(description=doc.string)
parser$add_argument('--count_table', '-c',
                    required=TRUE,
                    help='Mothur count table')
parser$add_argument('--tsv', '-t',
                    required=TRUE,
                    help='taxonomy2tsv tsv file')
parser$add_argument('--output', '-o',
                    required=TRUE,
                    help='name of image file to output')
parser$add_argument('-m', '--merge_threshold',
                    type='double',
                    default=2.5,
                    help=paste('Defines the upper threshold for low-abundance', 
                               'sequences. Any sequences below this threshold',
                               'in all samples are merged into "Low-abundance"'))
parser$add_argument('--level', '-l',
                    default=2,
                    type='double',
                    help='taxonomic level to analyze [Default: 2]')
parser$add_argument('--file_type', '-f',
                    default='png',
                    choices=c(
                      'bmp',
                      'jpeg',
                      'png',
                      'tiff'
                    ),
                    help='image file type to write [Default: png]')
parser$add_argument('-n', '--no_legend',
                    action='store_true',
                    help='do not add a legend to the bar graph')
args = parser$parse_args()

# Read in tables and format them for down-stream processing
otus = read.table(args$count_table, header=TRUE, row.names=1)
otus = otu_table(otus, taxa_are_rows=TRUE)
otus = subset(otus, select=-c(1:1)) # to delete the "total" column from the mothur count_table
tax = read.table(args$tsv, header=FALSE, row.names=1)
tax = tax_table(as.matrix(tax))
merged = phyloseq(otus, tax)
args$level = paste('V', args$level, sep='')
if (!is.element(args$level, rank_names(merged))){
    cat(args$level, 'is not a taxonomic level in', args$tsv)
    cat('Available options are:', rank_names(merged))
    stop()
}
merged.props = transform_sample_counts(merged, function(x) 100 * x / sum(x))
merged.props.level = tax_glom(merged.props, args$level)

# Find universally low-abundance OTUs to merge in order to improve graph readability
otus.to.merge = c()
sample.otus.to.check = list()
    
# Find all low-abundance OTUs 
for (col in 1:length(colnames(otu_table(merged.props.level)))) {
    sample.otus.to.check[[length(sample.otus.to.check)+1]] = c('seed')
    for (row in 1:length(otu_table(merged.props.level)[, col])) {
        if (otu_table(merged.props.level)[row, col] < args$merge_threshold) {
            sample.otus.to.check[[length(sample.otus.to.check)]] =  c(sample.otus.to.check[[length(sample.otus.to.check)]], rownames(otu_table(merged.props.level))[row])
        }
    }
}

# See if OTUs are universally low in abundance
for (i in 1:length(sample.otus.to.check)) {
    sample = sample.otus.to.check[[i]] 
    for (n in 2:length(sample)) {
        otu = sample[n]
        add = TRUE
        for (c in 1:length(sample.otus.to.check)) {
            check.against = sample.otus.to.check[[c]]
            if (!is.element(otu, check.against)) {
                add = FALSE
            }
        }
        if (add==TRUE) {
            otus.to.merge = c(otus.to.merge, otu)
        }
    }   
}

# Merge low-abundance OTUs and re-label the merge as low-abundace
otus.to.merge = unique(otus.to.merge)
merged.props.level = invisible(merge_taxa(merged.props.level, otus.to.merge, archetype=otus.to.merge[1]))
index = which(rownames(tax_table(merged.props.level))==otus.to.merge[1])
tax_table(merged.props.level)[index,] = rep(c(paste('Low-abundance (<', args$merge_threshold, '%)', sep='')), 
                                        times=length(tax_table(merged.props.level)[index,]))

# Calculate image proprotions
scale = length(get_taxa_unique(merged.props.level, args$level))
if (args$no_legend==FALSE) {
    image.height = 1000 + 8 * scale
} else {
    image.height = 1000
}

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

# Set parameters and make the barplot
cbPalette = c('#999999', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7')
if (args$no_legend==TRUE) {
    legend.position = 'none'
} else {
    legend.position = 'bottom'
}
plot_bar(merged.props.level, y='Sample', x='Abundance (count)', args$level,
         title=paste('Taxonomy by Sample (', args$level, ')', sep='')) + 
         coord_flip() + ylab('Percent Sequences') + theme(legend.position=legend.position) +
         guides(fill=guide_legend(keywidth=1, keyheight=1, ncol=4, label.position="right",
                title.position="top", title.hjust=0.5)) + 
         scale_fill_manual(values=colorRampPalette(cbPalette)(scale))

invisible(dev.off())
