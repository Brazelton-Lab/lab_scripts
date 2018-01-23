#! /usr/bin/env python
#! /usr/bin/env perl

"""
This script is part of the Metagenomics Assembly workflow
It is designed to be used after assembly of reads with idba_ud in preparation for tetranucleotide binning with ESOM

Copyright:

    contigs-to-esom  tetranucleotide binning with ESOM

    Copyright (C) 2016  William Brazelton <comma-separated list of authors>

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

## load needed modules
use strict;
use warnings;
use Getopt::Long; # for accepting command line options use Cwd;
use Cwd;
Getopt::Long::Configure ("bundling");

## define global variables
my $smpl;
my $k = 40;
my $help = 0;
my $path = cwd();

## Set command line arguments
GetOptions( 
        'help|h' => \$help,
        'sample|s=s' => \$smpl,
        'kvalue|k=i' => \$k,
        'path|p=s' => \$path,
);

## Main Program
HelpMsg() if ($help);

die "Sample name required. See $0 --help for usage." unless ($smpl);

my $outfile = "all.$smpl.contig.cat";
my $logfile = "all.$smpl.contig.log";
my $contig = "contig-$k.fa";
my @files = <$path/*.idba.ud/$contig>;

print "Starting $0\n";
CatFile();
LogFile();
print "\ncontigs-to-esom.pl has finished processing contig files\noutput is in $outfile and notes on the processed files can be found in $logfile\n\n";

exit 0;

## Subroutines
sub HelpMsg {
die "contigs-to-esom.pl, version 1.0
contigs-to-esom.pl takes all contig files of a specified k value from partition 
groups assembled with idba_ud, concatenates them, and formats the concatenated file 
for use in esomWrapper.pl

usage: contigs-to-esom.pl -sample [sample_name] [option]...  

  -s, --sample [sample_name] :  name of the sample (required argument)
  -k, --k_value [k_value]    :  k value to use (default is 40). The matching 
                                contig-[k_value].fa file must exist for all assembled 
                                partition groups
  -p, --path [path_to_dirs]  :  path to the *.idba.ud directories (default is the 
                                current working directory)
  -h, --help                 :  display this help message
";
}
sub CatFile {
	open my $out, "> $outfile"
		or die $!;
	my $i = 0;
	print "\nFiles to be combined:\n";
	foreach my $file (@files) {
		print "\t$file\n";
		open my $in, "$file"
			or die $!;
                while (my $line = <$in>) {
                        if ($line =~ /contig/) {
                                $line =~ s/contig-\d+_\d+/$i++/ge;
                        }
                        print $out $line;
                }
                close $in;
        }
        close $out;
}
sub LogFile {
	open my $log, "> $logfile"
                or die $!;
	my $length = scalar(@files);
	my $logmessage = "## contigs-to-esom.pl ##\nFile Name: $outfile\nk value: $k (-k)\nFiles Processed: $length\nFiles combined:\n";
	print $log $logmessage;
        foreach my $file (@files) { 
                open my $in, "$file"
                        or die $!;
		my $firstline = <$in>;
		my @a = split('_', $firstline);
		$a[0] =~ m/>contig-\d+/;
		$a[0] =~ s/^.//;
		my $file_k_value = $a[0];
		print $log "\t$file\n\t\tcontig k value: $file_k_value\n";
		close $in;
	}
	close $log;
}
