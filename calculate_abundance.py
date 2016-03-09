#!/usr/bin/env python
"""
This script takes an alignment file in SAM/BAM format and a
feature file in GFF format and calculates for each feature
the number of reads mapping to it.
"""
import sys, optparse, itertools, warnings, traceback, os.path
import HTSeq

class UnknownChrom( Exception ):
   pass
   
def invert_strand(iv):
   iv2 = iv.copy()
   if iv2.strand == "+":
      iv2.strand = "-"
   elif iv2.strand == "-":
      iv2.strand = "+"
   else:
      raise ValueError, "Illegal strand"
   return iv2

def scale_abundance_estimate(count, feature_len, read_count, method='rpk'):
    scale_factor = {
        'rpk': 1,
        'rpkm': 1000000,
        }
    count = (count / (feature_len / 1000)) * scale_factor[method]
    return count

def count_reads_in_features(sam_filename, gff_filename, samtype, order, 
      overlap_mode, feature_type, id_attribute, quiet, minaqual, normal):

    features = HTSeq.GenomicArrayOfSets("auto", False)
    counts = {}

    # Try to open samfile to fail early in case it is not there
    if sam_filename != "-":
        open(sam_filename).close()
      
    gff = HTSeq.GFF_Reader(gff_filename)
    i = 0
    try:
        for f in gff:
            if f.type == feature_type:
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    continue
                features[f.iv] += feature_id
                counts[f.attr[id_attribute]] = 0
            i += 1
            if i % 100000 == 0 and not quiet:
                sys.stderr.write( "%d GFF lines processed.\n" % i )
    except:
        sys.stderr.write( "Error occured when processing GFF file ({}):\n".format(gff.get_line_number_string()))
        raise

    if not quiet:
        sys.stderr.write( "{!s} GFF lines processed.\n".format(i))

    if len(counts) == 0:
        sys.stderr.write("Warning: No features of type '{}' found.\n".format(feature_type))

    if samtype == "sam":
        align_reader = HTSeq.SAM_Reader
    elif samtype == "bam":
        align_reader = HTSeq.BAM_Reader
    else:
        raise ValueError, "Unknown input format {} specified.".format(samtype)

    try:
        if sam_filename != "-":
            read_seq_file = align_reader(sam_filename)
            read_seq = read_seq_file
            first_read = iter(read_seq).next()
        else:
            read_seq_file = align_reader(sys.stdin)
            read_seq_iter = iter(read_seq_file)
            first_read = read_seq_iter.next()
            read_seq = itertools.chain([first_read], read_seq_iter)
        pe_mode = first_read.paired_end
    except:
        sys.stderr.write( "Error occured when reading beginning of SAM/BAM file.\n" )
        raise

    try:
        if pe_mode:
            if order == "name":
                read_seq = HTSeq.pair_SAM_alignments(read_seq)
            elif order == "pos":
                read_seq = HTSeq.pair_SAM_alignments_with_buffer(read_seq)
            else:
                raise ValueError, "Illegal order specified."
        empty = 0
        ambiguous = 0
        notaligned = 0
        lowqual = 0
        nonunique = 0
        i = 0   
        for r in read_seq:
            if i > 0 and i % 100000 == 0 and not quiet:
                sys.stderr.write("{!s} SAM alignment record{} processed.\n" % (i, "s" if not pe_mode else " pairs"))

            i += 1
            if not pe_mode:
                if not r.aligned:
                    notaligned += 1
                    continue
                try:
                    if r.optional_field("NH") > 1:
                        nonunique += 1
                        continue
                except KeyError:
                    pass
                if r.aQual < minaqual:
                    lowqual += 1
                    continue
                iv_seq = (invert_strand(co.ref_iv) for co in r.cigar if co.type == "M" and co.size > 0)
            else:
                if r[0] is not None and r[0].aligned:
                    iv_seq = (invert_strand(co.ref_iv) for co in r[0].cigar if co.type == "M" and co.size > 0)
                else:
                    iv_seq = tuple()
                if r[1] is not None and r[1].aligned:            
                    iv_seq = itertools.chain( iv_seq, 
                        (co.ref_iv for co in r[1].cigar if co.type == "M" and co.size > 0))
                else:
                    if (r[0] is None) or not (r[0].aligned):
                        notaligned += 1
                        continue         
                try:
                    if (r[0] is not None and r[0].optional_field("NH") > 1 ) or \
                            (r[1] is not None and r[1].optional_field("NH") > 1):
                        nonunique += 1
                        continue
                except KeyError:
                    pass
                if (r[0] and r[0].aQual < minaqual) or (r[1] and r[1].aQual < minaqual):
                    lowqual += 1
                    continue         
         
            try:
                if overlap_mode == "union":
                    fs = set()
                    for iv in iv_seq:
                         if iv.chrom not in features.chrom_vectors:
                             raise UnknownChrom
                         for iv2, fs2 in features[iv].steps():
                             fs = fs.union(fs2)
                elif overlap_mode == "intersection-strict" or overlap_mode == "intersection-nonempty":
                    fs = None
                    for iv in iv_seq:
                        if iv.chrom not in features.chrom_vectors:
                            raise UnknownChrom
                        for iv2, fs2 in features[ iv ].steps():
                            if len(fs2) > 0 or overlap_mode == "intersection-strict":
                                if fs is None:
                                    fs = fs2.copy()
                                else:
                                    fs = fs.intersection(fs2)
                else:
                    sys.exit("Illegal overlap mode.")
                if fs is None or len(fs) == 0:
                    empty += 1
                elif len( fs ) > 1:
                    ambiguous += 1
                else:
                    counts[list(fs)[0]] += 1
            except UnknownChrom:
                empty += 1

    except:
        sys.stderr.write("Error occured when processing SAM input ({}):\n".format(read_seq_file.get_line_number_string()))
        raise

    if not quiet:
        sys.stderr.write("{!s} SAM {} processed.\n" % (i, "alignments " if not pe_mode else "alignment pairs"))

    feature_abund = counts[fn] if formula == 'count' else normalize(counts[fn], feature_len, read_count, formula)

    for fn in sorted(counts.keys()):
        print "{}\t{!s}".format(fn, feature_abund)
    print "__no_feature\t{!s}".format(empty)
    print "__ambiguous\t{!s}".format(ambiguous)
    print "__too_low_aQual\t{!s}".format(lowqual)
    print "__not_aligned\t{!s}".format(notaligned)
    print "__alignment_not_unique\t{!s}".format(nonunique)

      
def main():
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('alignment_file', type=str, dest='alignment',
        help="file containing alignments. Must be in SAM or BAM format")

    parser.add_argument('feature_file', type=str, dest='feature',
        help="file containing feature annotations. Must be in GFF3 format")

    parser.add_argument('-f', '--format', metavar='FORMAT', dest='aformat',
        choices=['bam', 'sam'], default='bam',
        help="type of <alignment_file> data [default: bam]. Choice are sam or "
            "bam")

    parser.add_argument('-o', '--order', metavar='ORDER',
       choices=["position", "name"], default='position',
       help = "sorting order of <alignment_file> [default: position]. "
           "Paired-end sequencing data must be sorted either by position or "
           "by read name, and the sorting order must be specified. Ignored "
           "for single-end data. Choices are 'position' or 'name'" )

    parser.add_argument('-t', '--type', metavar='FEATURETYPE', dest='ftype',
        default= 'CDS', 
        help="feature type (3rd column in GFF file) to be used [default: CDS]. "
            "All features of other type are ignored")
         
    parser.add_argument('-a', '--attr', metavar='ATTRIBUTE',
        default="gene_id", 
        help="GFF attribute to be used as feature ID [default: gene_id]")

    parser.add_argument('-m', '--mode', metavar='MODE',
        choices=["union", "intersection-strict", "intersection-nonempty"], 
        default="union", 
        help="mode to handle reads overlapping more than one feature "
            "[default: union]")

    parser.add_argument('-n', '--norm', metavar='METHOD',
        choices=['count', 'rpkb'], default='count',
        help="normalization method to use [default: None]. Choices are "
            "reads per kilobase (rpk), reads per base (rpb), or none")

    parser.add_argument('-q', '--minqual', metavar='QUAL',
        type=int, default=10,
        help="skip all reads with alignment quality lower than the given "
         "minimum value [default: 10]")

    parser.add_argument( '-q', '--quiet', action='store_true',
        help="suppress progress report") # and warnings

    args = parser.parse_args()

    warnings.showwarning = my_showwarning
    try:
        count_reads_in_features(args.alignment, args.feature, args.aformat, args.order,
            args.mode, args.ftype, args.attr, args.quiet, args.minqual, args.norm)
    except:
        sys.stderr.write("  {}\n".format(sys.exc_info()[1]))
        sys.stderr.write("  [Exception type: {}, raised in {}:{}]\n"
            .format(sys.exc_info()[1].__class__.__name__, 
            os.path.basename(traceback.extract_tb(sys.exc_info()[2])[-1][0]), 
            traceback.extract_tb( sys.exc_info()[2] )[-1][1]))
        sys.exit(1)

def my_showwarning(message, category, filename, lineno=None, line=None):
   sys.stderr.write("Warning: {}\n".format(message))

if __name__ == "__main__":
   main()
   sys.exit(0)

