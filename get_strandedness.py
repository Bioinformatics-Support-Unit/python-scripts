import pysam
import argparse

#chromosomal locations of gene in different annotations
#currently, Ensembl human annotation...
gapdh_location = {'grch37': ("12", 6643093, 6647481, '+'),
    'grch38': ("12", 6533927, 6538374, '+'),
    'hg19': ("chr12", 6643585, 6647537, '+'),}

def main(sam, genome):
    """Finds the strandedness of an RNA-Seq sample based on the reads
    that are mapped to the GAPDH gene in a range of supported annotation."""
    location = gapdh_location[genome]
    #assumes coordinate-sorted BAM file with index(!)
    samfile = pysam.AlignmentFile(sam, "rb")
    #fetch the GAPDH reads
    itersam = samfile.fetch(location[0], location[1], location[2])
    #counters for different read types
    f1, f2, r1, r2, other = 0, 0, 0, 0, 0
    for read in itersam:
        if read.flag & 208 == 80: #first in pair, reverse
            r1 += 1
        elif read.flag & 208 == 64: #first in pair, forward
            f1 += 1
        elif read.flag & 208 == 144: #second in pair, reverse
            r2 += 1
        elif read.flag & 208 == 128: #first in pair, forward
            f2 += 1
        else:
            other += 1 # should never get here (as all reads are mapped)
    total = float(r1+r2+f1+f2)
    if location[3] == '+':
        percent_firststrand = ((r1+f2)/total)*100
        percent_secondstrand = ((r2+f1)/total)*100
    elif location[3] == '-':
        percent_firststrand = ((f1+r2)/total)*100
        percent_secondstrand = ((f2+r1)/total)*100
    if percent_firststrand > 99:
        print("First Stranded")
    elif percent_secondstrand > 99:
        print("Second Stranded")
    else:
        print("Unstranded")
    print("Percent of reads supporting FR-firststrand (1=AS, 2=S): {0}".format(str(percent_firststrand)))
    print("Percent of reads supporting FR-secondstrand (1=S, 2=AS): {0}".format(str(percent_secondstrand)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Determine RNA-Seq Strandedness from an indexed BAM file')
    parser.add_argument('genome', type=str,
                       help='Reference genome of alignment, choose from: {0}'.format(', '.join(gapdh_location.keys())))
    parser.add_argument('bamfile', type=str,
                       help='BAM file of aligned RNA-Seq reads (must be coordinate sorted and indexed)')
    args = parser.parse_args()
    if args.genome not in gapdh_location.keys():
        print("Error: Unsupported genome")
        parser.print_help()
        exit()
    main(args.bamfile, args.genome)
