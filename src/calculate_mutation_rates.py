#!/usr/bin/env python
"""
Precompute the background de novo mutation rate for every position in a given
BED file if specified or the entire genome
"""

import argparse
import sys
import pyfaidx
from ConfigParser import ConfigParser
from operator import le, lt
from OrionGlobals import *

def main(intervals, reference_genome, mutation_rate_by_trinucleotide,
         genome_mu_scaling_factor, missing_mu_rate, zero_based, output_fh):
    """output the mutation rate at each base in the intervals if specified,
    otherwise go through entire genome
    """
    for trinucleotide, mu in mutation_rate_by_trinucleotide.iteritems():
        mutation_rate_by_trinucleotide[trinucleotide] = mu * genome_mu_scaling_factor
    adjustment = 0 if zero_based else 1
    genome = pyfaidx.Fasta(reference_genome)
    if intervals:
        for chromosome, start, end in intervals:
            if 0 < start < end < CHROMOSOME_LENGTHS[chromosome]:
                # 0-based, end point not included
                sequence = genome[chromosome][start - 1:end + 1].seq
                for x in xrange(end - start):
                    trinucleotide = sequence[x:x + 3]
                    if trinucleotide in mutation_rate_by_trinucleotide:
                        mu = mutation_rate_by_trinucleotide[trinucleotide]
                    else:
                        mu = missing_mu_rate
                    # 1-based
                    output_fh.write("{CHROM}\t{POS}\t{mu}\n".format(
                        CHROM=chromosome, POS=start + x + adjustment, mu=mu))
            else:
                raise ValueError("invalid interval: {}:{}-{}".format(
                    chromosome, start, end))
    else:
        for chromosome in CHROMOSOME_LENGTHS.iterkeys():
            sequence = genome[chromosome][:].seq
            output_fh.write("{CHROM}\t{first}\t{missing_mu_rate}\n".format(
                CHROM=chromosome, first= adjustment,
                missing_mu_rate=missing_mu_rate))
            for x in xrange(len(sequence) - 2):
                trinucleotide = sequence[x:x + 3]
                if trinucleotide in mutation_rate_by_trinucleotide:
                    mu = mutation_rate_by_trinucleotide[trinucleotide]
                else:
                    mu = missing_mu_rate
                # 1-based
                output_fh.write("{CHROM}\t{POS}\t{mu}\n".format(
                    CHROM=chromosome, POS=x + 1 + adjusment, mu=mu))
            output_fh.write("{CHROM}\t{last}\t{missing_mu_rate}\n".format(
                CHROM=chromosome, last=len(sequence) - 1 + adjustment,
                missing_mu_rate=missing_mu_rate))
    genome.close()
    if output_fh is not sys.stdout:
        output_fh.close()

if __name__ == "__main__":
    cfg = get_cfg()
    genome_mu_rate = cfg.getfloat("sfs", "mu")
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter)
    parser.add_argument("-i", "--input_file", type=bed_file,
                        help="an optional BED file to get mutation rates for")
    parser.add_argument("-r", "--reference_genome", type=file_exists,
                        default=cfg.get("ref", "genome"),
                        help="the reference genome to use")
    parser.add_argument("-m", "--mutation_rate_matrix",
                        default=cfg.get("data", "mrt_matrix"),
                        type=mutation_rate_matrix,
                        help="the mutation rate by trinucleotide matrix")
    mu_group = parser.add_mutually_exclusive_group()
    mu_group.add_argument("--genome_mu_rate", type=non_negative_float,
                          default=genome_mu_rate,
                          help="the assumed overall rate of mutation per base "
                          "(mutation rate matrix values are scaled by the "
                          "specified value / the default value)")
    mu_group.add_argument("--genome_mu_scaling_factor", type=non_negative_float,
                          help="scale the mutation rate matrix by this factor")
    parser.add_argument("--missing_mu_rate", default=genome_mu_rate,
                        type=non_negative_float, help="the mu value to "
                        "default to for trinucleotides not in the matrx")
    parser.add_argument("-o", "--output_file", default=sys.stdout,
                        type=argparse.FileType("w"), help="the output file")
    parser.add_argument("-z", "--zero_based", default=False,
                        action="store_true", help="output zero-based scores")
    args = parser.parse_args()
    if args.genome_mu_scaling_factor:
        genome_mu_scaling_factor = args.genome_mu_scaling_factor
    else:
        genome_mu_scaling_factor = args.genome_mu_rate / genome_mu_rate,
    main(args.input_file, args.reference_genome, args.mutation_rate_matrix,
         genome_mu_scaling_factor, args.missing_mu_rate,
         args.zero_based, args.output_file)
