from __future__ import print_function
import os
import sys
import argparse
import cyvcf
import pgstats as pg
import numpy as np
import datetime
from iterartors import vcf_site_iterator


parser = argparse.ArgumentParser(description="Program to calculate population genetic statistics from a region in VCF "
                                             "file")
parser.add_argument('-i', '--vcf', required=True, dest='vcf_infile',
                    help="VCF input file (Needs to be filtered in some way. SNPs without PASS in filter field will be"
                         " excluded by default, but may be included if the -u flag is used (see below))")
parser.add_argument('-o', '--out', required=True, dest='outfile', help="Outfile to write results")
parser.add_argument('-d', '--min_mean_dp', required=False, dest='min_dp', type=int,
                    help="Minimum mean depth across genotypes")
parser.add_argument('-D', '--max_mean_dp', required=False, dest='max_dp', type=int,
                    help="Maximum mean depth across genotypes")
parser.add_argument('-L', '--chromosome', required=False, dest='chrom', type=str,
                    help="Chromsome to analyse")
parser.add_argument('-u', '--unifiltered_sites', required=False, action='store_true', dest='unfiltered',
                    help="Option to include sites without a PASS in the filter field of the VCF (default; False)")

# parser.add_argument('-x', '--pop_1', required=False, dest='pop1', help="Samples belonging to population 1")
# parser.add_argument('-y', '--pop_2', required=False, dest='pop2', help="Samples belonging to population 2")
# parser.add_argument('-r', '--region', required=False, dest='region', help="Region to calculate stats.
#  Format is chromosome:start-end or file listing regions (or sites) (1-based coordinate system).
#  If region is not given the stats will be calculated on the whole vcf file")

args = parser.parse_args()

min_dp = 1
if args.min_dp:
    min_dp = args.min_dp

max_dp = None
if args.max_dp:
    max_dp = args.max_dp

vcf_infile = cyvcf.Reader(open(args.vcf_infile, 'r'))

chrom_str = args.chrom + ','
header = vcf_infile.metadata
chrom_length = [int(header[x][:-1]) for x in header if chrom_str in x][0]

samples = vcf_infile.samples
n_total = 2 * len(samples)  # assumes diploid ge

#pop1_index = []
#pop2_index = []

#n1 = len(pop1_index)
#n2 = len(pop2_index)

vcf_region = vcf_infile.fetch(args.chrom, 0, chrom_length)

filtered = True
if args.unfiltered:
    filtered = False

sites = vcf_site_iterator(vcf_region, min_dp, max_dp, filtered)

rac_list = sites[0]
valid_sites = sites[1]

S = len(rac_list)
pi = pg.pi(n_total, rac_list)

with open(args.outfile, 'a') as outfile:
    if os.stat(args.outfile).st_size == 0:
        print('chromosome', 'chromosome_length', 'valid_sites', 'n', 'S', 'pi', sep='\t', file=outfile)
    print(args.chrom, chrom_length, valid_sites, n_total, S, pi, sep='\t', file=outfile)
