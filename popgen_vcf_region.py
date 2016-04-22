from __future__ import print_function
import os
import sys
import argparse
import cyvcf
import pgstats as pg
import numpy as np
import datetime


def calc_rac(snp):
    rac_num = (snp.num_hom_ref * 2) + snp.num_het

    return rac_num


parser = argparse.ArgumentParser(description="Program to calculate population genetic statistics from a region "
                                             "in VCF file")
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
n_total = 2 * len(samples)  # assumes diploid genotypes.

#pop1_index = []
#pop2_index = []

#n1 = len(pop1_index)
#n2 = len(pop2_index)

indels = 0
extreme_depth = 0
low_call_rate = 0
repeat_sites = 0
spanning_deletion = 0
snp_not_S = 0
multiallelic_snp = 0
valid_sites = 0
failed_snp = 0
ref_N = 0

rac_list = []  # list to hold reference allele counts from SNP sites

site_counter = 0

# print("Chromosome", "Position", "%_sites_analysed", "Date time")


for site in vcf_infile.fetch(args.chrom, 0, chrom_length):

    site_counter += 1

    if site_counter % 100000 == 0:
        print(site.CHROM, site.POS, round((site_counter / float(chrom_length)) * 100, 2), datetime.datetime.now())

    if site.REF == 'N':
        #ref_N += 1
        continue

    if site.is_indel and site.aaf != 0.0:
        # indels += 1
        continue

    if len(site.ALT) >= 1 and site.ALT[-1] == '*':  # SNPs at spanning deletion
        #spanning_deletion += 1
        continue

    all_DP = [x['DP'] for x in site.samples]  # get the genotype depths

    if None in all_DP:  # Exclude sites where a genotype DP field is set to '.'
        #low_call_rate += 1
        continue

    mean_DP = np.nanmean(all_DP)

    if mean_DP < args.min_dp or mean_DP > args.max_dp: # depth filter
        #extreme_depth += 1
        continue

    if 0 in all_DP or site.call_rate < 1.0:  # only consider sites where all samples have coverage
        #low_call_rate += 1
        continue

    if site.FILTER == "REPEAT" or "REPEAT" in site.FILTER: # exclude sites in repeat regions
        #repeat_sites += 1
        continue

    if site.is_monomorphic:
        valid_sites += 1
        continue

    if site.is_snp:

        if len(site.ALT) > 1:
            multiallelic_snp += 1
            continue

        if site.aaf == 1.0 or site.aaf == 0.0: # only want SNPs polymorphic in our sample
            valid_sites += 1
            # if site.aaf == 1.0:
            #     snp_not_S += 1 # A SNP different from reference, but not segregating in the samples
            # continue

        if not args.unfiltered:
            if site.FILTER == 'PASS':
                rac = calc_rac(site)
                rac_list.append(rac)
                valid_sites += 1
                continue
            else:
                failed_snp += 1
                continue
        else:
            rac = calc_rac(site)
            rac_list.append(rac)
            valid_sites += 1
            continue
    else:
        problem_site = site.CHROM + '\t' + str(site.POS)
        error_message = 'Could not assign ' + problem_site + ' to site type'
        sys.exit(error_message)

S = len(rac_list)
pi = pg.pi(n_total, rac_list)

with open(args.outfile, 'a') as outfile:
    if os.stat(args.outfile).st_size == 0:
        print('chromosome', 'chromosome_length', 'valid_sites', 'n', 'S', 'pi', sep='\t', file=outfile)
    print(args.chrom, chrom_length, valid_sites, n_total, S, pi, sep='\t', file=outfile)
