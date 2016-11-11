from __future__ import print_function, division
import sys
import argparse
from pysam import
import pgstats as pg
import numpy as np


parser = argparse.ArgumentParser(description="Program to calculate population genetic statistics from a region "
                                             "in VCF file")
parser.add_argument('-i', '--vcf', required=True, dest='vcf_infile',
                    help="VCF input file (Needs to be filtered)")
parser.add_argument('-o', '--out', required=True, dest='outfile', help="Outfile to write results")

# parser.add_argument('-x', '--pop_1', required=False, dest='pop1', help="Samples belonging to population 1")
# parser.add_argument('-y', '--pop_2', required=False, dest='pop2', help="Samples belonging to population 2")
# parser.add_argument('-b', '--bed', required=False, dest='region', help="Bed file specifying regions to calculate stats
# (e.g., bedfile specifying introns")

args = parser.parse_args()

vcf_infile = Variant(vcf_infile)


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

print("Chromosome", "Position", "sites_analysed", "Date time")

for site in vcf_infile:
        
    if site.CHROM not in chrom_list:
        continue

    site_counter += 1

    if site_counter % 100000 == 0:
        print(site.CHROM, site.POS, site_counter, datetime.datetime.now())

    if site.REF == 'N':
        ref_N += 1
        continue

    if site.is_indel:
        indels += 1
        continue

    if len(site.ALT) >= 1 and site.ALT[-1] == '*':  # SNPs at spanning deletion
        spanning_deletion += 1
        continue

    all_DP = [x['DP'] for x in site.samples]  # get the genotype depths

    if None in all_DP:  # Exclude sites where a genotype DP field is set to '.'
        low_call_rate += 1
        continue

    mean_DP = np.nanmean(all_DP)

    if mean_DP < args.min_dp or mean_DP > args.max_dp: # depth filter
        extreme_depth += 1
        continue

    if 0 in all_DP or site.call_rate < 1.0:  # only consider sites where all samples have coverage
        low_call_rate += 1
        continue

    if site.FILTER == "REPEAT" or "REPEAT" in site.FILTER:
        repeat_sites += 1
        continue

    if site.is_monomorphic:
        valid_sites += 1
        continue

    if site.is_snp:
        if len(site.ALT) > 1:
            multiallelic_snp += 1
            continue
        if site.aaf == 1.0 or site.aaf == 0.0:  # only want SNPs polymorphic in our sample
            valid_sites += 1
            if site.aaf == 1.0:
                snp_not_S += 1 # A SNP different from reference, but not segregating in the samples
            continue

        # if not args.unfiltered:
        if site.FILTER == 'PASS':
            rac = calc_rac(site)
            rac_list.append(rac)
            valid_sites += 1
            continue
        else:
            failed_snp += 1
            continue
        # else:
        #     rac = calc_rac(site)
        #     rac_list.append(rac)
        #     valid_sites += 1
        #     continue
    else:
        problem_site = site.CHROM + '\t' + str(site.POS)
        error_message = 'Could not assign ' + problem_site + ' to site type'
        sys.exit(error_message)

S = len(rac_list)
theta_w = pg.thetaW(n_total, S)
pi = pg.pi(n_total, rac_list)
tajD = round(pg.TajimasD(n_total, S, theta_w, pi), 5)

theta_w_site = round(theta_w / float(valid_sites), 5)
theta_w_percent = theta_w_site * 100

pi_site = round(pi / float(valid_sites), 5)
pi_percent = pi_site * 100


with open(args.outfile, 'w') as outfile:
    print(args.vcf_infile, file=outfile)
    print('S:', S, file=outfile)
    print('thetaW per site (%):', theta_w_site, '({0})'.format(str(theta_w_percent)), file=outfile)
    print('pi per site (%):', pi_site, '({0})'.format(str(pi_percent)), file=outfile)
    print("Tajima's D:", tajD, file=outfile)
    print('valid sites:', valid_sites, file=outfile)
    print('ref N:', ref_N, file=outfile)
    print('indels:', indels, file=outfile)
    print('snp spanning deletion:', spanning_deletion, file=outfile)
    print('low call rate:', low_call_rate, file=outfile)
    print('Sites with extreme depth', extreme_depth, file=outfile)
    print('Sites in repeat regions', repeat_sites, file=outfile)
    print('multiallelic snps:', multiallelic_snp, file=outfile)
    print('failed biallelic snps:', failed_snp, file=outfile)
    print('total sites considered:', valid_sites + ref_N + indels + spanning_deletion + low_call_rate + extreme_depth +
          repeat_sites + multiallelic_snp + failed_snp, file=outfile)
