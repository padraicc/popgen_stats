from __future__ import print_function, division
import sys
import argparse
from pysam import VariantFile, FastaFile
import pgstats as pg


def calc_stats(ac_list, chrom, callable_file, out):

    callable_f = FastaFile(callable_file)

    seq = callable_f.fetch(reference=chrom)
    chrom_length = len(seq)
    callable_sites = seq.count('0')

    S = len(ac_list)
    theta_w = pg.thetaW(n, S)
    theta_w_site = round(theta_w / callable_sites, 5)

    pi = pg.pi_tajima(n, ac_list)
    pi_site = round(pi / callable_sites, 5)

    tajd = round(pg.TajimasD(n, S, theta_w, pi), 5)

    print(chrom, chrom_length, callable_sites, S, theta_w_site, pi_site, tajd, sep='\t', file=out)

    return chrom_length, callable_sites


parser = argparse.ArgumentParser(description="Program to calculate population genetic statistics from a region "
                                             "in VCF file")
parser.add_argument('-i', '--vcf', required=True, dest='vcf_infile',
                    help="VCF input file containing biallelic SNPs (Needs to be filtered)")
parser.add_argument('-o', '--out', required=True, dest='outfile', help="Outfile to write results")
parser.add_argument('-p', '--ploidy', required=True, type=int, dest='ploidy', help="Ploidy of sample. 1 for haploid and 2 for "
                                                                          "diploid")
parser.add_argument('-f', '--callable', required=False, dest='callable', help="Callable sites in a fasta fomrmat. "
                                                                            "Callable is 0 and not-callable is 1")
parser.add_argument('-e', '--exclude', required=False, dest='exclude', help="File listing chromosome/Scaffolds to "
                                                                            "exclude")
# parser.add_argument('-x', '--pop_1', required=False, dest='pop1', help="Samples belonging to population 1")
# parser.add_argument('-y', '--pop_2', required=False, dest='pop2', help="Samples belonging to population 2")
# parser.add_argument('-b', '--bed', required=False, dest='region', help="Bed file specifying regions to calculate stats
# (e.g., bedfile specifying introns")

args = parser.parse_args()

vcf_infile = VariantFile(args.vcf_infile)
sample_num = len(vcf_infile.header.samples)

contigs = set(vcf_infile.header.contigs)


if args.ploidy == 2:
    n = 2 * sample_num
elif args.ploidy == 1:
    n = sample_num
else:
    sys.exit("Specify ploidy as 1 or 2")

if args.exclude:
    with open(args.exclude, 'r') as exclude_file:
        exclude_chr = set([i.rstrip() for i in exclude_file])

        contigs = contigs.difference(set(exclude_chr))

ac = []
chrom_list = []
total_ac = []
total_sites = 0
total_callable = 0
site_num = 0


with open(args.outfile, 'w') as outfile:
    print('Chromosome', 'Chromosome_length', 'Sites', 'S', 'thetaW', 'pi', 'tajd', sep='\t', file=outfile)

    for c in contigs:
        for site in vcf_infile.fetch(c):
            ac.append(site.info['AC'][0])

        total_ac += ac

        sites = calc_stats(ac, c, args.callable, outfile)
        total_sites += sites[0]
        total_callable += sites[1]

        del ac[:]

    S = len(total_ac)
    thetaw = pg.thetaW(n, S)
    thetaw_site = round(thetaw / total_callable, 5)
    pi_taj = pg.pi_tajima(n, total_ac)
    pisite = round(pi_taj / total_callable, 5)
    taj_d = round(pg.TajimasD(n, S, thetaw, pi_taj), 5)

    print('Total', total_sites, total_callable, S, thetaw_site, pisite, taj_d, sep='\t', file=outfile)
