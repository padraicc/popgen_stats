from __future__ import print_function, division
import sys
import argparse
from pysam import VariantFile, FastaFile
import pgstats as pg
import gzip


def get_feature_name(feature_rec):

    """
    retrieves the gene name from the fifth column of the bed fiel

    :param feature_rec: string
    :return: feat_name
    """

    atts = feature_rec.split(';')

    for f in atts:
        if f.startswith('ID='):
            feat_name = f.split('=')[1]

    return feat_name


def count_callable(callable_file, chrom, start=None, end=None):

    callable_file = FastaFile(callable_file)

    seq = callable_file.fetch(chrom, start, end)
    chrom_length = len(seq)
    callable_sites = seq.count('0')

    return chrom_length, callable_sites


def calc_stats_chr(ac_list, chrom, callable_file, out):

    callable_sites = count_callable(callable_file, chrom)

    segs = len(ac_list)
    theta_w = pg.thetaW(n, segs)
    theta_w_site = round(theta_w / callable_sites[1], 5)

    pi = pg.pi_tajima(n, ac_list)
    pi_site = round(pi / callable_sites[1], 5)

    tajd = pg.TajimasD(n, segs, theta_w, pi)

    if tajd is None:
        tajd = 'NA'

    print(chrom, callable_sites[0], callable_sites[1], segs, theta_w_site, pi_site, tajd, sep='\t', file=out)

    return callable_sites[0], callable_sites[1]


def calc_stats_region(chrom, feature_id, feature_type, ac_list, callable_sites, out):

    segs = len(ac_list)
    theta_w = pg.thetaW(n, segs)
    theta_w_site = round(theta_w / callable_sites, 5)

    pi = pg.pi_tajima(n, ac_list)
    pi_site = round(pi / callable_sites, 5)

    tajd = pg.TajimasD(n, segs, theta_w, pi)

    if tajd is None:
        tajd = 'NA'

    print(chrom, feature_id, feature_type, callable_sites, segs, theta_w_site, pi_site, tajd, sep='\t', file=out)


parser = argparse.ArgumentParser(description="Program to calculate population genetic statistics from a region "
                                             "in VCF file")
parser.add_argument('-i', '--vcf', required=True, dest='vcf_infile',
                    help="VCF input file containing biallelic SNPs (Needs to be filtered)")
parser.add_argument('-o', '--out', required=True, dest='outfile', help="Outfile to write results")
parser.add_argument('-p', '--ploidy', required=True, type=int, dest='ploidy',
                    help="Ploidy of samples. 1 for haploid and 2 for diploid")
parser.add_argument('-f', '--callable', required=True, dest='callable', help="Callable sites in a fasta format. "
                                                                             "Callable is 0 and not-callable is 1")
parser.add_argument('-e', '--exclude', required=False, dest='exclude', help="File listing chromosome/Scaffolds to "
                                                                            "exclude (Can not be use when -b specified")
parser.add_argument('-b', '--bed', required=False, dest='bed', help="Bed file specifying regions to calculate stats"
                                                                    "(e.g., bedfile specifying introns). Columns are "
                                                                    "Chromosome\tstart\tend\tstrand\tfeature_name\t"
                                                                    "feature_type. Can not be "
                                                                    "used when -e is specified")
parser.add_argument('-t', '--type', required=True, dest='type', help="Name for the feauture given in the bed file "
                                                                     "specified by b. Needed for printing the output."
                                                                     "(e.g., -t intron in the case that the bed file "
                                                                     "gives the coordinates for introns")
parser.add_argument('--min', required=False, dest='min_sites', help="Minimum number of sites for region to have stat"
                                                                    "calculated")
parser.add_argument('--weak_strong', required=False, dest='wwss', action='store_true',
                    help="This flag will mean only weak-to-weak (A/T) or strong-to-strong (G/C) sites are used")
# parser.add_argument('-x', '--pop_1', required=False, dest='pop1', help="Samples belonging to population 1")
# parser.add_argument('-y', '--pop_2', required=False, dest='pop2', help="Samples belonging to population 2")


args = parser.parse_args()

if args.bed and args.exclude:
    sys.exit('\nCan not use the -e and -b options together\n')

if args.wwss and not args.bed:
    sys.exit('\nNeed to specify a bed file of regions with -b option\n')

if args.bed and not args.min_sites:
    sys.exit('\nThe --min option must be specified with -b option\n')


vcf_infile = VariantFile(args.vcf_infile)
sample_num = len(vcf_infile.header.samples)

if args.ploidy == 2:
    n = 2 * sample_num
elif args.ploidy == 1:
    n = sample_num
else:
    sys.exit("\nSpecify ploidy as 1 or 2\n")

ac = []

if args.exclude:
    with open(args.exclude, 'r') as exclude_file:
        exclude_chr = set([i.rstrip() for i in exclude_file])

        contigs = set(vcf_infile.header.contigs)

        contigs = contigs.difference(set(exclude_chr))

    chrom_list = []
    total_ac = []
    total_sites = 0
    total_callable = 0

    with open(args.outfile, 'w') as outfile:
        print('Chromosome', 'Chromosome_length', 'Sites', 'S', 'thetaW', 'pi', 'tajd', sep='\t', file=outfile)

        for c in contigs:
            for site in vcf_infile.fetch(c):
                ac.append(site.info['AC'][0])

            total_ac += ac

            sites = calc_stats_chr(ac, c, args.callable, outfile)
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

if args.bed:
    compressed = False
    bed_dict = {}

    if args.wwss:
        wwss_alleles = [('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')]
        non_wwss_poly = 0

    if args.bed[-3:] == '.gz':
        bed_file = gzip.open(args.bed, 'r')
    elif args.bed[-3:] == '.gz':
        bed_file = open(args.bed, 'r')
    else:
        sys.exit("\nIs this a bed file? Is it compressed?\n")

    callable_f = FastaFile(args.callable)

    for rec in bed_file:
        col = rec.split()
        if ';' in col[4]:
            feature_name = get_feature_name(col[4])
        else:
            feature_name = col[4]

        if feature_name not in bed_dict:
            bed_dict[feature_name] = []

        bed_dict[feature_name].append((col[0], int(col[1]), int(col[2])))

    with open(args.outfile, 'w') as outfile:
        print('Chromosome', 'Feature_name', 'Feature_type', 'Sites', 'S', 'thetaW', 'pi', 'tajd', sep='\t',
              file=outfile)

        for g in bed_dict:
            sites = 0
            for r in bed_dict[g]:
                sites += callable_f.fetch(r[0], r[1], r[2]).count('0')
                vcf_chunk = vcf_infile.fetch(r[0], r[1], r[2])
                for site in vcf_chunk:
                    if args.wwss:
                        if site.alleles in wwss_alleles:
                            # print(site.alleles)
                            ac.append(site.info['AC'][0])
                        else:
                            non_wwss_poly += 1
                    else:
                        ac.append(site.info['AC'][0])

            # print(sites, non_wwss_poly)
            if args.wwss:
                sites = sites - non_wwss_poly
                non_wwss_poly = 0

            # print(sites, non_wwss_poly)

            if sites == 0 or site < int(args.min_sites):
                # print(g, sites, len(ac))
                del ac[:]
                continue
            else:
                calc_stats_region(r[0], g, args.type, ac, sites, outfile)
                del ac[:]

    bed_file.close()
