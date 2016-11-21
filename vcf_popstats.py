from __future__ import print_function, division
import sys
import argparse
from pysam import VariantFile, FastaFile
import pgstats as pg
import gzip
import numpy as np
import time


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


def calc_stats_chr(pop, af_lst, chrom, callable_file, out):

    callable_sites = count_callable(callable_file, chrom)

    segs = len(af_lst)
    theta_w = pg.thetaW(n, segs)
    theta_w_site = theta_w / callable_sites[1]

    pi = pg.pi_tajima(n, af_lst)
    pi_site = pi / callable_sites[1]

    tajd = pg.TajimasD(n, segs, theta_w, pi)

    if tajd is None:
        tajd = 'NA'

    print(pop, chrom, callable_sites[0], callable_sites[1], segs, theta_w_site, pi_site, tajd, sep='\t', file=out)

    return callable_sites[0], callable_sites[1]


def calc_stats_region(pop, chrom, feature_id, feature_type, af_lst, callable_sites, out):

    segs = len(af_lst)
    theta_w = pg.thetaW(n, segs)
    theta_w_site = theta_w / callable_sites

    pi = pg.pi_tajima(n, af_lst)
    pi_site = pi / callable_sites

    tajd = pg.TajimasD(n, segs, theta_w, pi)

    if tajd is None:
        tajd = 'NA'

    print(pop, chrom, feature_id, feature_type, callable_sites, segs, theta_w_site, pi_site, tajd, sep='\t', file=out)


def calc_stats_total(af_lst, callable_sites):

    segs = len(af_lst)
    theta_w = pg.thetaW(n, segs)
    theta_w_site = theta_w / callable_sites

    pi = pg.pi_tajima(n, af_lst)
    pi_site = pi / callable_sites

    tajd = pg.TajimasD(n, segs, theta_w, pi)

    if tajd is None:
        tajd = 'NA'

    return callable_sites, segs, theta_w_site, pi_site, tajd


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
parser.add_argument('--pop_id', required=True, dest='pop_id', help="Identifier for population/species to be written in"
                                                                   "the results file")
parser.add_argument('--min', required=False, dest='min_sites', help="Minimum number of sites for region to have stat"
                                                                    "calculated")
parser.add_argument('--weak_strong', required=False, dest='wwss', action='store_true',
                    help="This flag will mean only weak-to-weak (A/T) or strong-to-strong (G/C) sites are used")
parser.add_argument('--sfs', required=False, dest='sfs', help="Compute the folded sfs")
parser.add_argument('--summed_sfs', required=False, dest='sum_sfs', help="Compute the polystats summed accross regions "
                                                                         "in bed")
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

if args.sfs:
    sfs = np.zeros((n // 2) + 1, dtype=int)

if args.exclude:

    af = []
    ac = []
    with open(args.exclude, 'r') as exclude_file:
        exclude_chr = set([i.rstrip() for i in exclude_file])

        contigs = set(vcf_infile.header.contigs)
        contigs = contigs.difference(set(exclude_chr))

    chrom_list = []
    total_af = []
    total_ac = []
    total_sites = 0
    total_callable = 0

    with open(args.outfile, 'w') as outfile:
        print('Population', 'Chromosome', 'Chromosome_length', 'Sites', 'S', 'thetaW', 'pi', 'tajd', sep='\t',
              file=outfile)

        for c in contigs:
            for site in vcf_infile.fetch(c):
                ac = site.info['AC'][0]
                af.append(ac / float(n))

            total_ac += ac
            total_af += af

            sites = calc_stats_chr(args.pop_id, af, c, args.callable, outfile)
            total_sites += sites[0]
            total_callable += sites[1]

            del af[:]

        if args.sfs:
            for c in total_ac:
                sfs[min(n - c, c)] += 1

            sfs[0] += (total_callable - len(total_ac))

        S = len(total_af)
        thetaw = pg.thetaW(n, S)
        thetaw_site = thetaw / total_callable
        pi_taj = pg.pi_tajima(n, total_af)
        pisite = pi_taj / total_callable
        taj_d = pg.TajimasD(n, S, thetaw, pi_taj)

        print(args.pop_id, 'Total', total_sites, total_callable, S, thetaw_site, pisite, taj_d, sep='\t', file=outfile)

if args.bed:
    compressed = False
    bed_dict = {}

    wwss_alleles = [('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')]

    ac_dict = {}
    af_dict = {}
    sites_dict = {}
    ac_list = []
    af_list = []
    sites_list = []

    if args.bed[-3:] == '.gz':
        bed_file = gzip.open(args.bed, 'r')
    elif args.bed[-4:] == '.bed':
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
        print('Population', 'Chromosome', 'Feature_name', 'Feature_type', 'Sites', 'S', 'thetaW', 'pi', 'tajd', sep='\t',
              file=outfile)

        for g in bed_dict:
            sites_dict[g] = 0
            ac_dict[g] = []
            af_dict[g] = []
            for r in bed_dict[g]:
                sites_dict[g] += callable_f.fetch(r[0], r[1], r[2]).count('0')
                vcf_chunk = vcf_infile.fetch(r[0], r[1], r[2])
                for site in vcf_chunk:
                    if args.wwss:
                        if site.alleles in wwss_alleles:
                            ac = site.info['AC'][0]
                            ac_dict[g].append(site.info['AC'][0])
                            af_dict[g].append(ac / float(n))

                    else:
                        ac = site.info['AC'][0]
                        ac_dict[g].append(ac)
                        af_dict[g].append(ac / float(n))

            af_list.append(af_dict[g])
            sites_list.append(sites_dict[g])

            if args.sfs:
                for c in ac_dict[g]:
                    sfs[min(n - c, c)] += 1

                sfs[0] += (sites_dict[g] - len(ac_dict[g]))

            if sites_dict[g] == 0 or sites_dict[g] < int(args.min_sites):
                continue
            else:
                calc_stats_region(args.pop_id, r[0], g, args.type, af_dict[g], sites_dict[g], outfile)

    bed_file.close()

if args.sfs:
    with open(args.sfs, 'w') as sfs_out:
        sfs_str = ' '.join([str(i) for i in sfs])
        print(sfs_str, file=sfs_out)

if args.sum_sfs:

    with open(args.sum_sfs, 'w') as total_outfile:
        print('pop_id', 'feature', 'sites', 'gene_count', 'sites', 'S', 'pi', 'lower_pi', 'upper_pi',
              'thetaW', 'lower_thetaw', 'upper_thetaw', 'tajd', 'lower_tajd', 'upper_tajd', sep='\t',
              file=total_outfile)

        total_af = list(np.concatenate(af_list, axis=0))
        total_sites = sum(sites_list)
        total_polystats = calc_stats_total(total_af, total_sites)

        gene_count = len(af_list)

        total_sites = total_polystats[0]
        total_S = total_polystats[1]
        total_thetaW = total_polystats[2]
        total_pi = total_polystats[3]
        total_tajd = total_polystats[4]

        reps = 10000  # bootstrap reps to perform

        bootstrap_reps = {'Sites': [], 'S_reps': [], 'pi_reps': [], 'thetaw_reps': [], 'tajd_reps': []}

        idx = np.random.randint(0, len(af_list), (reps, len(af_list)))

        af_reps = np.array(af_list)[idx]
        sites_reps = np.array(sites_list)[idx]

        rep_count = 0
        for r, s in zip(af_reps, sites_reps):

            rep_count += 1
            print(time.localtime(), 'Calculating bootstrap rep:', rep_count)

            rep_total_ac = np.concatenate(r, axis=0)
            rep_total_sites = sum(s)

            boot_rep_stats = calc_stats_total(rep_total_ac, rep_total_sites)

            bootstrap_reps['Sites'].append(boot_rep_stats[0])
            bootstrap_reps['S_reps'].append(boot_rep_stats[1])
            bootstrap_reps['thetaw_reps'].append(boot_rep_stats[2])
            bootstrap_reps['pi_reps'].append(boot_rep_stats[3])
            bootstrap_reps['tajd_reps'].append(boot_rep_stats[4])

        pi_upper = np.percentile(bootstrap_reps['pi_reps'], 97.5, axis=0)
        pi_lower = np.percentile(bootstrap_reps['pi_reps'], 2.5, axis=0)
        thetaw_upper = np.percentile(bootstrap_reps['thetaw_reps'], 97.5, axis=0)
        thetaw_lower = np.percentile(bootstrap_reps['thetaw_reps'], 2.5, axis=0)
        tajd_upper = np.percentile(bootstrap_reps['tajd_reps'], 97.5, axis=0)
        tajd_lower = np.percentile(bootstrap_reps['tajd_reps'], 2.5, axis=0)

        if args.wwss:
            site_type = 'wwss'
        else:
            site_type = 'all'

        print(args.pop_id, args.type, site_type, gene_count, total_sites, total_S, total_pi, pi_lower, pi_upper,
              total_thetaW, thetaw_lower, thetaw_upper, total_tajd, tajd_lower, tajd_upper, sep='\t',
              file=total_outfile)
