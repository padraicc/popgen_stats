from pgstats import calc_rac
import sys
import numpy as np


def vcf_site_iterator(vcf, min_dp, max_dp, filtered):

    """
    Iterator that goes site by site through vcf and return
    a list with the reference allele count at SNP sites

    Returns:
        object: list
    """
    rac_list = []

    # indels = 0
    # extreme_depth = 0
    # low_call_rate = 0
    # repeat_sites = 0
    # spanning_deletion = 0
    # snp_not_S = 0
    # multiallelic_snp = 0
    valid_sites = 0
    # failed_snp = 0
    # ref_N = 0

    for site in vcf:

        if site.REF == 'N':
            # ref_N += 1
            continue

        if site.is_indel and site.aaf != 0.0:
            # indels += 1
            continue

        if site.is_indel and site.aaf == 0.0:
            valid_sites += 1
            continue

        if len(site.ALT) >= 1 and site.ALT[-1] == '*':  # SNPs at spanning deletion
            # spanning_deletion += 1
            continue

        all_dp = [x['DP'] for x in site.samples]  # get the genotype depths

        if None in all_dp:  # Exclude sites where a genotype DP field is set to '.'
            # low_call_rate += 1
            continue

        mean_dp = np.mean(all_dp)

        if mean_dp < min_dp or mean_dp > max_dp:  # depth filter
            # extreme_depth += 1
            continue

        if 0 in all_dp or site.call_rate < 1.0:  # only consider sites where all samples have coverage
            # low_call_rate += 1
            continue

        if site.FILTER is not None:
            if site.FILTER == "REPEAT" or "REPEAT" in site.FILTER:  # exclude sites in repeat regions
                # repeat_sites += 1
                continue

        if site.is_monomorphic:
            valid_sites += 1
            continue

        if site.is_snp:

            if len(site.ALT) > 1:
                # multiallelic_snp += 1
                continue

            if site.aaf == 1.0 or site.aaf == 0.0:  # only want SNPs polymorphic in our sample
                valid_sites += 1
                # if site.aaf == 1.0:
                #     snp_not_S += 1 # A SNP different from reference, but not segregating in the samples
                continue

            if filtered:
                if site.FILTER == 'PASS':
                    rac = calc_rac(site)
                    rac_list.append(rac)
                    valid_sites += 1
                    continue
                else:
                    # failed_snp += 1
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

    return rac_list, valid_sites
