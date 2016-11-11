import math


def thetaW(n, s):
    """ Calculate's Waterson's 1975 estimator of theta 
    Takes n the sample number and rac, a list of reference allele 
    counts at each segregating site. For a fasta file one of the alleles
    at a biallelic site is chosen as the ref

    n: sample size
    s: number of segregating sites"""

    if s == 0:
        tw = 0.0
    else:        
        a1 = sum(1.0 / i for i in range(1, n))
        tw = s / a1
    
    return tw


def pi_tajima(n, ac):

    """Equation 11 and 12 from Tajima 1989
    Only consider biallelic sites

    n: sample size
    ac: list of alternate allele count """

    s = len(ac)
    if s == 0:
        pi = 0.0
    else:
        pi = 0.0
        for i in ac:
            p_squared = (i / float(n))**2
            q_squared = (1 - (i / float(n)))**2
            pi += (n * (1 - sum([p_squared, q_squared]))) / float(n - 1)
                                 
    return pi


def pi_sfs(sfs):

    """Calculate pi using the folded sfs

    sfs is a list containing the minor allele counts
    """

    s = sum(sfs)
    n = len(sfs) * 2
    if s == 0:
        pi_total = 0.0
    else:
        pi = 0.0
        i = 1
        for Si in sfs:
            pi += Si * i * (n - i)
            i += 1

    pi_total = (2 / float(n * (n - 1))) * pi

    return pi_total


def TajimasD(n, s, tw, pi):
    
    if s == 0:
        D = None
    else:
        rawD = pi - tw
        
        a1 = 0.0
        a2 = 0.0
        for i in range(1, n):
            a1 += 1. / i
            a2 += 1. / (i * i)
        b1 = (n + 1.) / (3. * (n - 1.))
        b2 = 2. * (n * n + n + 3.) / (9. * n * (n - 1.))
        c1 = b1 - 1. / a1        
        c2 = b2 - (n + 2.) / (a1 * n) + a2 / (a1 * a1)
        e1 = c1 / a1
        e2 = c2 / (a1 * a1 + a2)
        V = e1 * s + e2 * s * (s - 1.)

        D = rawD / math.sqrt(V)
    
    return D


def calc_delta_pi(S, pi, n):

    """
    Calculation of delta_pi from Langley et al. (2014) PLoS Genet 10(7): e1004457.

    :param S: Segregating sites
    :param pi: nucleotide diversity
    :param n: Sample size
    :return: float

    """

    if S > 0:
        assert isinstance(S, int)
        delta_pi = (pi / float(S)) - (1 / sum(1.0 / i for i in range(1, n)))
    else:
        delta_pi = 'NA'

    return delta_pi


def sfs(ac, n):
    """Returns the folded site frequency spectrum as a list"""
    pass


def pib(ac_1, ac_2, n1, n2):
    """Calculates the pi between (pib) populations 1 and 2 (aka Dxy) in a window
    by calculating piB at each site where the alternate allele frequency in pop1 is p1
    and pop2 p2, using the equation p1*(1 - p2) + p2*(1 - p1)."""
    
    if sum(ac_1) + sum(ac_2) == 0:
        pi_b = 0.0
    else:        
        pi_b = 0.0
        for s in xrange(len(rac_1)):
            p1 = ac_1[s] / float(n1)
            p2 = ac_2[s] / float(n2)
            pi_b += p1 * (1 - p2) + p2 * (1 - p1)
                
    return pi_b
            

def fst(pi_b, ac_1, ac_2, n1, n2):
    """ Calculates the Hudson et al. (1992) unweighted version of Fst
        Fst = (piB - piS) / piS"""
    
    ##TODO add weithted Fst

    if sum(ac_1) + sum(ac_2) == 0:
        Fst = 0.0
    else:
        Fst = 0.0
        for s in xrange(len(ac_1)):
            p1 = ac_1[s] / n1
            p2 = ac_1[s] / n2
        
            piS = p1*(1-p1) + p2*(1-p2)
    
            Fst = (pi_b - piS) / piS

    return Fst

    

##TODO add some LD stats (e.g., r^2, D', Kelly's ZnS)


