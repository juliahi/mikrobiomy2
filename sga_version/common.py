# widely used methods methods
import math

epsilon = 0.000001


def compl(s):
    if s == '' or s == []: return ''
    if s is None: return None

    def c(x):
        if x == 'A':
            return 'T'
        elif x == 'C':
            return 'G'
        elif x == 'G':
            return 'C'
        elif x == 'T':
            return 'A'
        return 'N'

    return ''.join(map(c, s))


def cons_pairs(input_list):
    for i in xrange(len(input_list)-1):
        yield (input_list[i], input_list[i+1])


# def foldchange(n1, n2):
#         if n2 == 0: return float("inf")
#         return 1. * n1 / n2

def foldchange(n1, n2):
    if n2 == 0: n2 = epsilon
    if n1 == 0: n1 = epsilon
    return 1. * n1 / n2


# def abslog2foldchange(n1, n2):
#         return abs(math.log(foldchange(n1, n2), 2)) if n1 != 0 else float("-inf")
def abslog2foldchange(n1, n2):
        return abs(math.log(foldchange(n1, n2), 2))


def foldchange_compare(n1, n2, min_fc):
        if n1 >= n2*min_fc: return True
        if n2 >= n1*min_fc: return True
        return False


# check if foldchange for counts c1 stays enriched in the same direction aftr adding counts c2
def foldchange_dir(c1, c2, min_fc):
        if c1[0] >= c1[1]*min_fc and c1[0]+c2[0] >= (c1[1]+c2[1])*min_fc:
            return True
        if c1[1] >= c1[0]*min_fc and c1[1]+c2[1] >= (c1[0]+c2[0])*min_fc:
            return True
        return False