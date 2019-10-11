#!/usr/bin/env python2

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('T', type=str)
parser.add_argument('P', type=str)
parser.parse_args()

# burrows-wheeler matrix
T = "AATTGCGCGG"  # 10 bps


def bwIndex(seq):
    seq += "$"
    l = len(seq)
    bwm = [0] * l
    bwm_index = [0] * l

    #  moves last letter of sequence to start
    for i in range(1, l+1):
        bwm[i-1] = seq[l-i:l] + seq[0:l-i]
        bwm_index[i-1] = l-i

    bwm_alpha = sorted(bwm)    # lists in alphabetical order
    alpha_in_bwm = [0] * l
    bwm_alpha_index = [0] * l

    for i in range(1, l+1):
        alpha_in_bwm[i-1] = bwm.index(bwm_alpha[i-1])
        bwm_alpha_index[i-1] = bwm_index.index(alpha_in_bwm[i-1])

    k = len(bwm_alpha)
    first = ""
    last = ""
    temp = ""
    #   takes first/ last letter of each string
    for i in range(k):
        temp += bwm_alpha[i]
        first += temp[0]
        last += temp[l-1]
        temp = ""

    print ("first %s last \n%s index \n%s" % (first, last, bwm_alpha_index))


bwIndex(T)
