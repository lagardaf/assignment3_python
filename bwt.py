#!/usr/bin/env python2

# burrows-wheeler matrix
T = "AATTGCGCGG"  # 10 bps


def bwt(seq):
    seq += "$"
    l = len(seq)
    rotations = [0] * l

    #  moves last letter of sequence to start
    for i in range(1, l+1):
        rotations[i-1] = seq[l-i:l] + seq[0:l-i]
    rotations_alpha = sorted(rotations)    # lists in alphabetical order
    k = len(rotations_alpha)
    bwt = ""
    temp = ""

    #   takes last letter of each string
    for i in range(k):
        temp += rotations_alpha[i]
        bwt += temp[l-1]
        temp = ""

    #print bwt


bwt(T)
