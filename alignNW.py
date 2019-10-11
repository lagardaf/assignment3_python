#!/usr/bin/env python2
import numpy as np

seq1 = "TATCGCGCTTT"  # 11bp
seq2 = "ATTACCGCCGTT"  # 12 bp


# function checks  if mismatch or match ; Match = +1; Mismatch = -1; Gap = -2
def matching(s1, s2):
    g = -2
    a, b = len(s1) + 1, len(s2) + 1  # need one extra cell in row/column
    for i in range(max(a, b)):
        if s1[i] == s2[i]:
            s = 1
            return s
        elif s1[i] != s2[i]:
            s = -1
            return s
        else:
            return g


# function fills in matrix with highest scores and direction
def alignNW(s1, s2):
    g = -2
    a = len(s1) + 1
    b = len(s2) + 1  # need one extra cell in row/column
    NWmatrix = np.zeros((b, a))  # matrix with b rows and a columns since we want a on first row and b in first column
    NWdirection = np.zeros((b, a), dtype=str)

    # fill in first row and column with gap penalties
    for i in range(a):
        NWmatrix[0][i] = g * i
    for j in range(b):
        NWmatrix[j][0] = g * j

    # find diag, vert, and hor scores
    for i in range(1, b):  # rows
        for j in range(1, a):  # columns
            d = NWmatrix[i-1][j-1] + matching(s1[j-1], s2[i-1])  # diagonal, subtract 1 from sequence index since +1 for blank first cell
            v = NWmatrix[i-1][j] + g  # vertical
            h = NWmatrix[i][j-1] + g  # horizontal
            NWmatrix[i][j] = max(d, v, h)
            # saves direction in matrix
            if d == max(d, v, h):
                NWdirection[i][j] = "d"
            elif v == max(d, v, h):
                NWdirection[i][j] = "v"
            else:
                NWdirection[i][j] = "h"
    score = NWmatrix[b-1][a-1]
# traceback, only one alignment
# diagonal is aligned, horizontal is a gap of left sequence, vertical is gap on top sequence
    i = len(s2)
    j = len(s1)
    s1align = ""
    s2align = ""
    while j > 0 and i > 0:
        if "d" == NWdirection[i][j]:  # gives preference to diagonal high scores
            s1align += s1[j-1]  # columns
            s2align += s2[i-1]  # rows
            i -= 1
            j -= 1
        elif "v" == NWdirection[i][j]:
            s1align += s1[j-1]
            s2align += "-"
            j -= 1
        elif "h" == NWdirection[i][j]:
            s1align += "-"
            s2align += s2[i-1]
            i -= 1

    pipes = ""
    s1align = s1align[::-1]
    s2align = s2align[::-1]
# pipes matching bps
    for i in range(len(s1align)):
        if s2align[i] == s1align[i]:
            pipes += "|"
        elif s2align[i] == "-" or s1align[i] == "-":
            pipes += " "
        elif s2align[i] != s1align[i]:
            pipes += " "

    print ("score: %s" % score)
    print ("%s \n%s \n%s" % (s1align, pipes, s2align))


if __name__ == '__main__':
    alignNW(seq1, seq2)
