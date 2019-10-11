#!/usr/bin/env python2
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('T', type=str)
parser.add_argument('P', type=str)
parser.parse_args()

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
def alignSW(s1, s2):
    g = -2
    a = len(s1) + 1
    b = len(s2) + 1  # need one extra cell in row/column
    # matrix with b rows and a columns ; a on first row and b in first column
    SWmatrix = np.zeros((b, a))  # fill in first row and column with zeros instead of gap penalties
    SWdirection = np.zeros((b, a), dtype=str)

    # find diag, vert, and hor scores
    for i in range(1, b):  # rows
        for j in range(1, a):  # columns
            d = SWmatrix[i-1][j-1] + matching(s1[j-1], s2[i-1])  # diagonal, subtract 1 from seq index since +1 for blank first cell
            v = SWmatrix[i-1][j] + g  # vertical
            h = SWmatrix[i][j-1] + g  # horizontal
            o = 0  # other
            SWmatrix[i][j] = max(d, v, h, o)
            # saves direction in matrix
            if d == max(d, v, h, o):
                SWdirection[i][j] = "d"
            elif v == max(d, v, h, o):
                SWdirection[i][j] = "v"
            elif h == max(d, v, h, o):
                SWdirection[i][j] = "h"
            else:
                SWdirection[i][j] = "o"
    # print SWmatrix
    # print SWdirection

    # find high score of matrix and locations for mult alignments
    high_score = np.amax(SWmatrix)
    high_score_index = np.where(SWmatrix == high_score)
    listOfcoordinates = list(zip(high_score_index[0], high_score_index[1]))

    output1 = []
    output2 = []

    # makes i and j coordinates based on high score, starts off traceback at that coordinate
    for k in range(1, len(listOfcoordinates)+1):
        start_coordinate = listOfcoordinates[k-1]
        s_i = start_coordinate[0]  # start point
        s_j = start_coordinate[1]
        # traceback , can have multiple alignments
        s1align = ""
        s2align = ""
        while SWmatrix[s_i][s_j] != 0:  # end alignment when zero is encountered
            if "d" == SWdirection[s_i][s_j]:  # gives preference to diagonal high scores, bp aligned
                s1align = s1[s_j - 1] + s1align  # columns
                s2align = s2[s_i - 1] + s2align  # rows
                s_i -= 1
                s_j -= 1
            elif "v" == SWdirection[s_i][s_j]:  # gap on top sequence
                s1align = s1[s_j - 1] + s1align
                s2align = "-" + s2align
                s_j -= 1
            elif "h" == SWdirection[s_i][s_j]:  # gap on left sequence
                s1align = "-" + s1align
                s2align = s2[s_i - 1] + s2align
                s_i -= 1
        output1.append(s1align)
        output2.append(s2align)


# calculates alignment scores and pipes matching bps
    pipes = []
    for i in range(len(output1)):
        sequence1 = output1[i]
        sequence2 = output2[i]
        pipe = ""
        for j in range(len(sequence1)):
            char1 = sequence1[j]
            char2 = sequence2[j]
            if char1 == char2:
                pipe += "|"
            else:
                pipe += " "
        pipes.append(pipe)

    print ("score: %s" % high_score)
    print ("%s \n%s \n%s" % (output1, pipes, output2))


if __name__ == '__main__':
    alignSW(seq1, seq2)
