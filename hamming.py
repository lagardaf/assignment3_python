#!/usr/bin/env python2
#	Two sequences are: 1.	TTATCGCGCTTTCTTCCAA  2. ATACCGCGCGTTCGACCAA
# minimum number of changes needed to change seq into the other


def hamming(s1, s2):
    distance = 0
    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            distance = distance + 1
    print distance
    return distance


seq1 = "TTATCGCGCTTTCTTCCAA"
seq2 = "ATACCGCGCGTTCGACCAA"
hamming(seq1, seq2)
