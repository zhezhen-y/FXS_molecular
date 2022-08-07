#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 15:20:52 2022

@author: zhezhen
"""

import random
from Bio import SeqIO
from Bio.Seq import Seq

'''
Load fasta file
'''
# Load fasta file using SeqIO
filename = '/Users/zhezhen/Desktop/220322 muSeq Fragile X/Fragile X experiments/fmr1_primer.fasta'
fullseq1 = SeqIO.read(filename, "fasta")
fullseq = str(fullseq1.seq)

# Load fasta file using file IO
#fasta = open(filename,'r')
#lines = fasta.readlines()[1:]
#for i in range(len(lines)):
#    lines[i] = lines[i].rstrip()
#fullseq2 = "".join(lines)
#fasta.close()

'''
A function to partially mutate C to U in a given sequence
'''
def partial_CU(seq, freq=0.5):
    new_seq = ""
    for letter in seq:
        # use random.random() to flip the coin
        if letter == "C" and random.random() < freq:
            new_seq += "U"
        else:
            new_seq += letter
    return new_seq

seq1 = "GCCGCCGCCGCCGCC"
seq2 = "CGGCGGCGGCGGCGG"
partial_CU(seq1)
partial_CU(seq2)




'''
A function to generate reverse complementary sequence 
'''

base_dic ={"G":"C", 
           "C":"G",
           "A":"T",
           "T":"A",
           "U":"A"}

def complementary(seq):
    new_seq = ""
    for letter in seq:
        new_letter = base_dic[letter]
        new_seq += new_letter
    return new_seq

def rev_complementary(seq):
    new_seq = ""
    for letter in seq:
        new_letter = base_dic[letter]
        new_seq += new_letter
    rev_seq = new_seq[::-1]
    return rev_seq