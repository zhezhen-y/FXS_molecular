#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 14:52:19 2022

@author: zhezhen
"""
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import random

'''
Load a piece of the fasta sequence of the CGG repeat region
https://www.ncbi.nlm.nih.gov/nuccore/1022943339
>NG_007529.2:3000-10000 Homo sapiens FMRP translational regulator 1 (FMR1), RefSeqGene (LRG_762) on chromosome X

'''

# Load fasta file using SeqIO
filename = '/Users/zhezhen/Desktop/fmr1_primer.fasta'
fmrseq = SeqIO.read(filename, "fasta")
fmrseq = str(fmrseq.seq)

'''
Input primer sequences
'''
forward_asu = "TCAGGCGCTCAGCTCCGTTTCGGTTTCA"
reverse_asu = "AAGCGCCATTGGAGCCCCGCACTTCC"

'''
Find the forward primer start index in the sequence
'''
start = fmrseq.find(forward_asu) + 1 

'''
Reverse complement the reverse primer
'''
reverse_rc = str(Seq(reverse_asu).reverse_complement())

'''
Find the reverse primer start index in the sequence
'''
end = fmrseq.find(reverse_rc) + len(reverse_asu)
'''
Calculate the product length
'''
product_len = end - start + 1

'''
Write a function
'''
def primer_design(forward,reverse):
    start = fmrseq.find(forward) + 1
    reverse_rc = str(Seq(reverse).reverse_complement())
    end = fmrseq.find(reverse_rc) + len(reverse)
    product_len = end - start + 1
    forward_dic = {"start":start,"len":len(forward)}
    reverse_dic = {"end":end,"len":len(reverse)}
    return (forward_dic,reverse_dic,product_len)

def cpoor_primer_design(forward,reverse_rc):
    start = fmrseq.find(forward) + 1
    end = fmrseq.find(reverse_rc) + len(reverse_rc)
    product_len = end - start + 1
    reverse = str(Seq(reverse_rc).reverse_complement())
    forward_dic = {"start":start,"len":len(forward),"seq":forward}
    reverse_dic = {"end":end,"len":len(reverse),"seq":reverse}
    return (forward_dic,reverse_dic,product_len)

'''
Load sequences for c-poor primers
'''
# define a function that removes the whitespaces in a sequence
def rm_spaces(seq):
    return seq.strip().replace(' ','')

cpoor_left = 'CCT TCC ACT CCA CCT CCC '.strip().replace(' ','')
cpoor_right = 'CCCCATCTTCGCCCTTCCTTCCCTCCC'.strip().replace(' ','')
cpoor_primer_design(cpoor_left, cpoor_right)

'''
Convert the fasta sequences into 0s and 1s
'''
target_bases = {'plus strand': 'C', 'minus strand': 'G'}
target_base = 'G'

def binary(base):
    if base == target_base:
        base = 0
    else:
        base = 1
    return base

def binary_list(seq):
    binary_seq = []
    for i in range(len(seq)):
        base = seq[i]
        binary_seq.append(binary(base))
    return binary_seq

binary_fmr = binary_list(fmrseq)

def list_to_string(list):
    string_list = [str(i) for i in list]
    return "".join(string_list)
        
print(list_to_string(binary_fmr))

def score_fragments(binary_seq):
    var = []
    var.append(binary_seq[0])
    for i in range(len(binary_seq)):
        if i == 0:
            continue
        base_score = binary_seq[i]
        if base_score == 0:
            var.append(0)
        else:
            var.append(var[i-1] + 1)
    return var

fmr_seq_score = score_fragments(binary_fmr)
        
        


'''
Pick protection oligos for the primer
'''