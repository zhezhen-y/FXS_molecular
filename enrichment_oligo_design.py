#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 17:37:28 2022

@author: zhezhen
"""

'''
This script is to facilitate designing oligos for the enrichment experiment.

This script contain two parts.

The first part designs the capture oligo sequence, which contains a 5' G end that allows incorporation of cytosine-N3, a spacer and a universal primer, and a blocked 3' end that captures the genomic fragment. 
The genomic DNA is supposed to be cut by restriction enzymes first, and several candidate restriction enzymes are considered.

The second part designs control templates to help troubleshooting and evaluation of the experimental steps.
'''

import random
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from pandas import DataFrame

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

# complementary("GCCGCCGUCGUCGCU")
# complementary(complementary("GCCGCCGUCGUCGCU"))
# rev_complementary("GCCGCCGUCGUCGCU")

# '''
# Design primer pair 3 using UP1
# '''
# # concatenate UP1 sequence and the fragment cut by BclI
# # to check on Primer3 if the primer pairs work or not
# UP1 = "CCATCACTCTATCACATCTACAA"
# seq_cut_captured = UP1 + fullseq[1876:2930]
# seq_captured = UP1 + fullseq[1876:]


'''
Design the primer and spacer sequence
They should not have any Cs in them so that cytosine-N3 would not be added at wrong sites.
'''
def no_G_seq(length):
    seq = [random.choice(["A", "T", "C"]) for i in range(length)]
    return "".join(seq)

def no_G_seq_weighted(length):
    seq = [random.choice(["C", "A", "T", "C"]) for i in range(length)]
    return "".join(seq)

'''
Generate a random c-poor primer, to use after bisulfite-conversion 
'''
def no_C_seq(length):
    seq = [random.choice(["A", "T", "G"]) for i in range(length)]
    return "".join(seq)

# no_C_seq(25)
# 'GAGGTAAGGTAAATGTTTTTTGTGA'
# 'AGTGTTATAAGTGGATTTGGGATAT'

'''
Functions to generate capture oligo annealing sequences for each restriction enzyme
'''
# capture oligo anneals with the 3' of the cut genomic DNA
# note that the cut position and recognization start position are different for 
# a restriction enzyme 
# and the cut position in the bottom strand are also different from the top strand 
# the actual sequence would require some careful check.

# define a function that returns fmr1 gene fragment
def anneal(position, length):
    start = position - 1
    end = position + length -1
    return fullseq[start:end]

# define a function that returns capture oligo sequences
Gs="GGG"
Spacer="TTTACTCCTATCACTTCTCA"
UP1="CCATCACTCTATCACATCTACAA"
block_3prime="/3Phos/"

def capture(position, length):
    seq = Gs + Spacer + UP1 + anneal(position, length) + block_3prime
    return seq

# return the name of the capture oligos
def capture_oligo_name(enzyme):
    oligo_name= f"end_fmr_{enzyme}"
    return oligo_name
 
'''
Read restriction enzyme info from an excel file 
'''
# use python pandas to read from an xlsx file
excel_dir = "/Users/zhezhen/Desktop/220322 muSeq Fragile X/Fragile X experiments/220622 click chemistry/220701 Restriciton Enzyme.xlsx"
enzyme_info_df = pd.read_excel(excel_dir)
print(enzyme_info_df)
enzyme_info_df.columns
enzyme_info_df.index

'''
Design annealing fragment for each enzyme
A test and check block designed for the capture oligos
'''
# get the annealing fragment for each enzyme
anneal_frag_dic = {}
for ind in enzyme_info_df.index:
    # print the annealing seq for each enzyme
    key = enzyme_info_df.loc[ind,'RE']
    if key == "BamHI":
        continue
    position = enzyme_info_df.loc[ind,'Fragment_Start_Bottom_Strand']
    anneal_seq = anneal(position,30)
    anneal_frag_dic["anneal_"+str(key)+"_30"] = anneal_seq

# check that the annealing fragment for each enzyme is correctly designed    
for ind in enzyme_info_df.index:
    key = enzyme_info_df.loc[ind,'RE']
    if key == "BamHI":
        continue
    position = enzyme_info_df.loc[ind,'Fragment_Start_Bottom_Strand']
    anneal_seq = anneal(position,30)
    anneal_key = "anneal_"+str(key)+"_30"
    print(anneal_key,anneal_seq)
    # check that the cut sites match the genomic positions in the reference genome
    rec_seq = enzyme_info_df.loc[ind,'Cut_Site']
    rec_pos = int(enzyme_info_df.loc[ind,'Sequence_Start_Position'])-1
    query = fullseq[rec_pos:rec_pos + 6]
    # print out the anneal sequences, added with the cut bases to double check 
    base_cut = enzyme_info_df.loc[ind,'Bases_Cut_Bottom_Strand']
    rec_pos = enzyme_info_df.loc[ind,'Sequence_Start_Position']
    print(key,rec_seq,query, anneal(rec_pos,30+base_cut),"\n")

# output each part of the oligos into an excel file

'''
Output the full capture oligo sequences
'''
# Store the full capture sequences in a dictionary
capture_oligo_dic = {}
for ind in enzyme_info_df.index:
    key = enzyme_info_df.loc[ind,'RE']
    if key == "BamHI":
        continue
    position = enzyme_info_df.loc[ind,'Fragment_Start_Bottom_Strand']
    capture_oligo_dic[capture_oligo_name(key)+"_30"] = capture(position,30)
    capture_oligo_dic[capture_oligo_name(key)+"_45"] = capture(position,45)

# Convert the dictionary into dataframe for exporting
anneal_frag_df = DataFrame({"name":anneal_frag_dic.keys(),"seq":anneal_frag_dic.values()})
capture_oligo_df = DataFrame({"oligo":capture_oligo_dic.keys(),"seq":capture_oligo_dic.values()})

# Export the oligo sequences into an excel file for ordering on IDT website
with pd.ExcelWriter('220705_oligo_design.xlsx') as writer:
    capture_oligo_df.to_excel(writer, sheet_name='capture oligo full seq', index=False)
    anneal_frag_df.to_excel(writer, sheet_name='anneal fragments', index=False)
   
'''
Design control templates
'''
# Design an artificial cpoor primer
# In the future, can change this into genomic c_poor regions
cpoor_primer='GAGTAAGTAAATGTTGTGAAGTGTTATAAGTGATTTGGATAT'[0:25]

# A function for generating the varietal tag sequence
# for quantifing the efficiency of each experimental step
def make_vt_seq(length):
    return "N" *length  
make_vt_seq(20) 

# Add some sequence that contains Cs so I can test whether bisulfite works
# fullseq[1420:1440]
# 'CACCCTGGTTGCCACTGTTC'
short_repeat='CCG'*5

# A function for generating control templates
def control_seq(anneal_seq):
    seq = "/5Phos/" + cpoor_primer + make_vt_seq(20) + short_repeat + anneal_seq
    return seq

# Generate control templates for each enzyme
anneal_rc_dic = {}
control_seq_dic = {}
control_primer_dic = {}
for ind in enzyme_info_df.index:
    # loop through the enzymes
    key = enzyme_info_df.loc[ind,'RE']
    if key == "BamHI":
        continue
    position = enzyme_info_df.loc[ind,'Fragment_Start_Bottom_Strand']
    # anneal_rc_seq is the 3'end of the control template, sequence from 5' to 3'
    # 30bp should be enough for annealing two oligos
    anneal_rc_seq = rev_complementary(anneal(position,30))
    # Generate primers to amplify single-stranded control templates into double-strands
    # Control templates need to be amplified into double-strand templates for a good quantification of the input amount. They come as single-stranded when ordered from IDT.
    control_primer = "/5Phos/" + rev_complementary(anneal_rc_seq[-18:])
    # Generate control template full sequences
    seq = control_seq(anneal_rc_seq)
    # Store the information in dictionaries
    anneal_rc_dic["anneal_rc_"+str(key)+"_30"] = anneal_rc_seq
    control_primer_dic["ctr_primer_"+str(key)+"_30"] = control_primer
    control_seq_dic["control_"+str(key)+"_30"] = seq


    
# Convert the dictionary into dataframe for exporting
control_template_df = DataFrame({"Name":control_seq_dic.keys(),"Sequence":control_seq_dic.values()})
control_primer_df = DataFrame({"Name":control_primer_dic.keys(),"Sequence":control_primer_dic.values()})

# Export the oligo sequences into an excel file for ordering on IDT website
with pd.ExcelWriter('220705_control_template_design.xlsx') as writer:
    control_template_df.to_excel(writer, sheet_name='control template full seq', index=False)
    control_primer_df.to_excel(writer, sheet_name='control primer seq', index=False)



