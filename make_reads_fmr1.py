#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 13:03:38 2021

@author: zhezhen
"""

# Museq Simulation
# Simulate Fastq Reads
# See if the museq pipeline can assemble the contigs, phase the haplotypes, 
# and identify two linked SNPs from the simulated reads


from Bio import SeqIO
from Bio.Seq import Seq
from subprocess import check_call
import numpy as np
import random
import os
import sys

random.seed(42)


# set the output directory, total samples(initial molecules) and coverage
# in the future: add parameters such as [repeat units] and [C-T conversion rate]
TESTING = True
if TESTING:
    dir             = "/Volumes/data/safe/zhezhen/FMR1/test1_100samples/input"
    total_sample    = 10
    coverage        = 1000

else:
    dir             = sys.argv[1]
    total_sample    = int(sys.argv[3])
    coverage        = int(sys.argv[4])
    transition_prob = int(sys.argv[2])

# make the input directory if not already exist
if not os.path.isdir(dir):
    os.makedirs(dir)

'''
PART 0: Make the templates for expanded FMR1 alleles
'''
# Input the top strand sequence of the CGG repeat region 
seq_upstream = "CCTTCCCGCCCTCCACCAAGCCCGCGCACGCCCGGCCC\
GCGCGTCTGTCTTTCGACCCGGCACCCCGGCCGGTTCCCAGCAGCGCGCATGCGCGCGCTCCCAGGCCAC\
TTGAAGAGAGAGGGCGGGGCCGAGGGGCTGAGCCCGCGGGGGGAGGGAACAGCGTTGATCACGTGACGTG\
GTTTCAGTGTTTACACCCGCAGCGGGCCGGGGGTTCGGCCTCAGTCAGGCGCTCAGCTCCGTTTCGGTTT\
CACTTCCGGTGGAGGGCCGCCTCTGAGCGGGCGGCGGGCCGACGGCGAGCGCGGGCGGCGGCGGTGACGG\
AGGCGCCGCTGCCAGGGGGCGTGCGGCAGCG"

seq_downstream = "CTGGGCCTCGAGCGCCCGCAGCCCACCTCTCGGGGGCGGGCTCCCGGCG\
CTAGCAGGGCTGAAGAGAAGATGGAGGAGCTGGTGGTGGAAGTGCGGGGCTCCAATGGCGCTTTCTACAA\
GGTACTTGGCTCTAGGGCAGGCCCCATCTTCGCCCTTCCTTC"

repeat_seq = "CGG"

# Define a function that gives the fullseq of the primer flanked region with variable repeats
def make_fullseq(seq_upstream,seq_downstream,repeat_seq,repeat_number):
    repeats = repeat_seq * repeat_number
    fullseq = seq_upstream + repeats + seq_downstream
    return fullseq

# Make the top strands of the FMR1 gene
normal = make_fullseq(seq_upstream, seq_downstream, repeat_seq, 30)
premutation = make_fullseq(seq_upstream, seq_downstream, repeat_seq, 100)
fgx = make_fullseq(seq_upstream, seq_downstream, repeat_seq, 300)

# Define a function that returns the reverse complementary strands using Bio.Seq
def rc(seq):
    return str(Seq(seq).reverse_complement())

'''
PART 1: Make fake haplotype 1 and haplotype 2 templates for expanded FMR1 alleles
        the current museq pipeline uses SNPs to differentiate haplotypes
        add the SNPs so that the pipeline can go through normally
'''

# Make the template for haplotype 1
# Use the bottom strand of the allele
# The bottom strand is what we target in the experiment
haplotype1 = rc(normal)

# Make the template for haplotype 2
# Pick the positions of the fake SNPs
snp1_pos = 100
snp2_pos = -100
# Get the nucleotides at those positions
snp1_before = haplotype1[snp1_pos]
snp2_before = haplotype1[snp2_pos]
# Define a function to randomly mutate nucleotides
def mutate_base(old_base):
    new_base = random.choice(["A","C","G"])
    if new_base == old_base:
        new_base = "T"
    return new_base
# Pick the two SNPs
snp1_after = mutate_base(snp1_before)
snp2_after = mutate_base(snp2_before)
# Make the haplotype2
haplotype2 = list(haplotype1)
haplotype2[snp1_pos] = snp1_after
haplotype2[snp2_pos] = snp2_after
haplotype2 = "".join(haplotype2)

# get the template length
template_length = len(haplotype1)

'''
PART 2: Sample from the haplotype 1 and haplotype 2
'''
# Frequency of haplotype 1 
p = 0.5

# Simulate the sampling process of a binomial distribution
# If the result is 1, then a wildtype sample(haplotype1) is drawn; 
# if the result is 0, then a non-wildtype(haplotype2) sample is drawn 
sample_flag = np.random.binomial(1, p, size = total_sample)
# The list of sequences of the 50 samples 
samples_ori = [haplotype1 if i==1 else haplotype2 for i in sample_flag]

'''
PART 3: Add the C-T mutation footprints to the templates
'''
# Assume the probability of C-T transitions
transition_prob = 0.5
# Mutate the templates
samples_mu = []            
for template in samples_ori:
    new_template = ""
    for letter in template:
        # use random.random() to flip the coin
        if letter == "C" and random.random() < transition_prob:
            new_template += "T"
        else:
            new_template += letter
    samples_mu.append(new_template)
    
# Import the mutated templates in a fasta file 
fa = open(os.path.join(dir, "demo.mutated.fasta") , "w")
for ind, template in enumerate(samples_mu):
    seq_label = ['haplotype1' if sample_flag[ind]==1 else 'haplotype2'][0]
    seq_name = "_".join([str(ind),seq_label])
    seq_block =[">" + seq_name, template]
    fa.write("\n".join(seq_block) + "\n")
fa.close()

'''
PART 4: Fragment the templates
'''
# # Simulate the coverage of the sequencing data
# coverage = 50

# Set the read length for the read pairs (150bp*2)
read_length = 150
# Calculated the total number of reads based on the coverage
# coverage = total number of read pairs * read lengths / (number of templates * template length)
total_reads = coverage * total_sample * template_length // (read_length * 2)  

# Set the parameters for the fake sequencing data
# MINL and MAXL are the smallest and longest half-length of the fragments 
parameters = {"total_reads" : total_reads,
              "template_length" : template_length,
              "total_sample" : total_sample,
              "MINL" : 150,
              "MAXL" : 250
              }
 # Define a function to fragment the templates  
def fragmentation_pattern(parameters, SORT_DATA = True):
    total_reads = parameters["total_reads"]
    template_length = parameters["template_length"]
    total_sample =  parameters["total_sample"]
    # Pick a random template from the samples
    template = np.random.randint(0, total_sample, total_reads) 
    # Pick a random center for the fragment (uniformly over the length)
    midpoints = np.random.randint(0, template_length, total_reads)
    # Pick a random insert length for the fragment (uniformly from 300bp to 500 bp)
    MINL = parameters["MINL"]
    MAXL = parameters["MAXL"]
    half_length = np.random.randint(MINL, MAXL, total_reads)
    # If the fragment goes off the end, make it shorter.
    starts = np.maximum(midpoints - half_length, 0)
    ends =  np.minimum(midpoints + half_length, template_length)
    # Pick orientation of read at random (r1 is forward, r2 reverse) or opposite.
    FLIP = np.random.binomial(1, 0.5, size = total_reads)
    # Sort the reads by template, start and end positions
    if SORT_DATA:
     	order = np.lexsort([ends, starts, template])
     	template = template[order]
     	starts = starts[order]
     	ends  = ends[order]
    # Output the random choices in an array
    patterns = np.array([template,starts,ends,FLIP])
    return patterns

# Generate fragmentation patterns for the mutated and unmutated templates
mutated_patterns = fragmentation_pattern(parameters)
unmutated_patterns = fragmentation_pattern(parameters)


'''
PART 5: Generate reads and write to fastq files
'''	

# Set read qualities. The MuSeq pipeline does not check read qualities
# "K" represents quality score 42
reads_qual = "K"*150

# Add read error to the sequence
def add_read_error(seq, freq=0.005):
    flag = np.random.binomial(1, freq, size = len(seq))
    seq = list(seq)
    for ind, base in enumerate(seq):
        if flag[ind] == 1:
            seq[ind] = mutate_base(base)
    seq = "".join(seq)
    return seq

# Define a function to output paired-end reads to fastq files
# folder = "mutated" or "unmutated"
def make_reads(folder, samples_seq, patterns, ADD_ERROR=True):   
    # Make a new directory for the output files
    working_dir = os.path.join(dir, folder)
    if not os.path.exists(working_dir):
    	os.mkdir(working_dir)
    # Open the read1 and read2 fastq files 
    fo1 = open(os.path.join(working_dir,"r1.fastq"), "w")
    fo2 = open(os.path.join(working_dir,"r2.fastq"), "w")
    # Set read qualities. "K" represents quality score 42
    # MuSeq pipeline does not check read qualities so they do not matter
    reads_qual = "K"*150
    # Iterate through the reads and output to fastq files
    for i in range(patterns.shape[1]):
        (temp, start, end, flip) = patterns[:,i]
        # Give the reads a unique name
     	# Include the following info: haplotype1 or haplotype2, template index, fragment start, end, flipped
        sample_label = ['haplotype1' if sample_flag[temp]==1 else 'haplotype2'][0]
        template_label = "template" + str(temp)
        start_label = "start" + str(start)
        end_label = "end" + str(end)
        read_name_short = '_'.join([sample_label,template_label,start_label,end_label])
        read_name = read_name_short + "_unflipped"
        # Get the read sequences
        # Make read1 from the start of the fragment
        read1_seq = samples_seq[temp][start : start + read_length]
        # Make read2 from the reverse complement end
        read2_rc = samples_seq[temp][end - read_length : end]
        read2_seq = str(Seq(read2_rc).reverse_complement())
        # Add read errors
        if ADD_ERROR:
            read1_seq = add_read_error(read1_seq)
            read2_seq = add_read_error(read2_seq)
        # Randomize the orientations of the reads
        if flip:
        	read_name = read_name_short + "_flipped"
        	flipped = read1_seq
        	read1_seq = read2_seq
        	read2_seq = flipped  
        # Save the read1 info in a list  
        read1_block = [read_name + "/1",
    						read1_seq,
    						"+",
    						reads_qual]
        # Save the read2 info in a list
        read2_block = [read_name+ "/2",
    						read2_seq,
    						"+",
    						reads_qual]
        # Write read1 and read2 to the fastq file
        fo1.write("\n".join(read1_block) + "\n")
        fo2.write("\n".join(read2_block) + "\n")
    fo1.close()
    fo2.close()	
    # end

# Generate paired-end fastq files for the mutated and unmutated samples
make_reads("mutated", samples_mu, mutated_patterns)
make_reads("unmutated", samples_ori, unmutated_patterns)

'''
PART 6: Zip the files
''' 
for folder in ["mutated","unmutated"]:
    working_dir = os.path.join(dir, folder)
    for file in ["r1.fastq","r2.fastq"]:
        check_call(['gzip',os.path.join(working_dir, file)])

'''
PART 7: Write the settings and parameters to a file
'''
# record 

# total length of the template molecule
# number of the repeat units
# number of bases upstream and downstream

# positions of the fake SNPs
# number of the haplotype1 and haplotype2 molecules




