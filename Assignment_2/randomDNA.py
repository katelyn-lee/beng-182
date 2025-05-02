# This program generates n random DNA sequence(s) (where of A=G=C=T=0.25) with a specified length
# Call the programm with the following arguments: python3 randomDNA.py <number of seq> <length of seq> 
# Output: random sequences, one per line. followed by a summary of the observed nucleotide frequencies in the set of sequences

import sys
import subprocess
import random

# initialize variables 
num_seqs = 0
seq_length = 0
forHistogram = False

if len(sys.argv) == 3:
    num_seqs = int(sys.argv[1])
    seq_length = int(sys.argv[2])
    print(f"Correct args input for randomDNA.py. Running program...")
elif len(sys.argv) == 4:
    num_seqs = int(sys.argv[1])
    seq_length = int(sys.argv[2])
    forHistogram = True
else:
    print(f"Incorrect args input. Try again with format: python randomDNA.py <number of seq> <length of seq>.")
    num_seqs = 0
    seq_length = 0

# generates collection of random sequences, and generates frequency at the same time
def random_dna_seq(num_seqs:int, seq_length:int):
    base_freq = {'A':0, 'C':0, 'G':0, 'T':0}
    nucleotides = ['A','T','C','G']
    seq_collection = []
    print(f"Generated {num_seqs} random DNA sequences of length {seq_length}:")
    for i in range(int(num_seqs)):
        seq = ''
        for j in range(int(seq_length)):
            rand_base = random.choice(nucleotides)
            seq = seq + rand_base
            if rand_base == 'A':
                base_freq['A'] = base_freq['A'] + 1
            if rand_base == 'C':
                base_freq['C'] = base_freq['C'] + 1
            if rand_base == 'G':
                base_freq['G'] = base_freq['G'] + 1
            if rand_base == 'T':
                base_freq['T'] = base_freq['T'] + 1
        seq_collection.append(seq)
        print(seq)
    return seq_collection, base_freq

def random_dna_seq_noFreq(num_seqs:int, seq_length:int):  #for internal use for other question
    nucleotides = ['A','T','C','G']
    seq_pair = ""
    for i in range(int(num_seqs)):  #num seqs will always be 2 to get a pair
        seq = ''
        for j in range(int(seq_length)):
            rand_base = random.choice(nucleotides)
            seq = seq + rand_base
        seq_pair = seq_pair + seq + '\n'
    return seq_pair

def frequencies(base_freq: dict[str, int]):
    total_bases = num_seqs*seq_length
    print(f"Frequencies per DNA base: A {(base_freq['A'])/total_bases}%, C {(base_freq['C'])/total_bases}%, G {(base_freq['G'])/total_bases}%, T {(base_freq['T'])/total_bases}%")

if forHistogram==False:
    seq_collection, base_freq = random_dna_seq(num_seqs,seq_length)
    frequencies(base_freq)
else:
    # just for generating random seq collection, no need for frequency for histogram question
    print(random_dna_seq_noFreq(num_seqs,seq_length))
