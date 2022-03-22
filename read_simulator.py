#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:55:57 2022

@author: silvia
"""

# =============================================================================
# Load modules, load genome
# =============================================================================
import os
import numpy as np

os.chdir("/home/silvia/AAA/2022-02-21_TnSeq_simul")

f = open("filtered_F113.fna", "r")
fna = [line.strip() for line in f.readlines()]
f.close()

# final headers we want: 
# '> 1077 reference=gi|378947941|ref|NC_016830|pseudocap|266 position=complement(3371764..3371793) description="[Pseudomonas fluorescens F113 chromosome, complete genome.]"'
header      = fna[0][1:]
reference   = header.split(" ",   maxsplit=1)[0]
description = header.split(" ",   maxsplit=1)[1]
genome      = fna[1]

# =============================================================================
# Find all instances of "TA" and save the 26 bp area around them
# =============================================================================
# https://stackoverflow.com/a/34445090/15477489
def findall(p, s, l=0):
    '''Yields all the positions of the pattern p in the
    string s. l is the size of the desired area. If l!=0,
    yields a tuple delimiting this area around the pattern
    position'''
    i = s.find(p)
    if l == 0:
        while i != -1:
            yield i
            i = s.find(p, i+1)
    else:
        while i != -1:
            yield (i - (l - len(p))/2, i + l/2)
            i = s.find(p, i+1)
            
mylen = 26

hits     = [i for i in findall('TA', genome, l = mylen) if (i[0] * i[1] > 0)]
# hits_rev = [i for i in findall('TA', genome[::-1], l = mylen) if (i[0] * i[1] > 0)]



# =============================================================================
# Create FASTA entries
# =============================================================================
# Randomize with desired coverage
cov = 30
pool = np.array(range(len(hits)))
choices = np.array(hits)[np.random.choice(pool, size=cov * len(hits), replace=True)]
# pool = np.array(range(len(hits_rev)))
# choices_rev = np.array(hits)[np.random.choice(pool, size=cov * len(hits), replace=True)]

# save file
with open("experimental_reads.fa", "w") as f:
    for n, c in enumerate(choices):
        start = int(c[0]); end = int(c[1])
        f.write("> " + str(n) + \
                " reference=" + header + \
                " position=" + str(start) + ".." + str(end) + \
                " description=" + description + \
                "\n")
        f.write(genome[start:end] + "\n")