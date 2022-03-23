#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:55:57 2022

@author: silvia
"""

# =============================================================================
# load modules, reference genome and annotations
# =============================================================================
import pandas as pd
import os
import numpy as np

os.chdir("/home/silvia/AAA/2022-02-21_TnSeq_simul")

in_file = "F113newannot_2020.csv" # sep = ","

tab_annot = pd.read_csv(in_file, sep=",", 
                        usecols=lambda x: x in ["length", "start", "end", "strand", "productfeat"])

f = open("F113.fna", "r")
fna = [line.strip() for line in f.readlines()]
f.close()

# final headers we want: 
# '> 1077 reference=gi|378947941|ref|NC_016830|pseudocap|266 position=complement(3371764..3371793) description="[Pseudomonas fluorescens F113 chromosome, complete genome.]"'
header      = fna[0][1:]
reference   = header.split(" ",   maxsplit=1)[0]
description = header.split(" ",   maxsplit=1)[1]
genome      = fna[1]


# =============================================================================
# choose excluded genes
# =============================================================================
lst = ["Motility protein A",
       "Motility protein B",
       "Methyl-accepting chemotaxis protein CtpL",
       "Methyl-accepting chemotaxis protein McpS",
       "Methyl-accepting chemotaxis protein McpU",
       "Methyl-accepting chemotaxis protein McpP",
       "Methyl-accepting chemotaxis protein PctA",
       "Chemotaxis protein CheY",
       "Chemotaxis protein CheW"]


excluded = tab_annot.query('productfeat in @lst')

print("Genes to exclude")
print("-----------------")
print(excluded)


# =============================================================================
# save the included ranges
# =============================================================================
s_e = []
for index, row in excluded.iterrows():
    s_e.append(range(row["start"]-1, row["end"])) # -1 at the start because python

to_include = r = set(range(len(genome)))
for r in s_e:
    to_include -= set(r)
    
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

# select valid ones
hits     = [i for i in findall('TA', genome, l = mylen) 
            if ((i[0] * i[1] > 0) and (i[1]-1 in to_include))]
hits_rev = [i for i in findall('AT', genome, l = mylen) 
            if ((i[0] * i[1] > 0) and (i[1]-1 in to_include))]


# =============================================================================
# Create FASTA entries
# =============================================================================
# Randomize with desired coverage
cov = 30
pool = np.array(range(len(hits)))
choices = np.array(hits)[np.random.choice(pool, size=int(cov/2) * len(hits), replace=True)]
pool = np.array(range(len(hits_rev)))
choices_rev = np.array(hits_rev)[np.random.choice(pool, size=int(cov/2) * len(hits_rev), replace=True)]

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
        
    for n, c in enumerate(choices_rev):
        start = int(c[0]); end = int(c[1])
        f.write("> " + str(n) + \
                " reference=" + header + \
                " position=complement(" + str(start) + ".." + str(end) + ")" \
                " description=" + description + \
                "\n")
        f.write(genome[start:end] + "\n")