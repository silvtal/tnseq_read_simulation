#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 13:02:33 2022

@author: silvia
"""

# =============================================================================
# load modules, data and choose excluded genes
# =============================================================================
import pandas as pd
import os

os.chdir("/home/silvia/AAA/2022-02-21_TnSeq_simul")

in_file = "F113newannot_2020.csv" # sep = ","

tab_annot = pd.read_csv(in_file, sep=",", 
                        usecols=lambda x: x in ["length", "start", "end", "strand", "productfeat"])

# I want to exclude these genes (because they are essential, they won't have 
# inserted elements in them)

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
# exclude them from the whole genome .fna
# =============================================================================
f = open("F113.fna", "r")
fna = [line.strip() for line in f.readlines()]
f.close() 

header = fna[0]
genome = fna[1]

# get the sequences we want to remove
s_e = []
for index, row in excluded.iterrows():
    s_e.append(range(row["start"]-1, row["end"])) # -1 at the start because python

# and remove them
to_include = r = set(range(len(genome)))
for r in s_e:
    to_include -= set(r)
    
genomenew = "".join([letter for idx, letter in enumerate(genome) if idx in to_include])

print("length before excluding: "+str(len(genome)))
print("length AFTER excluding: "+str(len(genomenew)))

# save file
with open("filtered_F133.fna", "w") as f:
    f.write(header+"\n")
    f.write(genomenew)
