#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 12:10:14 2022

@author: silvia

@description: this script changes a .wig file from fixedstep format to 
              variablestep format
"""

# =============================================================================
# Load data
# =============================================================================

fixed_wig = "out/experimental_reads.wig"
var_wig   = "out/experimental_reads_var.wig"

read_header = False
pos = None

# =============================================================================
# Convert
# =============================================================================

# input
f = open(fixed_wig, "r")

# temp output
with open(var_wig, "w") as out:
    
    for line in f.readlines():
        
        l = line.strip()
        
        # if it's a "range" line, write one line for each position
        if l[0] == "g":            
            l = l.split("\t")
            rng = range(int(float(l[1])), 
                        int(float(l[2]))+1)
            for p in rng:
                out.write(str(p) + "\t" + l[3] + "\n")
                pos = p
        
        # or if it's a header we ignore it
        elif l[0] == "f":
            pass
        
        # if it isn't just write the line preceeded by its position
        else:
            if pos == None:
                pos = 1
            else:
                pos += 1
            out.write(str(pos) + "\t" + l + "\n")
    
    f.close()