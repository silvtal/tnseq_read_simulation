# tnseq_read_simulation
Reads simulation and analysis following the alternative method proposed at https://journals.asm.org/doi/full/10.1128/mSystems.00976-20# (Fig 1, right)

## Brief description

### read_simulator.py
Samples short sequences with a given coverage and read length from a given genome .fna file. Gives the same .fasta output that [Grinder](https://github.com/zyxue/biogrinder) would do with "random primers", but this script has the additional option of only generating reads that include a desired subsequence in the middle. In our case, it's a "TA" sequence, so the simulated reads are like this

`NNNNNNNNNNNNTANNNNNNNNNNNN`

Notably, this script takes a list of "genes to exclude" (to simulate essentiality). If the base pairs after "TA" are included in the "genes to exclude", no reads will be generated from that "TA".

The essentiality analysis can be carried out by aligning with a short read aligner, formatting to .wig and using [Transit](https://github.com/mad-lab/transit). 

### rm_selected_genes.py
Removes specific genes from a .fna file given a .tsv annotation file. In our example, we remove chemotaxis and motility genes. The filtered genome will serve as a reference for visualization on software like IGV.
