#!/bin/bash

# Automated mafft step
# Iterate over all .fasta files in the directory
for fasta_file in *.fasta; do
    # Run MAFFT on the .fasta file and produce a corresponding .mafft file
    mafft "$fasta_file" > "${fasta_file%.fasta}.mafft"
done
