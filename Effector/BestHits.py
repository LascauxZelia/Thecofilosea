#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 15:02:52 2024

@author: Zelia
"""
## modules
from pathlib import Path
import pandas as pd
import glob

# Function to extract best hits from pisblast file
def extract_best_hits(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t', header=None, usecols=[0, 1, 2, 10, 11])
        df.columns = ['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore']
        df_sorted = df.sort_values('evalue')
        best_hit = df_sorted.iloc[0]
        return best_hit
    except pd.errors.EmptyDataError:
        print(f"Le fichier {file_path} est vide.")
        return None

# Psiblast files directories path
psiblast_dir = "/home/zelia/Desktop/Nanopore/flagella/psiblast/"
output_file = "/home/zelia/Desktop/Nanopore/flagella/psiblast/besthits.tab"

# Initialization of a dictionnary to store the best hits
best_hits_dict = {'effector':[], 'qseqid':[], 'sseqid':[], 'pident':[], 'evalue':[], 'bitscore':[]}

# Read all the pisblast files
for file_path in glob.glob(f"{psiblast_dir}/*.psiblast"):
    best_hit = extract_best_hits(file_path)
    if best_hit is not None:
        effector_name = Path(file_path).stem
        best_hits_dict['effector'].append(effector_name)
        for key in ['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore']:
            best_hits_dict[key].append(best_hit[key])

# Dataframe creation from the dictionnary
best_hits_df = pd.DataFrame.from_dict(best_hits_dict)

# Writing the df in the output file
best_hits_df.to_csv(output_file, sep='\t', index=False)
