#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 10:21:16 2024

@author: Zelia
"""
## modules
import pandas as pd

# Load data from files, removing the header from the TSV file
besthits_df = pd.read_csv('besthits_id.txt', header=None, names=['Protein_Name'])
annotation_df = pd.read_csv('Fisci.tsv', sep='\t', header=None, usecols=[0, 6], names=['Protein_Name', 'Annotation'])

# Merge the two DataFrames on the 'Protein_Name' column
merged_df = pd.merge(besthits_df, annotation_df, on='Protein_Name', how='left')

# Save the result to a new file
merged_df.to_csv('besthits_with_annotation.txt', sep='\t', index=False)
