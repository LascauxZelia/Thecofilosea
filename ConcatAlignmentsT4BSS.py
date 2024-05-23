#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 17:40:18 2022

@author: nina
"""
## modules
import pandas as pd
import glob
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil

### POSITION OF ORGANISMS IN ALIGNMENT ###
with open('/home/zelia/Desktop/Nanopore/T4BSS/concat_aln/organisms') as oL:
    OrgL = [line.rstrip("\n") for line in oL.readlines()]
print(OrgL)
######

### ORGANISM DICTIONARY ###
oL_df = pd.read_csv('/home/zelia/Desktop/Nanopore/T4BSS/concat_aln/organismList',sep='\t', header = None)
print(oL_df)

organismList = oL_df[0].tolist()
organismList = list(dict.fromkeys(organismList))
#print(organismList)

R64_values = []
LegPhi_values = []
LegLon_values = []
LegOak_values = []
Tara121_values = []
AquLus_values = []
RicIso_values = []
CoxRSA93_values = []
BerAqu_values = []
PisSal_A_values = []
PisSal_B_values = []
PisSal_C_values = []
FanHon_A_values = []
FanHon_B_values = []
FanHon_D_values = []
AciFer_values = []
Poke_values = []
Fisci_values = []

for index, row in oL_df.iterrows():
    if "R64" in row[0]:
        R64_values.append(row[4])
    if "LegPhi" in row[0]:
        LegPhi_values.append(row[4])
    if "LegLon" in row[0]:
        LegLon_values.append(row[4])
    if "LegOak" in row[0]:
        LegOak_values.append(row[4])
    if "Tara121" in row[0]:
        Tara121_values.append(row[4])
    if "AquLus" in row[0]:
        AquLus_values.append(row[4])
    if "RicIso" in row[0]:
        RicIso_values.append(row[4])
    if "CoxRSA93" in row[0]:
        CoxRSA93_values.append(row[4])
    if "BerAqu" in row[0]:
        BerAqu_values.append(row[4])
    if "PisSal_A" in row[0]:
        PisSal_A_values.append(row[4])
    if "PisSal_B" in row[0]:
        PisSal_B_values.append(row[4])
    if "PisSal_C" in row[0]:
        PisSal_C_values.append(row[4])
    if "FanHon_A" in row[0]:
        FanHon_A_values.append(row[4])
    if "FanHon_B" in row[0]:
        FanHon_B_values.append(row[4])
    if "FanHon_D" in row[0]:
        FanHon_D_values.append(row[4])
    if "AciFer" in row[0]:
        AciFer_values.append(row[4])
    if "Poke" in row[0]:
        Poke_values.append(row[4])
    if "Fisci" in row[0]:
        Fisci_values.append(row[4])
print(Fisci_values)

valuesList = []
def concats(lists):
    for i in lists:
        valuesList==valuesList.append(i)

concats([R64_values, LegPhi_values,LegLon_values, LegOak_values, Tara121_values, AquLus_values, RicIso_values, CoxRSA93_values, BerAqu_values, PisSal_A_values, PisSal_B_values, PisSal_C_values, FanHon_A_values, FanHon_B_values, FanHon_D_values, AciFer_values, Poke_values, Fisci_values,])
print(valuesList)

organsim_dict = {organismList[i]: valuesList[i] for i in range(len(organismList))}
print(organsim_dict)
######

### DOWOH ###
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def Diff(li1, li2):
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif

dowoh = {'dotA':[],'dotB':[],'dotC':[],'dotD':[],'icmB':[],'icmD':[],'icmH':[],'icmK':[],'icmO':[],'icmQ':[],'icmS':[],'icmV':[],'icmC':[],'icmE':[],'icmG':[],'icmJ':[],'icmL':[],'icmN':[],'icmP':[],'icmT':[],'icmW':[],'icmX':[],'icmR':[],'icmF':[],'icmM':[]}

dotA = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/dotA.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in dotA])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['dotA'] = Diff(organismList, org_flat)

dotB = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/dotB.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in dotB])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['dotB'] = Diff(organismList, org_flat)

dotC = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/dotC.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in dotC])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['dotC'] = Diff(organismList, org_flat)

dotD = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/dotD.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in dotD])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['dotD'] = Diff(organismList, org_flat)

icmB = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmB.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmB])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmB'] = Diff(organismList, org_flat)

icmC = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmC.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmC])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmC'] = Diff(organismList, org_flat)

icmD = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmD.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmD])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmD'] = Diff(organismList, org_flat)

icmE = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmE.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmE])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmE'] = Diff(organismList, org_flat)

icmF = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmF.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmF])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmF'] = Diff(organismList, org_flat)

icmG = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmG.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmG])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmG'] = Diff(organismList, org_flat)

icmH = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmH.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmH])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmH'] = Diff(organismList, org_flat)

icmJ = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmJ.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmJ])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmJ'] = Diff(organismList, org_flat)

icmK = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmK.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmK])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmK'] = Diff(organismList, org_flat)

icmL = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmL.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmL])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmL'] = Diff(organismList, org_flat)

icmM = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmM.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmM])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmM'] = Diff(organismList, org_flat)

icmN = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmN.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmN])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmN'] = Diff(organismList, org_flat)

icmO = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmO.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmO])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmO'] = Diff(organismList, org_flat)

icmP = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmP.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmP])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmP'] = Diff(organismList, org_flat)

icmQ = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmQ.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmQ])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmQ'] = Diff(organismList, org_flat)

icmR = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmR.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmR])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmR'] = Diff(organismList, org_flat)

icmS = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmS.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmS])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmS'] = Diff(organismList, org_flat)

icmT = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmT.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmT])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmT'] = Diff(organismList, org_flat)

icmV = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmV.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmV])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmV'] = Diff(organismList, org_flat)

icmW = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmW.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmW])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmW'] = Diff(organismList, org_flat)

icmX = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/mafft/icmX.mafft","fasta"))
org = []
for key in organsim_dict:
    intersect = intersection(organsim_dict[key], [rec.id for rec in icmX])
    for el in intersect:
        if el in organsim_dict[key]:
            org.append([k for k, v in organsim_dict.items() if v == organsim_dict[key]])
org_flat = [x for xs in org for x in xs]
dowoh['icmX'] = Diff(organismList, org_flat)
####

### ADDING GAP SEQUENCES ###
def seqlen(align):
    length = list(set([len(seq) for seq in align]))
    return length

files = glob.glob("/home/zelia/Desktop/Nanopore/T4BSS/mafft/*.mafft")
for f in files:
    T4BSS_seq = list(SeqIO.parse(f,"fasta"))
    T4BSS_len = seqlen(T4BSS_seq)
    for key in dowoh:
        if f[50:54] in key:
            if not dowoh[key]:
                name = f.split('/')[-1].split('.')[0]
                print(name)
                SeqIO.write(T4BSS_seq, name + "_ed.mafft", "fasta")
            for org in dowoh[key]:
                for i in OrgL:
                    if org in i:
                        gapseq = SeqRecord(Seq("-" * T4BSS_len[0]), id=i)
                        T4BSS_seq_list = []
                        T4BSS_seq_list.extend([T4BSS_seq[:OrgL.index(i)], T4BSS_seq[OrgL.index(i):]])
                        T4BSS_seq_list[0].append(gapseq)
                        T4BSS_seq = [x for xs in T4BSS_seq_list for x in xs]
                        name = f.split('/')[-1].split('.')[0]
                        SeqIO.write(T4BSS_seq, name + "_ed.mafft", "fasta")
                    else:
                        name = f.split('/')[-1].split('.')[0]
                        SeqIO.write(T4BSS_seq, name + "_ed.mafft", "fasta")
####

### COPY EDITED FILES TO MAFFT_EDITED ###
dst_path = "/home/zelia/Desktop/Nanopore/T4BSS/mafft_edited"
for file in glob.glob("/home/zelia/Desktop/Nanopore/T4BSS/concat_aln/*_ed.mafft"):
    shutil.move(file, dst_path)
####

### CONCATENATE ALIGNMENTS ###
MSAList = []
files = glob.glob("/home/zelia/Desktop/Nanopore/T4BSS/mafft_edited/*")
for f in files:
    MSAList.append(AlignIO.read(f, "fasta"))

concat_align = MSAList[0]
for aln in MSAList:
    concat_align = concat_align + aln
AlignIO.write(concat_align, "concatAlign.fasta", "fasta")
####

### ADDING HEADERS ###
concatAlign = list(SeqIO.parse("/home/zelia/Desktop/Nanopore/T4BSS/concat_aln/concatAlign.fasta", "fasta"))
print(str(len(concatAlign)))
orgindex = list(range(0,17))
orgindex_dict = {orgindex[i]: organismList[i] for i in range(len(orgindex))}

for index in orgindex_dict:
    concatAlign[index].id = orgindex_dict[index]
    concatAlign[index].name = orgindex_dict[index]
    concatAlign[index].description = ""
SeqIO.write(concatAlign, "concatAlign.fasta", "fasta")
####

# Vérifier l'ordre des séquences concaténées
for index, seq_record in enumerate(concatAlign):
    organism = organismList[index]
    if seq_record.id != organism:
        print(f"Erreur: L'en-tête {seq_record.id} ne correspond pas à l'organisme {organism} dans l'ordre.")

# Vérifier l'attribution des en-têtes
for index, seq_record in enumerate(concatAlign):
    if seq_record.id != organismList[index]:
        print(f"Erreur: L'en-tête {seq_record.id} ne correspond pas à l'organisme {organismList[index]}.")

