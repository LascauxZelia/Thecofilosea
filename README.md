# Cell fusion facilitates horizontal transmission of Legionellales independent of an environmental life cycle stage.  

Zélia Bontemps(1,†), Nina Pohl(1,2,†), Herve Nicoloff(1), Susanne Schlegel(3), Kenneth Dumack(4,†), Lionel Guy(1,†,*)  
 
- (1) Department of Medical Biochemistry and Microbiology, Science for Life Laboratory, Uppsala University, Uppsala, Sweden
- (2) Department of Organismal Biology, Program in Systematic Biology, Uppsala University, 752 36 Uppsala, Sweden
- (3) Department of Cell and Molecular Biology, Science for Life Laboratory, Uppsala University, P.O. Box 596, 75123 Uppsala, Sweden
- (4) Terrestrial Ecology, Institute for Zoology, University of Cologne, Zülpicher Str. 47b, 50674, Cologne, Germany
† These authors contributed equally to this work.

[DOI]()

## Abstract
Understanding the mechanisms of replication and transmission of endosymbiotic pathogens is crucial for human health. While Legionella pneumophila and its relatives all exhibit an obligatory endosymbiotic life history stage in which they replicate, horizontal transmission is a common and important dispersal mechanism. The common hosts for Legionellales are amoebae, where they can multiply, shed into the environment, and then infect a new host via a molecular machinery called the Type 4B Secretion System (T4BSS). Recently, Legionellales have been discovered in several species of the rhizarian class Thecofilosea - a taxon entirely unrelated to previously known hosts of Legionellales.  To characterize these novel endosymbionts of an evolutionary novel host group, genomes of 'Ca. Pokemonas kadabra' and 'Ca. Fiscibacter pecunius' were sequenced. We show that these two endosymbionts have a fragmented and reduced T4BSS and the relationship with their host is still unclear. However, the missing flagellar genes and the fragmented and reduced T4BSS let us investigate their dispersal strategies. We hypothesized and confirmed that next to vertical transmission, these Legionellales transmit horizontally during the frequent fusion of their host cells. Horizontal transmission during cell fusion is a novel dispersal mechanism for Legionellales. 

#### Keywords: Endosymbionts; Horizontal transmission; Legionellales; Mutualism-parasitism continuum; T4BSS

## Scripts
* Genome assembly steps (assembly.sh)
* Genome size (genome_size.R)
* T4BSS (find_T4BSS.sh, extract_CDS.py, mappintab_T4BSS.py, mafft_loop.sh, ConcatAlignementsT4BSS.py, plotMap.R)
* Effector (BestHits.py, add_annotation.py)
* Tree (create_mapping_tab.py, concatenateRenameAlignment.pl)

## R version
R 4.3.1
## Python version
Python 3.10.13

## R packages
* genoPlotR
* ggplot2
* dplyr
* readr
* ggrepel
* gridExtra

## Python modules
* AlignIO
* SeqIO
* Seq
* SeqRecord
* shutil
* pandas
* glob
* os
* Path

## Bioinformatics tools
* Filtlong
* Porechop
* FastQC
* Flye
* Bwa
* Samtools
* Ratatosk
* miComplete
* barrnap
* seqtk
* blast (psiblast, tblastx)
* mafft
* trimAl
* iqtree2
* ACT
