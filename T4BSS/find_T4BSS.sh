## Find T4BSS genes
# Make Legionella pneumophila db
makeblastdb -in GCA_001753085.1_ASM175308v1_genomic.fna -dbtype nucl -out Lpn_db
#Comparision with Legionella pneumophila genome
tblastx -query 3C_reordered.fasta -db Lpn_DB/Lpn_db -out results.txt -evalue 1e-5 -outfmt 6
# Visualization using ACT
act 3_reordered_prokka/3C_reordere_prokka.gbk results.txt Lpn_DB/Lpn_prokka/Lpn_prokka.gbk 
