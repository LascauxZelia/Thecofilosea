#Don't forget to activate conda env (conda activate Thecofilosea)

#### ASSEMBLY ####
#1. run fastQC on fastq.gz
fastqc -t 14 -o /home/zelia/Desktop/Nanopore/fastqc/ *.fastq.gz 
#2. run porechop on fastq.gz
porechop -i 3C.fastq.gz -t 14 -v 2 -o 3C_trimmed.fastq.gz |tee porechop.log
#3. Re-run fastQC on trimmed sequences
fastqc -t 14 -o /home/zelia/Desktop/Nanopore/trim/fastqc/ *.fastq.gz 
#4. filtlong
filtlong --min_length 1000 --keep_percent 95 3C_trimmed.fastq.gz | gzip > 3C_trimmed_q95.fastq.gz
# 5. Rakakost (corrected long read with short reads)
Ratatosk correct -v -G -c 16 -s ~/Desktop/Nanopore/Illumina_3C/3C_MiSeq.fastq -l 3C_trimmed_q95.fastq -o ratatosk_reads

# 6.Assembly using  Flye (nano-raw) 
flye --nano-raw /home/zelia/Desktop/Nanopore/trim/ratatosk_reads.fastq.gz --o 3C_assembly --thread 32 --meta
# 7. Run blast on 16S endosymbiont db (From the German lab)
#Before I have to make a db. Reference 16S + 5 best hit on blast and then:
makeblastdb -in 3C_16S_2023.fasta -dbtype nucl
blastn -query assembly.fasta -db /home/zelia/Desktop/Nanopore/DB/3C/3C_16S_2023.fasta -out blast_3C.txt
# 8. Extract contig with seqtk
seqtk subseq assembly.fasta subsetIDs.txt > contigs_subset.fasta
# 9. Run prokka on extracted contigs
prokka contig_subset.fasta -outdir prokka --prefix 3C_prokka
# 10. Map contigs with 16S RNA againt filtered reads from the start
minimap2 -t 14 -x map-ont -a -o /home/zelia/Desktop/Nanopore/mixed_reads/3C_all_assembly/mapping/3C_16S_mapping.sam /home/zelia/Desktop/Nanopore/mixed_reads/3C_all_assembly/contigs_subset.fasta /home/zelia/Desktop/Nanopore/mixed_reads/3C_all_trimmed_q95.fastq.gz 
# 11. Map contigs with 16S RNA against Illumina MiSeq reads
#Create an index with
bwa index contigs_subset.fasta
#And then execute bwa mem
bwa mem -t 14 /home/zelia/Desktop/Nanopore/mixed_reads/3C_all_assembly/mapping/contigs_subset.fasta /home/zelia/Desktop/Nanopore/old_reads/3C_MiSeq_R1.fastq.gz > 3C_16S_mapping.sam
11. Convert sam files to sorted and indexed bam files (samtools) - samtools works only in new conda env without other packages even isn't compatible with libopensl
Execute for ONP and Illumina
samtools view -b 3C_16S_mapping.sam | samtools sort - > 3C_16S_mapping.bam
samtools view -b 3C_Illumina_16S_mapping.sam | samtools sort - > 3C_Illumina_16S_mapping.bam
12. Then extract only mapped sequences from BAM to FASTA format
samtools fasta -F 4 3C_16S_mapping.bam > 3C_16S_mapping.fasta
samtools fasta -F 4 3C_Illumina_16S_mapping.bam > 3C_Illumina_16S_mapping.fasta
Want to transform in fastq ? 
seqtk seq -A 3C_16S_mapping.fasta > 3C_16S_mapping.fastq
13. cat mapping and contig subset fasta file delete step
cat 3C_16S_mapping.fasta contigs_subset.fasta > 3C_16S_final.fasta

# 16. Run miComplete
find -maxdepth 1 -type f -name "*.fna" | miCompletelist.sh > test_set.tab
miComplete test_set.tab --hmms Bact105 --weights Bact105 --threads 14 --hlist list > results.tab
# 17. Run barnnap, extract 16S and blast 16S
barrnap -o contig_1456_rrna.fa contig_1456.fasta
# 19. Binning metawrap
metawrap binning -o Initial_binnig -t 32 -a assembly.fasta --metabat2 --maxbin2 --concoct --single-end 3C_all_trimmed_q95.fastq
metawrap bin_refinement -o bin_refinement -t 32 -A Initial_binnig/metabat2_bins/ -B Initial_binnig/maxbin2_bins/ -C Initial_binnig/concoct_bins/ -c 50 -x 10
# 20. Medaka
medaka_consensus -i 3C_trimmed_q95.fastq -d 3C_assembly/assembly.fasta -o medaka_3C -t 4 -m r941_min_high_g303
# 6. Calculate the average of the CDs length for the prokka output
awk 'NR > 1 {sum += $3} END {print sum / (NR - 1)}' bin5.tsv
# 7. Calculate the sum of genome that is coding (from prokka output)
awk 'NR > 1 {sum += $3} END {print sum}' bin5.tsv



## RE-doing the assembly process pour Pokemonad avec Rataosk + mapping + medaka and see
Ratatosk correct -v -G -c 16 -s B5_MiSeq.fastq -l 20210826_B5_porechop_filtlong95.fastq -o ratatosk_B5_reads
flye --nano-raw ratatosk_B5_reads.fastq.gz --o B5_ratatosk_assembly --thread 14 --meta
gunzip ratatosk_B5_reads.fastq.gz
metawrap binning -o binning_ratatosk -t 12 -a assembly.fasta --concoct --single-end ratatosk_B5_reads.fastq