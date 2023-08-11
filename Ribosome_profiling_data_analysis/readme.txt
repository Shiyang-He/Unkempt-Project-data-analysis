# This folder contains the scripts for Ribosome Profiling data analysis
############### required softwares/packages ###############
# bowtie2 (version 2.4.5) (https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/)
# samtools (version 1.14) (http://github.com/samtools/)
# bedtools (version 2.29.2) (https://bedtools.readthedocs.io/en/latest/content/installation.html)
# UMI-tools (version 0.5.5) (https://github.com/CGATOxford/UMI-tools)

# gencode.v35.pc_transcripts.fa: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.pc_transcripts.fa.gz
# human.rRNA.combine.fa: it is generated from UCSC genome browser using Table tools: (https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1673605240_YmLhp8gHVrg69Cd5Fgsld27uG1jn)
# gencode.v35.lncRNA_transcripts.fa: it is a subset of gencode.v35.transcripts.fa: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.transcripts.fa.gz
###########################################################

# step 1: map the reads to rRNA. The rRNA sequences 
bowtie2-build human.rRNA.combine.fa human.rRNA.bowtie2.index 
bowtie2 human.rRNA.bowtie2.index -U Sample1.fq.gz -p 8 |samtools sort - -@ 8 -O BAM -o Sample1.rRNA.bam && samtools index Sample1.rRNA.bam
samtools view -f 4 -b Sample1.rRNA.bam > Sample1.rRNA.unmapped.bam && samtools index Sample1.rRNA.unmapped.bam 
bedtools bamToFastq -i Sample1.rRNA.unmapped.bam -fq Sample1.rRNA.unmapped.fq.gz 

# step 2: map the remaining reads to lncRNA
bowtie2-build gencode.v35.lncRNA_transcripts.fa hg38.lncRNA_transcripts.bowtie2.index
bowtie2 hg38.lncRNA_transcripts.bowtie2.index -U Sample1.rRNA.unmapped.fq.gz -p 8 |samtools sort - -@ 8 -O BAM -o Sample1.lncRNA.bam && samtools index Sample1.lncRNA.bam 
samtools view -f 4 -b Sample1.rRNA.bam > Sample1.lncRNA.unmapped.bam && samtools index Sample1.lncRNA.unmapped.bam 
bedtools bamToFastq -i Sample1.lncRNA.unmapped.bam -fq Sample1.lncRNA.unmapped.fq.gz 

# step 3: map the remaining reads to mRNA and remove PCR duplications
bowtie2-build gencode.v35.pc_transcripts.fa hg38.pc_transcripts.bowtie2.index 
bowtie2 hg38.pc_transcripts.bowtie2.index -U Sample1.lncRNA.unmapped.fq.gz -p 8 |samtools sort - -@ 8 -O BAM -o Sample1.mRNA.bam && samtools index Sample1.mRNA.bam 
umi_tools dedup -I Sample1.mRNA.bam -L Sample1.mRNA.dedup.log -S Sample1.mRNA.dedup.bam --extract-umi-method read_id --method unique && samtools index Sample1.mRNA.dedup.bam 

# step 4: count the reads for each mRNA
samtools view Sample1.mRNA.dedup.bam |cut -f 3 |uniq -c > Sample1.mRNA.dedup.bam.counts

# step 5: combine counts in all samples
python3 combine.read.count.py Ribo_seq.combine.counts.tsv Ribo_seq.3M_L11.mRNA.dedup.bam.counts Ribo_seq.3M_L3.mRNA.dedup.bam.counts Ribo_seq.dPAM2_L12.mRNA.dedup.bam.counts Ribo_seq.dPAM2_L4.mRNA.dedup.bam.counts Ribo_seq.Uninduced_CNOT9ko_1.mRNA.dedup.bam.counts Ribo_seq.Uninduced_CNOT9ko_2.mRNA.dedup.bam.counts Ribo_seq.Uninduced_L1.mRNA.dedup.bam.counts Ribo_seq.Uninduced_L9.mRNA.dedup.bam.counts Ribo_seq.WT_CNOT9ko_1.mRNA.dedup.bam.counts Ribo_seq.WT_CNOT9ko_2.mRNA.dedup.bam.counts Ribo_seq.WT_L10.mRNA.dedup.bam.counts Ribo_seq.WT_L2.mRNA.dedup.bam.counts

# step 6: Normalize with 311 UNK unbound genes
python3 fish.column.py unbound.list.protein_coding.fpkm.10.list,tab,1 Ribo_seq.combine.counts.tsv,tab,1 |awk '{a+=$2;b+=$3;c+=$4;d+=$5;e+=$6;f+=$7;g+=$8;h+=$9};END{print a,b,c,d,e,f,g,h}' 
    #75060 31292 39523 39300 40379 33566 16294 32616 
head -1 Ribo_seq.combine.counts.tsv > Ribo_seq.unk.unbound.normalized.xls
awk 'BEGIN{OFS="\t"};$1!="genes"{$2=sprintf("%.0f",40000*($2/75060));$3=sprintf("%.0f",40000*($3/31292));$4=sprintf("%.0f",40000*($4/39523));$5=sprintf("%.0f",40000*($5/39300));$6=sprintf("%.0f",40000*($6/40379));$7=sprintf("%.0f",40000*($7/33566));$8=sprintf("%.0f",40000*($8/16294));$9=sprintf("%.0f",40000*($9/32616));print $0}' Ribo_seq.combine.counts.tsv > Ribo_seq.normalized.counts.tsv  ## 40000 is an arbitary number close to the average of counts of the 311 genes in all conditions. 

# step 7: DEseq2 and volcano plots
Rscript deseq.R Ribo_seq.normalized.counts.tsv Ribo_seq.coldata Ribo_seq.contrast Ribo_seq
Rscript Ribo_seq.volcano.plot.R Ribo_seq.WT.vs.Uninduced.xls
Rscript Ribo_seq.volcano.plot.R Ribo_seq.3M.vs.Uninduced.xls
Rscript Ribo_seq.volcano.plot.R Ribo_seq.dPAM2.vs.Uninduced.xls

