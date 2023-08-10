# This folder demonstrates how RNA_seq data was analyzed
########## prerequisite ##########
# software/packages: 
# STAR-2.7.3a (https://github.com/alexdobin/STAR/releases/tag/2.7.3a) was used for mapping. Please refer to this page (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for installation instructions. 
# samtools 1.14 (http://github.com/samtools/)
# featureCounts v2.0.0 (https://subread.sourceforge.net/) 
# cufflinks v2.2.1 (http://cole-trapnell-lab.github.io/cufflinks/install/)
# DEseq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
# ggplot2 (https://www.rdocumentation.org/packages/ggplot2/versions/3.4.2)
########## prerequisite ##########

# step 1: building STAR mapping indexes
STAR --runThreadN 64 --genomeDir hg38.star.index.genecode --runMode genomeGenerate --genomeFastaFiles hg38.fa --sjdbGTFfile gencode.v35.annotation.gtf --sjdbOverhang 75 

#### step 2,3,4 can be done in parallel if you have multiple nodes. Here I only write one example (Sample1). In real analysis, there were 12 samples (3M-1, 3M-2, 3M-3, dPAM2-1, dPAM2-2, dPAM2-3, Uninduced-1, Uninduced-2, Uninduced-3, WT-1, WT-2, WT-3)
  # step 2: mapping
  STAR --genomeDir hg38.star.index.gencode --readFilesIn raw.data/Sample1_1.fq.gz raw.data/Sample1_2.fq.gz --outFileNamePrefix samplename1 --runThreadN 8 --limitBAMsortRAM 60000000000 --outBAMsortingThreadN 8 --outSAMattrRGline ID:Sample1 SM:Sample1 --outSAMtype BAM SortedByCoordinate -- 
outSAMunmapped Within --outSAMstrandField intronMotif --readFilesCommand zcat --chimSegmentMin 20 --genomeLoad NoSharedMemory && samtools index Sample1.Aligned.sortedByCoord.out.bam

  # step 3: count features for all sample 
  featureCounts -a gencode.v35.annotation.gtf -o Sample1.featureCounts -M -C -T 8 Sample1.Aligned.sortedByCoord.out.bam

  # step 4: compute FPKM
  cufflinks -u -p 8 -o Sample1.cufflinks --max-bundle-frags 1000000000 --no-effective-length-correction --compatible-hits-norm -G gencode.v35.annotation.gtf Sample1

# step 5: combine feature counts from all sample
python3 combine.featureCounts.py UNK.RNAseq.featureCounts.combine.xls 3M-1.featureCounts 3M-2.featureCounts 3M-3.featureCounts dPAM2-1.featureCounts dPAM2-2.featureCounts dPAM2-3.featureCounts Uninduced-1.featureCounts Uninduced-2.featureCounts Uninduced-3.featureCounts WT-1.featureCounts WT-2.featureCounts WT-3.featureCounts

# step 6: normalize counts with 311 UNK unbound genes (These genes has no iCLIP peaks in any of UNK_WT, UNK_3M, and UNK_dPAM2 but has high expression level (FPKM >=10))
python3 fish.column.py unbound.list.protein_coding.fpkm.10.list,tab,1 UNK.RNAseq.featureCounts.combine.xls,tab,1 |awk '{a+=$2;b+=$3;c+=$4;d+=$5;e+=$6;f+=$7;g+=$8;h+=$9;i+=$10;j+=$11;k+=$12;l+=$13};END{print a,b,c,d,e,f,g,h,i,j,k,l}' >UNK.RNAseq.featureCounts.combine.xls.unk.unbound.xls.read.count
#266561 190328 233278 226218 282778 230812 245494 234755 243683 274887 260257 300119
head -1 UNK.RNAseq.featureCounts.combine.xls >UNK.RNAseq.featureCounts.combine.xls.unk.unbound.normalized.xls
awk 'BEGIN{OFS="\t"};$1!="genes"{$2=sprintf("%.0f",250000*($2/266561));$3=sprintf("%.0f",250000*($3/190328));$4=sprintf("%.0f",250000*($4/233278));$5=sprintf("%.0f",250000*($5/226218));$6=sprintf("%.0f",250000*($6/282778));$7=sprintf("%.0f",250000*($7/230812));$8=sprintf("%.0f",250000*($8/245494));$9=sprintf("%.0f",250000*($9/234755));$10=sprintf("%.0f",250000*($10/243683));$11=sprintf("%.0f",250000*($11/274887));$12=sprintf("%.0f",250000*($12/260257));$13=sprintf("%.0f",250000*($13/300119));print $0}' UNK.RNAseq.featureCounts.combine.xls >>UNK.RNAseq.featureCounts.combine.xls.unk.unbound.normalized.xls  # 250000 is a arbitary number close to the average of UNK unbound gene counts. 

# step 7: perform DEG analysis with DEseq2
Rscript deseq.R UNK.RNAseq.featureCounts.combine.xls.unk.unbound.normalized.xls coldata contrast UNK.RNAseq 

# step 8: draw volcanoplots (this step requries iCLIP result file: UNK.iCLIP.gene.counts.xls.combined.xls.classified.xls)
Rscript volcano.R UNK.RNAseq.WT.vs.Uninduced.xls
Rscript volcano.R UNK.RNAseq.3M.vs.Uninduced.xls
Rscript volcano.R UNK.RNAseq.dPAM2.vs.Uninduced.xls


  





