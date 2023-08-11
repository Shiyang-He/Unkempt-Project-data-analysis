# This folder contains a example of how iCLIP data were analyzed.
############ required tools and R packages ############
# STAR-2.7.3a (https://github.com/alexdobin/STAR/releases/tag/2.7.3a) was used for mapping. Please refer to this page (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for installation instructions. 
# samtools 1.14 (http://github.com/samtools/)
# bedtools v2.29.2 (https://bedtools.readthedocs.io/en/latest/content/installation.html)
# kentUtils (version v365) (https://github.com/ENCODEDCC/kentUtils) 
# PureCLIP (version 1.3.0) (https://github.com/skrakau/PureCLIP)
# FASTX-Toolkit (version 0.0.14) (http://hannonlab.cshl.edu/fastx_toolkit/)
# seqtk (version 1.3) (https://github.com/lh3/seqtk)
# UMI-tools (version 0.5.5) (https://github.com/CGATOxford/UMI-tools)
# R packages: 
# GenomicFeatures (https://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
# rtracklayer  (https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)
# GenomicRanges (https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
# ChIPseeker (https://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
# clusterProfiler (https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
# AnnotationDbi (https://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
# org.Hs.eg.db (https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

# hg38 genome: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz  (unzip it after download)
# gencode V35 annotation gtf: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz (unzip it after download)
############################################


# step 1: trim and filter 
zcat flexbar_Out_barcode_Sample1.fastq.gz |fastx_trimmer -Q 33 -l 15 |fastq_quality_filter -Q 33 -q 10 -p 100 |awk 'FNR%4==1{print $1}' |sed 's/@//' >Sample1.qualFilteredIDs
seqtk subseq flexbar_Out_barcode_Sample1.fastq.gz Sample1.qualFilteredIDs |sed 's/ /#/g; s/\//#/g' |gzip > Sample1.fq.gz 

# step 2: mapping
(skip indexing if you have STAR index)
STAR --runThreadN 64 --genomeDir hg38.star.index.genecode --runMode genomeGenerate --genomeFastaFiles hg38.fa --sjdbGTFfile gencode.v35.annotation.gtf --sjdbOverhang 75  ### you may modify thread number based on your cluster node number
STAR --runMode alignReads --genomeDir hg38.star.index --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile gencode.v35.annotation.gtf --sjdbOverhang 75 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn Sample1.fq.gz --runThreadN 8 --outFileNamePrefix Sample1. --outTmpDir Sample1.tmp && samtools index Sample1.Aligned.sortedByCoord.out.bam 

# step 3: remove PCR duplicates 
umi_tools dedup -I Sample1.Aligned.sortedByCoord.out.bam -L Sample1.dedup.log -S Sample1.dedup.bam --extract-umi-method read_id --method unique && samtools index Sample1.dedup.bam

# step 4: generate bigwigfiles for crosslinkg events 
bedtools bamtobed -i Sample1.dedup.bam >Sample1.dedup.bam.bed
bedtools shift -m 1 -p -1 -i Sample1.dedup.bam.bed -g hg38.fa.fai >Sample1.dedup.bam.bed.shift.bed
bedtools genomecov -bg -5 -i Sample1.dedup.bam.bed.shift.bed -g hg38.fa.fai |sort -k1,1 -k2,2n >Sample1.bedgraph
bedGraphToBigWig Sample1.bedgraph hg38.fa.fai Sample1.bw

# step 5: call peaks with pureclip
pureclip -i Sample1.dedup.bam -bai Sample1.dedup.bam.bai -g hg38.fa -ld -nt 8 -o Sample1.PureCLIP.crosslink_sites.bed -or Sample1.PureCLIP.crosslink_sites.region.bed

# step 6: peak annotation
Rscript annotation.R Sample1.PureCLIP.crosslink_sites.bed 
