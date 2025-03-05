# This folder contains a example of how iCLIP data were analyzed.
############ required tools and R packages ############
# STAR-2.7.3a (https://github.com/alexdobin/STAR/releases/tag/2.7.3a) was used for mapping. Please refer to this page (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for installation instructions. 
# samtools 1.14 (http://github.com/samtools/)
# bedtools v2.29.2 (https://bedtools.readthedocs.io/en/latest/content/installation.html)
# kentUtils (version v365) (https://github.com/ENCODE-DCC/kentUtils) 
# PureCLIP (version 1.3.0) (https://github.com/skrakau/PureCLIP)
# FASTX-Toolkit (version 0.0.14) (http://hannonlab.cshl.edu/fastx_toolkit/)
# seqtk (version 1.3) (https://github.com/lh3/seqtk)
# UMI-tools (version 0.5.5) (https://github.com/CGATOxford/UMI-tools)
# R packages: 
# GenomicFeatures (version 1.50.3) (https://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
# rtracklayer (version 1.58.0) (https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)
# GenomicRanges (version 1.50.2) (https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
# ChIPseeker (version 1.34.1) (https://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
# AnnotationDbi  (version 1.60.0) (https://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
# org.Hs.eg.db (version 3.16.0) (https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

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

# normalize bedgraph to RPM bedgraph (the number 68821122, 66264322, and  78096442 are total reads from the merged samples)
awk 'BEGIN{OFS="\t"};{$4=$4/(68821122/1000000);print $0}' 3M.merge.bam.shifted.bedgraph >3M.merge.bam.shifted.bedgraph.rpm.bedgraph
bedGraphToBigWig 3M.merge.bam.shifted.bedgraph.rpm.bedgraph hg38.fa.fai 3M.rpm.bw
awk 'BEGIN{OFS="\t"};{$4=$4/(66264322/1000000);print $0}' dPAM2.merge.bam.shifted.bedgraph >dPAM2.merge.bam.shifted.bedgraph.rpm.bedgraph 
bedGraphToBigWig dPAM2.merge.bam.shifted.bedgraph.rpm.bedgraph hg38.fa.fai dPAM2.rpm.bw
awk 'BEGIN{OFS="\t"};{$4=$4/(78096442/1000000);print $0}' WT.merge.bam.shifted.bedgraph >WT.merge.bam.shifted.bedgraph.rpm.bedgraph 
bedGraphToBigWig WT.merge.bam.shifted.bedgraph.rpm.bedgraph hg38.fa.fai WT.rpm.bw

# step 5: call peaks with pureclip
pureclip -i Sample1.dedup.bam -bai Sample1.dedup.bam.bai -g hg38.fa -ld -nt 8 -o Sample1.PureCLIP.crosslink_sites.bed -or Sample1.PureCLIP.crosslink_sites.region.bed

# step 6: peak annotation
Rscript annotation.R Sample1.PureCLIP.crosslink_sites.bed 

# draw meta-gene profile
for i in 3M.rpm.bw WT.rpm.bw dPAM2.rpm.bw; 
  do python3 get.read.count.matrix.from.whole.genenome.mapping.py $i gencode.v35.pc_transcripts.longest.fa.table gencode.v35.annotation.gtf $i.full.matrix; 
     python3 average.depth.100.bins.py.with.up.down.stream.py $i.full.matrix $i.100percent.marix
     python3 sum.matrix.400.py $i.100percent.marix $i.100percent.marix.sum
  done
Rscript draw.meta.gene.100.percent.overlay.R WT.100percent.matrix.sum 3M.100percent.matrix.sum dPAM2.100percent.matrix.sum motif.density.TAG.100.percent.matrix.sum UNK.iCLIP.metagene.plot.with.UAG 

# draw motif enrichment heatmap 
Rscript draw.motif.enrichment.heatmap.R WT.slide.table WT.fa.table 3M.slide.table 3M.fa.table dPAM2.slide.table dPAM2.fa.table Motif.enrichment.heatmap.pdf 4

# draw motif enrichment overlay plot
Rscript draw.motif.enrichment.overlay.R WT.slide.table WT.fa.table 3M.slide.table 3M.fa.table dPAM2.slide.table dPAM2.fa.table Motif.enrichment.overlay.pdf 4
