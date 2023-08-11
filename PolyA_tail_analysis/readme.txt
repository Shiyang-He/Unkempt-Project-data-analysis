# This folder shows how poly(A) tail seq data were analyzed.
################# required softwares/packages #################
# nanopolish (version 0.14.0) (https://github.com/jts/nanopolish)
# guppy (version 6.3.7-GPU) (https://community.nanoporetech.com/downloads)
# minimap2 (version 2.26-r1175) (https://github.com/lh3/minimap2/releases)
# tailfindr (version 1.4)(https://github.com/adnaniazi/tailfindr)

# gencode.v35.transcripts.fa: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.transcripts.fa.gz (unzip it after download)
# hg38 genome: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz  (unzip it after download)
################################################################

# step 1: mapping reads to genome
minimap2 -ax splice -t 8  hg38.fa Sample1.fastq |samtools sort -o Sample1.genome.bam -T Sample1.tmp && samtools index Sample1.genome.bam

# step 2: polyA estimation with nanopolish
nanopolish index -d Sample1.fast5_pass/ Sample1.fastq
nanopolish polya --threads=8 --reads=Sample1.fastq --bam=Sample1.genome.bam --genome=hg38.fa > Sample1.nanopolish.polyA.result.tsv

# step 3: bascall with guppy 
guppy_basecaller --cpu_threads_per_caller 8 --input_path Sample1.fast5_pass/ --save_path Sample1.guppy.results --fast5_out --flowcell FLO-PRO002 --kit SQK-RNA002  --post_out

# step 4: poly(A) tail estimation with tailfindr
Rscript tailfinder.R Sample1.Tailfindr.output Sample1.guppy.results

# step 5: map reads to transcripts
minimap2 -ax map-ont -t 8 gencode.v35.transcripts.fa Sample1.fastq |samtools sort -o Sample1.transcript.bam 

# step 6: combine nanopolish and tailfindr results
python3 merge.tailfindr.and.nannopolish.results.py Sample1.transcript.bam Sample1.tailfindr.output/rna_tails.csv Sample1.nanopolish.polyA.result.tsv Sample1.RNA.joined.tsv
