suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(stringr))

argv <- commandArgs(TRUE)
txdb=''
if (grep("\\.gtf$",argv[2])){
	annotation=import(argv[2],format="GTF")
	txdb=makeTxDbFromGRanges(annotation)
}else if(grep("\\.sqlite$",argv[2])){
	txdb=loadDb(argv[2])
}else{
	message("try to load your \n",argv[2], " \nas a TxDb.sqlite object, if it's not generated from makeTxDbFromGRanges function, please provide a gtf file and retry")
	txdb=loadDb(argv[2])
}

peakAnno <- annotatePeak(argv[1], tssRegion=c(-150, 10),TxDb=txdb, annoDb="org.Hs.eg.db")
pA=as.data.frame(peakAnno)
anno.db=basename(argv[2])
anno.db=str_replace(anno.db,".TxDb.sqlite","")
anno.db=str_replace_all(anno.db,".gtf","")
write.table(pA,file=paste(argv[1],anno.db,"annotation.xls",sep="."),sep="\t",quote=F,col.names=T,row.names=F)
