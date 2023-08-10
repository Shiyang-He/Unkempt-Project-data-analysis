argv <- commandArgs(TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))

combine=read.table(argv[1],header=T,row.names=1)
sample_name=colnames(combine)
coldata=data.frame(read.table(argv[2],header=T,row.names=1))
contrast=read.table(argv[3],header=T)
prefix=argv[4]
normFactors=matrix(rep(rep(1,nrow(combine)),ncol(combine)),nrow=nrow(combine),ncol=ncol(combine))
colnames(normFactors)=c(1:ncol(combine))
rownames(normFactors)=c(1:nrow(combine))

dds <- DESeqDataSetFromMatrix(countData = combine,colData = coldata,design = ~ condition)
normalizationFactors(dds) <- normFactors
dds <- DESeq(dds)
normalized_counts <- round(counts(dds, normalized=TRUE))
for (i in rownames(contrast)){
	c=as.character(contrast[i,]$control)
	t=as.character(contrast[i,]$treat)
	res <- data.frame(results(dds,contrast=c("condition",t,c)));
	colnames(res)=c("baseMean",paste("log2FC",t,"vs",c,sep="."),"lfcSE","stat","pvalue","padj")
	sink(paste(prefix,t,"vs",c,"xls",sep="."));cat("gene\t");sink()
	write.table(cbind(normalized_counts,res),paste(prefix,t,"vs",c,"xls",sep="."),row.names=T,quote=F,sep="\t",append=T)
}
rld <- rlog(dds, blind=FALSE)
png(paste(prefix,".all.samples.heatmap.rld.png",sep=""),width=900,height=900)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:2000]

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE)
dev.off()

pdf(paste(prefix,".all.samples.PCA.pdf",sep=""),width=9)
pca=plotPCA(rlog(dds),returnData=TRUE)
percent=round(100*attr(pca,"percentVar"))
ggplot(pca,aes(x = PC1, y = PC2,color=condition,label=pca$name))+
	geom_point(size=2)+
       # scale_color_manual(values=colors)+
        xlab(paste0("PC1: ",percent[1],"% of Variance"))+
	geom_text(size=3,hjust=0.5,vjust=1)+
        ylab(paste0("PC2: ",percent[2],"% of Variance"))+
	ggtitle("PCA of all samples")+
	theme(plot.title = element_text(hjust = 0.5),
		panel.background = element_blank(),
		panel.border=element_rect(fill=NA),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		strip.background=element_blank(),
		axis.text.x=element_text(colour="black"),
		axis.text.y=element_text(colour="black"),
		axis.ticks=element_line(colour="black"),
		plot.margin=unit(c(1,1,1,1),"line")
	)
dev.off()


vsd <- vst(dds, blind=FALSE)
png(paste(prefix,".all.samples.heatmap.vsd.png",sep=""),width=900,height=900)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE)
dev.off()

