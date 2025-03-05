library(ComplexHeatmap)
library(circlize)
library(ggplot2)
argv=commandArgs(TRUE)
wt_full=read.table(argv[1],row.names=1)
wt_short=read.table(argv[2],row.names=1)
m3_full=read.table(argv[3],row.names=1)
m3_short=read.table(argv[4],row.names=1)
dpam2_full=read.table(argv[5],row.names=1)
dpam2_short=read.table(argv[6],row.names=1)

`%!in%` <- Negate(`%in%`)
urich=c("AUU","CUU","GUU","UUA","UUU","UUC","UAC","UAU","UCG","UCA","UCC","UAA","UUG")
process_wt=function(a,b){
        uag_df_long=subset(a,rownames(a) == "UAG")
        uag_df_short=subset(b,rownames(b) == "UAG")
        uag_median=median(as.numeric(uag_df_long))
        uag_df_short_enriched=uag_df_short/uag_median
	non_uag_df_long=subset(a,rownames(a) != "UAG")
	non_uag_short=subset(b,rownames(b) != "UAG")
	non_uag_df_long$max <- apply(non_uag_df_long, 1, max, na.rm=TRUE)
	non_uag_df_long=non_uag_df_long[order(-non_uag_df_long$max),]
	urich_df_long=head(non_uag_df_long,n=10)[1:31]
	add_on=c("UAA", "UAU", "AUA")
	urich_df_long=subset(urich_df_long,rownames(urich_df_long) %!in% add_on)
	add_on_df=subset(non_uag_df_long,rownames(non_uag_df_long) %in% add_on)[1:31]
	urich_df_long=rbind(urich_df_long,add_on_df)
	urich_df_short=subset(b,rownames(b) %in% rownames(urich_df_long))
	urich_df_short=urich_df_short[rownames(urich_df_long),]
	other_df_long=subset(non_uag_df_long,rownames(non_uag_df_long) %!in% rownames(urich_df_long)) 
	other_df_short=subset(non_uag_short,rownames(non_uag_short) %!in% rownames(urich_df_long))
	other_df_long_avg=colSums(other_df_long)/length(rownames(other_df_long))
	urich_df_long_avg=colSums(urich_df_long)/length(rownames(urich_df_long))
        urich_median=median(urich_df_long_avg)
        other_median=median(other_df_long_avg)
	other_df_short_enriched=(colSums(other_df_short)/length(rownames(other_df_short)))/other_median
        other_df_short_enriched=data.frame(t(other_df_short_enriched))
        rownames(other_df_short_enriched)="Other"
	urich_df_short_enriched=urich_df_short/urich_median
	urich_df_short_enriched_avg=colSums(urich_df_short_enriched)/length(rownames(urich_df_long))
        c=rownames(urich_df_short_enriched) ## a vector tells which group to show
	df=rbind(uag_df_short_enriched,urich_df_short_enriched,other_df_short_enriched)
	colnames(df)=c(rep("",5),-10,rep("",4),-5,rep("",4),0,rep("",4),5,rep("",4),10,rep("",5))
	cc=cbind(as.numeric(urich_df_short_enriched_avg),as.numeric(other_df_short_enriched),as.numeric(uag_df_short_enriched),c(-15:15))
        colnames(cc)=c("U_rich","Other","UAG","pos")
        cc=data.frame(cc)
	output=list(df,c,cc)
	# df:	a dataframe for heatmap
	# c:	a sorted column name for mutation to follow
	# cc:	a matrix for drawing lines
}
process_mut=function(a,b,c){
        uag_df_long=subset(a,rownames(a) == "UAG")
        uag_df_short=subset(b,rownames(b) == "UAG")
        uag_median=median(as.numeric(uag_df_long))
        uag_df_short_enriched=uag_df_short/uag_median

	urich_df_long=subset(a,rownames(a) %in% c)
	urich_df_long=urich_df_long[c,]
	urich_df_short=subset(b,rownames(b) %in% c)
	urich_df_short=urich_df_short[c,]
	other_df_long=subset(a,rownames(a) %!in% rownames(urich_df_long) & rownames(a) != "UAG") 
	other_df_short=subset(b,rownames(b) %!in% rownames(urich_df_long) & rownames(b) != "UAG")

	other_df_long_avg=colSums(other_df_long)/length(rownames(other_df_long))
	urich_df_long_avg=colSums(urich_df_long)/length(rownames(urich_df_long))
        urich_median=median(urich_df_long_avg)
        other_median=median(other_df_long_avg)
	other_df_short_enriched=(colSums(other_df_short)/length(rownames(other_df_short)))/other_median
        other_df_short_enriched=data.frame(t(other_df_short_enriched))
        rownames(other_df_short_enriched)="Other"
	urich_df_short_enriched=urich_df_short/urich_median
	urich_df_short_enriched_avg=colSums(urich_df_short_enriched)/length(c)
	df=rbind(uag_df_short_enriched,urich_df_short_enriched,other_df_short_enriched)
	colnames(df)=c(rep("",5),-10,rep("",4),-5,rep("",4),0,rep("",4),5,rep("",4),10,rep("",5))
	cc=cbind(as.numeric(urich_df_short_enriched_avg),as.numeric(other_df_short_enriched),as.numeric(uag_df_short_enriched),c(-15:15))
        colnames(cc)=c("U_rich","Other","UAG","pos")
        cc=data.frame(cc)
	output=list(df,cc)
	# df:	a dataframe for heatmap
	# c:	a sorted column name for mutation to follow
	# cc:	a matrix for drawing lines
}
loess_df=function(cc){
	U_rich_lo=loess(U_rich ~ pos,cc,span=0.2)
	UAG_lo=loess(UAG ~ pos,cc,span=0.2)
	Other_lo=loess(Other ~ pos, cc, span=0.2)
	cc=data.frame(U_rich=predict(U_rich_lo),UAG=predict(UAG_lo),Other=predict(Other_lo))
	#cc$U_rich=predict(U_rich_lo)
	#cc$UAG=predict(UAG_lo)
	#cc$Other=predict(Other_lo)
	return(cc)
}
draw_line_plot=function(cc,ylim=c(0,4)){
	p=ggplot(cc,aes(x=Pos,y=WT))+
		geom_smooth(span = 0.2,se=F,col="red")+
		geom_smooth(data=cc,aes(x=Pos,y=M3),span=0.2,se=F,col="green")+
		geom_smooth(data=cc,aes(x=Pos,y=dPAM2),span=0.2,se=F,col="blue")+
		scale_y_continuous(limits=ylim,expand=c(0,0))+
		ylab("Enrichment")+
		xlab("Distance to Peak (nt)")+
		facet_grid(~condition)+
		theme(
			axis.ticks.length=unit(.25, "cm"),
        	        axis.ticks= element_line (linewidth=1.2,colour = "black"),
	                axis.line = element_line(colour = "black",linewidth=1.2),
			text = element_text(size=30),
			axis.text.x = element_text(size = 30,colour="black"),
			axis.text.y = element_text(size = 30,colour="black"),
			plot.title = element_text (hjust = 0.5),
           		panel.background=element_blank(),
			panel.border=element_blank(),
			panel.grid.major=element_blank(),
           		panel.grid.minor=element_blank(),
			plot.background=element_blank())
	return(p)
}

draw_mut=function(df,cc,title,ylim){
	ha_top = HeatmapAnnotation(" " = anno_lines(cc, gp = gpar(col = c("#4488C7","red","gray"),lwd=3),height = unit(4, "cm"),ylim = c(0, ylim)),annotation_name_side="left")
	split=c(1,rep(2,length(rownames(df))-2),3)
	col_fun= colorRamp2(c(0,2,4,6,8), c("#FFFFCC", "#FECB76","#FD9040", "#FC4E2A","#A91A27"))
	p=Heatmap(df,cluster_rows=F,cluster_columns=F,column_names_rot=0,name="Enrichment",col=col_fun,row_split = split,top_annotation = ha_top,row_names_side="left",column_title = title,column_title_side="top")
	return(p)
}
draw_wt=function(df,cc,ylim){
	ha_top = HeatmapAnnotation(Enrichment = anno_lines(cc, gp = gpar(col = c("#4488C7","red","gray"),lwd=3),height = unit(4, "cm"),ylim = c(0, ylim)),annotation_name_side="left")
	split=c(1,rep(2,length(rownames(df))-2),3)
	col_fun= colorRamp2(c(0,2,4,6,8), c("#FFFFCC", "#FECB76","#FD9040", "#FC4E2A","#A91A27"))
	#ha_left=rowAnnotation(empty = anno_empty(border = FALSE)," "=anno_block(gp = gpar(fill = c("red","#4488C7","gray")),labels=c("","","")))
	#ha_left=rowAnnotation(" "=anno_block(gp = gpar(col=c("red","#4488C7","gray"),fill = c("red","#4488C7","gray")),labels=c("","","")))
	ha_left=rowAnnotation(" "=anno_block(gp = gpar(col=c("red","#4488C7","gray"),size=0.5,fill = c("red","#4488C7","gray")),labels=c("","","")))
	p=Heatmap(df,cluster_rows=F,cluster_columns=F,column_names_rot=0,name="Enrichment",col=col_fun,row_split = split,left_annotation=ha_left,top_annotation = ha_top,row_names_side="left",column_title = "WT",column_title_side="top")
	return(p)
}
wt_output=process_wt(wt_full,wt_short)
cc_wt=loess_df(wt_output[[3]])
m3_output=process_mut(m3_full,m3_short,wt_output[[2]])
cc_3m=loess_df(m3_output[[2]])
dpam2_output=process_mut(dpam2_full,dpam2_short,wt_output[[2]])
cc_dpam2=loess_df(dpam2_output[[2]])
ylim=as.numeric(argv[8])

U_rich_df=data.frame(WT=cc_wt$U_rich,M3=cc_3m$U_rich,dPAM2=cc_dpam2$U_rich,Pos=c(-15:15))
UAG_df=data.frame(WT=cc_wt$UAG,M3=cc_3m$UAG,dPAM2=cc_dpam2$UAG,Pos=c(-15:15))
Other_df=data.frame(WT=cc_wt$Other,M3=cc_3m$Other,dPAM2=cc_dpam2$Other,Pos=c(-15:15))
U_rich_df$condition="U/A rich"
UAG_df$condition="UAG"
Other_df$condition="Other"
dffff=rbind(UAG_df,U_rich_df,Other_df)
dffff$condition=factor(dffff$condition,levels=unique(dffff$condition))

pdf(argv[7],width=12,height=6)
draw_line_plot(dffff,ylim=c(0,ylim))
lgd=Legend(labels=c("WT","3M","dPAM2"),title="Condition",type = "lines",legend_gp = gpar(col = c("red","green","blue"),lty=1,lwd=3,fontsize=20),background ="white")
pushViewport(viewport(width = 0.85, height = 0.68))
draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right"))
dev.off()
