library(ComplexHeatmap)
library(circlize)
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
draw_line_plot=function(cc,ylim,ylab=""){
	p=ggplot(cc,aes(x=pos,y=U_rich))+
		geom_smooth(span = 0.2,se=F)+
		geom_smooth(data=cc,aes(x=pos,y=UAG),span=0.2,se=F,col="red")+
		geom_smooth(data=cc,aes(x=pos,y=Other),span=0.2,se=F,col="black")+
		ylim(ylim)+
		theme(
			axis.line=element_line(colour = "black"),
			axis.text.x=element_blank(),
			axis.title.x=element_blank(),
			legend.position="none",
           		panel.background=element_blank(),
			panel.border=element_blank(),
			panel.grid.major=element_blank(),
           		panel.grid.minor=element_blank(),
			plot.background=element_blank())
		if (ylab !=""){
			p=p+ylab(ylab)
		}
	return(p)
}

draw_mut=function(df,cc,title,ylim){
	ha_top = HeatmapAnnotation(" " = anno_lines(cc, gp = gpar(col = c("orange","brown4","gray"),lwd=3),height = unit(4, "cm"),ylim = c(0, ylim)),annotation_name_side="left")
	split=c(1,rep(2,length(rownames(df))-2),3)
	col_fun= colorRamp2(c(0,2,4,6,8), c("#FFFFCC", "#FECB76","#FD9040", "#FC4E2A","#A91A27"))
	p=Heatmap(df,cluster_rows=F,cluster_columns=F,column_names_rot=0,name="Enrichment",col=col_fun,row_split = split,top_annotation = ha_top,row_names_side="left",column_title = title,column_title_side="top")
	return(p)
}
draw_wt=function(df,cc,ylim){
	ha_top = HeatmapAnnotation(Enrichment = anno_lines(cc, gp = gpar(col = c("orange","brown4","gray"),lwd=3),height = unit(4, "cm"),ylim = c(0, ylim)),annotation_name_side="left")
	split=c(1,rep(2,length(rownames(df))-2),3)
	col_fun= colorRamp2(c(0,2,4,6,8), c("#FFFFCC", "#FECB76","#FD9040", "#FC4E2A","#A91A27"))
	#ha_left=rowAnnotation(empty = anno_empty(border = FALSE)," "=anno_block(gp = gpar(fill = c("red","orange","gray")),labels=c("","","")))
	#ha_left=rowAnnotation(" "=anno_block(gp = gpar(col=c("red","orange","gray"),fill = c("red","orange","gray")),labels=c("","","")))
	ha_left=rowAnnotation(" "=anno_block(gp = gpar(col=c("brown4","orange","gray"),size=0.5,fill = c("brown4","orange","gray")),labels=c("","","")))
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
wt_p=draw_wt(wt_output[[1]],cc_wt,ylim)
dpam2_p=draw_mut(dpam2_output[[1]],cc_dpam2,"dPAM2",ylim)
m3_p=draw_mut(m3_output[[1]],cc_3m,"3M",ylim)


lgd=Legend(labels=c("U/A rich","UAG","Other"),title="Triplet",type = "lines",legend_gp = gpar(col = c("orange","brown4","gray"),lty=1,lwd=3),background ="white")
pdf(argv[7],width=12,height=6)
draw(wt_p+m3_p+dpam2_p,column_title = "Distance to Peak (nt)",column_title_side="bottom")
pushViewport(viewport(width = 0.78, height = 0.77))
draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right"))
dev.off()
