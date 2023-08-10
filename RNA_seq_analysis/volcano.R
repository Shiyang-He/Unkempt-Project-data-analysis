argv <- commandArgs(TRUE)
library(ggplot2)
library(grid)
library(stringr)
library(gridExtra)
library(ggforce)
library(dplyr)

deg=read.table(argv[1],header=T,row.names=1)
filename=str_replace(argv[1],".xls","")
b=strsplit(argv[1],"\\.vs\\.")[[1]][1]
c=strsplit(argv[1],"\\.vs\\.")[[1]][2]
b=strsplit(b,"\\.")[[1]][-1]
c=strsplit(c,"\\.")[[1]][1]

conditions=c(b,c)
df <- read.table(text = "",colClasses = c("character", "character", "character", "integer"),col.names=c("Condition","Group","Bound","Count"))
for (i in c(1:8)){
	a=read.table(paste(argv[1],".g",i,".bound",sep=""),header=T,sep="\t",row.names=1)
	for (j in conditions){
		k=j
		if (j=="3M"){
                	j="X3M_bound"
        	}else{
                	j=paste(j,"_bound",sep="")
        	}
		if (j %in% names(a)){
			tmp_df=data.frame(table(a[names(a)==j]))
			tmp_df=cbind(rep(k,length(tmp_df[,1])),rep(paste("g",i,sep=""),length(tmp_df[,1])),tmp_df)
			colnames(tmp_df)=c("Condition","Group","Bound","Count")
			df=rbind(df,tmp_df)
		}
	}
}

df_pie= left_join(df, 
		df %>% 
		group_by(Group,Condition) %>% 
		summarize( Count_total=sum(Count))) %>% 
		group_by(Group) %>%
		mutate(end_angle = 2*pi*cumsum(Count)/Count_total,      # ending angle for each pie slice
         	start_angle = lag(end_angle, default = 0),   		# starting angle for each pie slice
         	mid_angle = 0.5*(start_angle + end_angle))   		# middle of each pie slice, for the text label

df_pie$Bound<-factor(df_pie$Bound, levels=c("Highly bound","Moderate bound","Not/Lowly bound"))
df_pie$Percentage=df_pie$Count/df_pie$Count_total*100
df_pie$Group=paste(df_pie$Group,"(",df_pie$Count_total,")",sep="")

bar_chart=ggplot(data=df_pie, aes(x=" ", y=Percentage,group=Bound,colour=Bound,fill=Bound))+
	geom_bar(width = 1, stat = "identity") + 
	facet_grid(~ Group)+
	xlab("")+
	ylab("")+
	scale_fill_manual(values=c(`Moderate bound`="burlywood2",`Highly bound`="brown3",`Not/Lowly bound`="beige"))+
	scale_color_manual(values=c(`Moderate bound`="burlywood2",`Highly bound`="brown3",`Not/Lowly bound`="beige"))+
	scale_y_continuous(breaks=seq(0,100,by=50))+
	theme(axis.ticks = element_blank(), 
		strip.text.x = element_text(size = 12),
		legend.title = element_text(size=14),
		legend.text = element_text(size=10),
		text = element_text(size=rel(5.5)),
                axis.text.y = element_text(size = rel(6.2)),
		#axis.text.x = element_blank(),
		legend.key = element_blank(), 
		strip.background = element_rect(colour="white", fill="white"),
		axis.line = element_line(colour = "white"),
		panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())
legend=as.numeric(argv[2])
if (legend==0){bar_chart=bar_chart+theme(legend.position="none")}


xlab=paste("log2FC(",b,".vs.",c,")",sep="")
x=deg[,8]
y=-log10(deg$padj)
y[y>5]<-5 #bring all high value to 5
x[x< -3] <- -3
x[x> 3] <- 3
dat=data.frame(x,y)
colnames(dat)=c("x","y")

significant_dn_normal=filter(dat, x <= -1  & ( y >= -log10(0.01) & y<5)) #  
significant_up_normal=filter(dat, x >= 1  & ( y >= -log10(0.01) & y<5))

nochange=filter(dat,y < -log10(0.01) | (x>-1 & x<1)) #all gray dots. Under the p value cutoff

high_sig_dn=filter(dat, x <= -1 & y >=5)
high_sig_up=filter(dat, x >= 1  & y >=5)
high_nosig=filter(dat, (x > -1 & x < 1 ) & y>=5)
extreme_sig_dn=filter(dat, x <= -3 & y > -log10(0.01))
extreme_sig_up=filter(dat, x >= 3 & y > -log10(0.01))
extreme_nosig=filter(dat, (x <= -3 | x >= 3) & y < -log10(0.01))

volcano_plot = ggplot(nochange,aes(x,y))+geom_point(col=rgb(0.5,0.5,0.5,0.3))
if (length(significant_up_normal$x)>0){volcano_plot=volcano_plot+geom_point(data=significant_up_normal,mapping=aes(x=x,y=y),col=rgb(0.973, 0.463, 0.427,0.3))}
if (length(significant_dn_normal$x)>0){volcano_plot=volcano_plot+geom_point(data=significant_dn_normal,mapping=aes(x=x,y=y),col=rgb(0.486, 0.682, 0, 0.3))}
if (length(extreme_nosig$x)>0){volcano_plot=volcano_plot+geom_point(data=extreme_nosig,mapping=aes(x=x,y=y),shape=25,colour=rgb(0.5,0.5,0.5,0.3),fill=rgb(0.5,0.5,0.5,0.3))} 
if (length(extreme_sig_up$x)>0){volcano_plot=volcano_plot+geom_point(data=extreme_sig_up,mapping=aes(x=x,y=y),shape=25,colour=rgb(0.973, 0.463, 0.427,0.8),fill=rgb(0.973, 0.463, 0.427,0.3))}
if (length(extreme_sig_dn$x)>0){volcano_plot=volcano_plot+geom_point(data=extreme_sig_dn,mapping=aes(x=x,y=y),shape=25,colour=rgb(0.486, 0.682, 0,0.8),fill=rgb(0.486, 0.682, 0,0.3))}
if (length(high_nosig$x)>0){volcano_plot=volcano_plot+geom_point(data=high_nosig,mapping=aes(x=x,y=y),shape=25,colour=rgb(0.5,0.5,0.5,0.3),fill=rgb(0.5,0.5,0.5,0.3))} 
if (length(high_sig_up$x)>0){volcano_plot=volcano_plot+geom_point(data=high_sig_up,mapping=aes(x=x,y=y),shape=25,colour=rgb(0.973, 0.463, 0.427,0.8),fill=rgb(0.973, 0.463, 0.427,0.3))}
if (length(high_sig_dn$x)>0){volcano_plot=volcano_plot+geom_point(data=high_sig_dn,mapping=aes(x=x,y=y),shape=25,colour=rgb(0.486, 0.682, 0,0.8),fill=rgb(0.486, 0.682, 0,0.3))}
volcano_plot=volcano_plot+
	xlab(xlab)+ylab("-log10(padj)")+
	geom_segment(aes(x = 1, xend=1, y= -log10(0.01), yend=5.5), lty=4)+
        geom_segment(aes(x = 2, xend=2, y= -log10(0.01), yend=5.5), lty=4)+
        geom_segment(aes(x = 0.5849625, xend=0.5849625, y= -log10(0.01), yend=5.5), lty=4)+
        geom_segment(aes(x = 0, xend=0, y= -log10(0.01), yend=5.5), lty=4)+
        geom_segment(aes(x = -0.5849625, xend=-0.5849625, y= -log10(0.01), yend=5.5), lty=4)+
        geom_segment(aes(x = -2, xend=-2, y= -log10(0.01), yend=5.5), lty=4)+
        geom_segment(aes(x = -1, xend=-1, y= -log10(0.01), yend=5.5), lty=4)+
        geom_hline(yintercept=-log10(0.01),lty=4)+
        annotate("text",label="g1",x=-2.5,y=5.5,size=7)+
        annotate("text",label="g2",x=-1.5,y=5.5,size=7)+
        annotate("text",label="g3",x=-0.77,y=5.5,size=7)+
        annotate("text",label="g4",x=-0.29,y=5.5,size=7)+
        annotate("text",label="g5",x=0.29,y=5.5,size=7)+
        annotate("text",label="g6",x=0.77,y=5.5,size=7)+
        annotate("text",label="g7",x=1.5,y=5.5,size=7)+
        annotate("text",label="g8",x=2.5,y=5.5,size=7)+
	scale_x_continuous(breaks=seq(-3,3,by=1))+
	scale_y_continuous(breaks=seq(0,5,by=1))+
	theme(axis.line = element_line(colour = "black"),
    		panel.grid.major = element_blank(),
   		panel.grid.minor = element_blank(),
		text = element_text(size=rel(5.5)),
                axis.text.x = element_text(size = rel(6.2)),
                axis.text.y = element_text(size = rel(6.2)),
		panel.border = element_blank(),
	 	panel.background = element_blank())

width=8
if (legend==1){width=10}
png(paste(filename,".png",sep="."),width=width,height=7,res=350,units="in")
if (legend==1){
grid.arrange(bar_chart,arrangeGrob(volcano_plot,p,ncol=2,widths=c(21.5,3.9)),heights=c(0.6,1),top = textGrob(paste(b,".vs.",c,sep=""),gp=gpar(fontsize=20)))
}else{
grid.arrange(bar_chart,volcano_plot,ncol=1,heights=c(0.6,1),top = textGrob(paste(b,".vs.",c,sep=""),gp=gpar(fontsize=20)))
}
dev.off()

