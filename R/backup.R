
library(ggpubr)
library(pheatmap)
library(ape)
library(ggplot2)
library(phyloseq)
library(matrixStats)
library(vegan)

 #' Function to parse the content of OTU and KO files
 #'
 #' @param otutable_path The file path of otu
 #' @param otutaxa_path The file path of rdp classification
 #' @return OTU oject including otu tab
 #' @examples
 #' readingOtuData("C:/Microbiome/MMDRPipeline/data/ben1/ben1_16s_otu_table.tsv","C:/Microbiome/MMDRPipeline/data/ben1/ben1_16s_otu_taxa.tsv","C:/Microbiome/MMDRPipeline/data/ben1/ben1_16s_sp_metadata.tsv")
 
 readingOtuData = function(otutable_path, otutaxa_path,spmeta_path){
   # Read the actual count data
   otu<- NULL
   otu$otu.tab <- read.csv(otutable_path, header=TRUE,sep="\t",check.names=FALSE,row.names= 1)
   otu$otu.taxa <- read.csv(otutaxa_path, header=TRUE,sep="\t",check.names=FALSE,row.names= 1)
   otu$sp.meta <- read.csv(spmeta_path, header=TRUE,sep="\t",check.names=FALSE,row.names= 1)
	rownames(otu$otu.tab)=gsub("\\\\", "", rownames(otu$otu.tab))
	rownames(otu$otu.tab)=gsub(" ", "_", rownames(otu$otu.tab))
	rownames(otu$otu.tab)=gsub("\"_", "", rownames(otu$otu.tab))
	rownames(otu$otu.tab)=gsub("\"", "", rownames(otu$otu.tab))
   return(otu)
 }
 
 #' Function to parse the content of KO files
 #'
 #' @param kotable_path The file path of otu
 #' @param kometa_path The file path of ko description
 #' @return ko ojbect including ko.tab and ko.meta
 #' @examples
 #' readingFunctionData("C:/Microbiome/MMDRPipeline/data/ben1/ben1_16s_otu_table.tsv","C:/Microbiome/MMDRPipeline/data/ben1/ben1_16s_otu_taxa.tsv")
 
 readingFunctionData = function(kotable_path, kometa_path){
 
   # Read the actual count data
 
   ko <-NULL
   ko$ko.tab <- read.csv(kotable_path, header=TRUE,sep="\t",check.names=FALSE,row.names = 1)
   ko$ko.meta <- read.csv(kometa_path, header=TRUE,sep="\t",check.names=FALSE,row.names = 1)
   return(ko)
 }
 
 
 
 
 ################ PART B: PLOT : PCA, HEATMAP, BOXPLOT, ETETETEET  ###########################
 
 
 #' Function to generate PCA plot of abundance table : OTU,GENE,KO
 #'
 #' @param abund.tab otu/gene table as matrix
 #' @param sp.meta sample metadata as.matrix (ID,NICKNAME,CONDITION,SUBJECT,CATEGORY)
 #' @param groups list of selected biological condition for plot
 #' @param meanCountCutOff The mean count value a gene/otu must have to be considered for the plot.
 #' @return PCA plot
 #' @examples
 #' AbundPcaPlot(OTUDATA$otu.tab,OTUDATA$sp.meta,groups=c('ain','protein','cholic'))
 #' AbundPcaPlot(OTUDATA$otu.tab,OTUDATA$sp.meta,groups=c('protein','cholic'))
 #' AbundPcaPlot(KODATA$ko.tab,OTUDATA$sp.meta,groups=c('protein','cholic'))
 #' AbundPcaPlot(OTUDATA$otu.tab,OTUDATA$sp.meta,meanCountCutOff=1,groups='ALL',title="JBCC: OTU Abundance ",label="TREATMENT",color="GENDER",size_group=3,size_label=2)
 #' AbundPcaPlot(OTUDATA$otu.tab,OTUDATA$sp.meta,meanCountCutOff=1,groups='ALL',title="JBCC: OTU Abundance ")
 #' AbundPcaPlot(OTUDATA$otu.tab,OTUDATA$sp.meta,meanCountCutOff=1,groups='ALL',title="JBCC: OTU Abundance ",label="TREATMENT",size_group=3)
 #' AbundPcaPlot(OTUDATA$otu.tab,OTUDATA$sp.meta,meanCountCutOff=1,groups='ALL',title="JBCD: PCA by OTU Abundance ",color="TREATMENT",label='NICKNAME',size_group=3,size_label=2)


 AbundPcaPlot = function(abund.tab,sp.meta, meanCountCutOff=1,groups="ALL",title = "Title",label=NULL,color=NULL,size_group=2,size_label=2, distance="bray", prop=FALSE){
   
   listall1=rownames(sp.meta)
   listall2=colnames(abund.tab)
   listall=intersect(listall1,listall2)
 
   abund.tab=abund.tab[,listall]
   otu.mean <- apply(X = abund.tab, MARGIN = 1, FUN = mean)
   otu.filter <- abund.tab[otu.mean >= meanCountCutOff,]
   otu=t(otu.filter)
   sp.meta=sp.meta[listall,]

   if (prop==TRUE) {
	otu=prop.table(otu)
	print("proportion")
   }

   dist=vegdist(otu,method=distance)


   if (groups=="ALL") 
   {
  	condition=sp.meta$CONDITION
   	pca=prcomp(dist)
	eigs <- pca$sdev^2
	pc_1_axis=round((eigs[1] / sum(eigs))*100)
	pc_2_axis=round((eigs[2] / sum(eigs))*100)
   	data_df <- data.frame(pc_1 = pca$x[,1], pc_2 = pca$x[,2], shape = condition)
   }
   else
   {
	temp=sp.meta[sp.meta$CONDITION %in% groups,]
 	mylist=rownames(temp)
 	mycondition=temp$CONDITION
 	otu=otu[mylist,]
 	pca=prcomp(dist)
	eigs <- pca$sdev^2
	eigs[1] / sum(eigs)
	pc_1_axis=round((eigs[1] / sum(eigs))*100)
	pc_2_axis=round((eigs[2] / sum(eigs))*100)

	data_df <- data.frame(pc_1 = pca$x[,1], pc_2 = pca$x[,2], shape = mycondition)
 
   }

   if (is.null(label)) {
   	if (is.null(color)) {
		p <-ggplot(data_df, aes(pc_1,pc_2)) + geom_point(aes(shape=shape),size=size_group) + labs(title = title,x=paste0("pc_1:",pc_1_axis,"%"),y=paste0("pc_2:",pc_2_axis,"%"))
	}
	else {
		color_=sp.meta[listall,][,paste0(color)]
		p <-ggplot(data_df, aes(pc_1,pc_2)) + geom_point(aes(color=color_,shape=shape),size=size_group) + labs(title = title,x=paste0("pc_1:",pc_1_axis,"%"),y=paste0("pc_2:",pc_2_axis,"%"))
	}
   }
   else {
	meta_label=sp.meta[listall,][,paste0(label)]
	
	if (is.null(color)) {
   		p <-ggplot(data_df, aes(pc_1,pc_2)) + geom_point(aes(color=shape,shape=shape),size=size_group) + labs(title = title,x=paste0("pc_1:",pc_1_axis,"%"),y=paste0("pc_2:",pc_2_axis,"%")) +geom_text(mapping = aes(label=meta_label),size=size_label,hjust = -0.8)
	}	
	else {
		color_=sp.meta[listall,][,paste0(color)]
   		p <-ggplot(data_df, aes(pc_1,pc_2)) + geom_point(aes(color=color_,shape=shape),size=size_group) + labs(title = title,x=paste0("pc_1:",pc_1_axis,"%"),y=paste0("pc_2:",pc_2_axis,"%")) +geom_text(mapping = aes(label=meta_label),size=size_label,hjust = -0.8)
	}
   }
   p
 }
 
  AbundPcaPlotKO = function(abund.tab,sp.meta, meanCountCutOff=1,groups="ALL",title = "Title",label=NULL,color=NULL,size_group=2,size_label=2){
   
   listall1=rownames(sp.meta)
   listall2=colnames(abund.tab)
   listall=intersect(listall1,listall2)
 
   abund.tab=abund.tab[,listall]
   otu.mean <- apply(X = abund.tab, MARGIN = 1, FUN = mean)
   otu.filter <- abund.tab[otu.mean >= meanCountCutOff,]
   otu=t(otu.filter)
   sp.meta=sp.meta[listall,]

 

   if (groups=="ALL") 
   {
  	condition=sp.meta$CONDITION
   	pca=prcomp(log2(otu+1))

	eigs <- pca$sdev^2
	pc_1_axis=round((eigs[1] / sum(eigs))*100)
	pc_2_axis=round((eigs[2] / sum(eigs))*100)
   	data_df <- data.frame(pc_1 = pca$x[,1], pc_2 = pca$x[,2], shape = condition)
   }
   else
   {
	temp=sp.meta[sp.meta$CONDITION %in% groups,]
 	mylist=rownames(temp)
 	mycondition=temp$CONDITION

 	otu=otu[mylist,]
 	pca=prcomp(log2(otu+1))

	eigs <- pca$sdev^2
	eigs[1] / sum(eigs)
	pc_1_axis=round((eigs[1] / sum(eigs))*100)
	pc_2_axis=round((eigs[2] / sum(eigs))*100)

	data_df <- data.frame(pc_1 = pca$x[,1], pc_2 = pca$x[,2], shape = mycondition)
 
   }

   if (is.null(label)) {
   	if (is.null(color)) {
		p <-ggplot(data_df, aes(pc_1,pc_2)) + geom_point(aes(shape=shape),size=size_group) + labs(title = title,x=paste0("pc_1:",pc_1_axis,"%"),y=paste0("pc_2:",pc_2_axis,"%"))
	}
	else {
		color_=sp.meta[listall,][,paste0(color)]
		p <-ggplot(data_df, aes(pc_1,pc_2)) + geom_point(aes(color=color_,shape=shape),size=size_group) + labs(title = title,x=paste0("pc_1:",pc_1_axis,"%"),y=paste0("pc_2:",pc_2_axis,"%"))
	}
   }
   else {
	meta_label=sp.meta[listall,][,paste0(label)]
	
	if (is.null(color)) {
   		p <-ggplot(data_df, aes(pc_1,pc_2)) + geom_point(aes(color=shape,shape=shape),size=size_group) + labs(title = title,x=paste0("pc_1:",pc_1_axis,"%"),y=paste0("pc_2:",pc_2_axis,"%")) +geom_text(mapping = aes(label=meta_label),size=size_label,hjust = -0.8)
	}	
	else {
		color_=sp.meta[listall,][,paste0(color)]
   		p <-ggplot(data_df, aes(pc_1,pc_2)) + geom_point(aes(color=color_,shape=shape),size=size_group) + labs(title = title,x=paste0("pc_1:",pc_1_axis,"%"),y=paste0("pc_2:",pc_2_axis,"%")) +geom_text(mapping = aes(label=meta_label),size=size_label,hjust = -0.8)
	}
   }
   p
 }
 
 
 #' Function to create a heatmap  of abundance data (ROWs:OTU,GENE,KO, COLS: SAMPLE NAMES)  
 #'
 #' @param abund.tab abundance table as.matrix
 #' @param sp.meta sample metadata as.matrix (ID,NICKNAME,CONDITION,SUBJECT,CATEGORY)
 #' @param groups list of biological condition(Healthy vs Disease)
 #' @return Heatmap plot
 #' @examples
 #' AbundPlotHeatmap(OTUDATA$otu.tab,OTUDATA$sp.meta, meanCountCutOff = 50, groups=c('protein','cholic'),fontsize_row=5,fontsize_col=5)
 
 AbundPlotHeatmap <- function(abund.tab,sp.meta, meanCountCutOff = 1, groups="ALL",fontsize_row=5,fontsize_col=5, id_list=NULL) {
 
 	# merge samle names from abundance table and metadata
	listall1=rownames(sp.meta)
	listall2=colnames(abund.tab)
	listall=intersect(listall1,listall2)
	counts=abund.tab[,listall]

	# Remove counts with too low count
	counts.mean <- apply(X=counts, MARGIN = 1, FUN = mean)
	expressed <- counts.mean >= meanCountCutOff
	counts <- counts[expressed,]
 

	# if gene/otu list (id_list) is not null, then, create heatmap with only id_list
	if (!is.null(id_list)) {
		k1=rownames(abund.tab)
		k2=rownames(id_list)
		k=intersect(k1,k2)
		counts=counts[k,]
		id_list=id_list[k,]
		rownames(counts)=paste(rownames(counts),"_",id_list$GBM_NAME)
	}


	log2.counts <- log2(counts + 1)
 
	if (groups=="ALL") {
		annot <- data.frame(Group = sp.meta$CONDITION,Subject= sp.meta$SUBJECT)
		rownames(annot) <- rownames(sp.meta)
	} else
	{
		temp=sp.meta[sp.meta$CONDITION %in% groups,]
		mylist=rownames(temp)
		annot <- data.frame(Group =temp$CONDITION,Subject= sp.meta$SUBJECT)
		rownames(annot) <- rownames(temp)
		log2.counts=log2.counts[,mylist]
 
	}
	pheatmap(log2.counts, scale="row", show_rownames = F, show_colnames = T, annotation_col = annot,fontsize_row=fontsize_row,fontsize_col=fontsize_col)
 
 }
 
 #  VolcanoPlot(diff,fdr=0.05,fc=4,title="BEN1 KO Differential Abundance: cholic" %->% "protein ")
 
 
 VolcanoPlot <-function(diffdata,id_list=NULL, fdr=0.05, fc=2,title="Differential Abundance ",top=10,label_size=5) {
 	
	#diffdata$log2FoldChange=as.numeric(diffdata$log2FoldChange)
	#diffdata$pvalue=as.numeric(diffdata$pvalue)
	#diffdata$padj=as.numeric(diffdata$padj)

	NS_FC=paste0("FC = ",fc)
	NS_FDR=paste0("P_adj = ",fdr)
        title=paste0(title," ; ",NS_FC," ", NS_FDR, " top=",top)
 
 	ggmaplot(diffdata,
	id_list,
	main = title,
	fdr = fdr,
	fc =fc,
	size = 1,
	palette = c("#B31B21", "#1465AC", "darkgray"),
	genenames = as.vector(rownames(diffdata)),
	legend = "top", 
	top = top,
	font.label = c("bold", label_size),
	font.legend = "bold",
	font.main = "bold",
	xlab="Log mean abundance of two groups",
	ylab="Log2 fold-change of mean abundance ",
	col=rownames(diffdata),
	ggtheme = ggplot2::theme_minimal()
	)
 
 }

 VolcanoPlotCV <-function(diffdata, fdr=0.05, fc=2,title="Differential Abundance ",top=40,label_size=5) {
 	
	diffdata$fc=diffdata$fc_cv
	diffdata$log2FoldChange=as.numeric(diffdata$log2FoldChange_cv)

	diffdata$pvalue=as.numeric(diffdata$pvalue)
	diffdata$padj=as.numeric(diffdata$padj)

	NS_FC=paste0("FC = ",fc)
	NS_FDR=paste0("P_adj(for NS) = ",fdr)
        title=paste0(title," ; ",NS_FC," ", NS_FDR)
 
 	ggmaplot(diffdata, 
	main = title,
	fdr = fdr,
	fc =fc,
	size = 1,
	palette = c("#B31B21", "#1465AC", "darkgray"),
	genenames = as.vector(rownames(diffdata)),
	legend = "top", 
	top = top,
	font.label = c("bold", label_size),
	font.legend = "bold",
	font.main = "bold",
	xlab="Log mean abundance of two groups",
	ylab="Log2 Fold-Change of Coefficient of Variation ",
	col=rownames(diffdata),
	ggtheme = ggplot2::theme_minimal()
	)
 
 }
  
 
 ################ PART C: SUMMARY REPORTING : CSV FILE   ###########################
 
 
 
 #' Function for get richness and count summary for each sampleid using phyloseq  
 #'
 #' @param abund.tab Integer abundance table including only interger values (row names are OTU/KO/GENES)
 #' @return richness summary table
 #' @examples
 #' otu_richness=richness_phyloseq(OTUDATA$otu.tab)
 #' ko_richness=richness_phyloseq(round(KODATA$ko.tab))
 
 
 richness_phyloseq <- function(abund.tab) {
	OTU=otu_table(abund.tab,taxa_are_rows=TRUE)
	richness=estimate_richness(OTU)
	rownames(richness)=colnames(OTU)
	return(richness)
 }
 
 
 
 #' Function for statistical summary of Abundance Table (OTU,GENE,KO,e tc) using phyloseq  
 #'
 #' @param abund.tab otu table object from read.csv/table("otu.tsv"))
 #' @return Statistical summary  for each sample-id
 #' @examples
 #' otu_sp_stat=SampleStatSummary(OTUDATA$otu.tab)
 #' ko_sp_stat=SampleStatSummary(round(KODATA$ko.tab))
 
 SampleStatSummary<- function(abund.tab) {
 
 
  #sum=colSums(otu.tab)
  #mean=colMeans(otu.tab)
  summary=sapply(abund.tab, function(x) c( "SD" = sd(x), 
                          "Mean"= mean(x,na.rm=TRUE),
                           "Sum" =sum(x),
                          "Median" = median(x),
                          "CV" = sd(x)/mean(x,na.rm=TRUE),
                          "Minimum" = min(x),
                          "Maximun" = max(x),
                          "UpperQuantile" = quantile(x,1),
                          "LowerQuartile" = quantile(x,0)
                     )
   )
  summary=t(summary)
  richness=NULL
  try( 
 {
	richness=richness_phyloseq(abund.tab) 
  }
   )
 
  if (!is.null(richness)) {  
	summary=cbind(as.data.frame(richness),as.data.frame(summary)) 
  }
  return(summary)
 }
 
 #' Function for auto-generate statistical summary of any abudance table (OTU,GENE,KO...)
 #' @param abund.tab as.matrix/as.data.frame: otu,ko, rownames=OTU/GENE/KO names, colnames=sample name
 #' @param sp.meta sample metadata as.matrix (ID,NICKNAME,CONDITION,SUBJECT,CATEGORY)
 #' @examples
 #' AbundStatSummary(OTUDATA$otu.tab,OTUDATA$sp.meta)
 
 AbundStatSummary<- function(abund.tab,sp.meta) {
 
	listall1=rownames(sp.meta)
	listall2=colnames(abund.tab)
  	listall=intersect(listall1,listall2)
  	abund.tab=abund.tab[,listall]
 
  	sum=rowSums2(as.matrix(abund.tab))
  	prevalent=rowSums2(as.matrix(abund.tab)>0)
 	prevalent_pct=(prevalent/ncol(as.matrix(abund.tab)))*100
  	mean=rowMeans2(as.matrix(abund.tab))
  	sd=rowSds(as.matrix(abund.tab))
  	cv=rowSds(as.matrix(abund.tab))/rowMeans2(as.matrix(abund.tab),na.rm=TRUE)
  	min=rowMins(as.matrix(abund.tab))
  	max=rowMaxs(as.matrix(abund.tab))
  	summary=cbind(as.data.frame(sum),as.data.frame(prevalent),as.data.frame(prevalent_pct),as.data.frame(mean),as.data.frame(sd),as.data.frame(cv),as.data.frame(min),as.data.frame(max))
  	groups=levels(sp.meta$CONDITION)
  	if (length(groups)>1) {
  		for (g in groups) {
 		print (g)
   		list=rownames(sp.meta[sp.meta$CONDITION==g,])
 		sum_g=rowSums2(as.matrix(abund.tab[,list]))
  		prevalent_g=rowSums2(as.matrix(abund.tab[,list])>0)
  		mean_g=rowMeans2(as.matrix(abund.tab[,list]))
  		sd_g=rowSds(as.matrix(abund.tab[,list]))
  		cv_g=rowSds(as.matrix(abund.tab[,list]))/rowMeans2(as.matrix(abund.tab[,list]),na.rm=TRUE)
  		min_g=rowMins(as.matrix(abund.tab))
  		max_g=rowMaxs(as.matrix(abund.tab))
 		tab=cbind(as.data.frame(sum_g),as.data.frame(prevalent_g),as.data.frame(mean_g),as.data.frame(sd_g),as.data.frame(cv_g),as.data.frame(min_g),as.data.frame(max_g))
 		colnames(tab)=c(paste0("sum_",g),paste0("prevalent_",g),paste0("mean_",g),paste0("sd_",g),paste0("cv_",g),paste0("min_",g),paste0("max_",g))
 		summary=cbind(summary,as.data.frame(tab))
  		}
  	}
 
  	rownames(summary)=rownames(abund.tab)
  	return(summary)
 }
 
 
 
 
 #' Function for differential abundance  analysis for 2 groups  
 #' @param abund.tab as.matrix/as.data.frame: otu,ko, rownames=OTU/GENE/KO names, colnames=sample name
 #' @param sp.meta sample metadata as.matrix (ID,NICKNAME,CONDITION,SUBJECT,CATEGORY)
 #' @param groups vector of 2 biological condition/groups , ex. c('disease','healthy') 
 #' @examples
 #' AbundDiffAnalysis(OTUDATA$otu.tab,OTUDATA$sp.meta,groups)
 
 AbundDiffAnalysis<- function(abund.tab,sp.meta,groups) {
 
	## select only sample names avalaible in metadata and aband table.
 
	listall1=rownames(sp.meta)
	listall2=colnames(abund.tab)
	listall=intersect(listall1,listall2)
	abund.tab=abund.tab[,listall]
 
	## get abundance table for only samples in selected conditions(groups)
	ab=sp.meta[sp.meta$CONDITION %in% groups,]
	#drop unused levels 
	ab=droplevels(ab)
 
	abundab=abund.tab[,rownames(ab)]
 
	## AbundStat of only 2 selecyed conditions
 	stat=AbundStatSummary(abundab,ab)
 
 	lista=rownames(ab[ab$CONDITION==groups[1],])
	listb=rownames(ab[ab$CONDITION==groups[2],])
 
	abunda=as.matrix(abundab[,lista])
	abundb=as.matrix(abundab[,listb])
	mean_a=stat[,paste0("mean_",groups[1])]
	mean_b=stat[,paste0("mean_",groups[2])]

	cv_a=stat[,paste0("cv_",groups[1])]
	cv_b=stat[,paste0("cv_",groups[2])]
         
	# these varibales for MAplot and Valcano plot using 'ggmaplot'
	baseMean=(mean_a+mean_b)/2        
	fc=mean_a/(mean_b+0.000000001)
	log2FoldChange=log2(fc)

     
	fc_cv=cv_a/(cv_b+0.000000001)
	log2FoldChange_cv=log2(fc_cv)


 
	## wilcox compare: get p value, p value adj
 
 
	##OTU/KO ID LIST for testing
         ID.LIST=rownames(stat)
 
        abundab=cbind(as.data.frame(t(abundab)),as.data.frame(ab))
 
	pvalue=c()
 	padj=c()
 	i=1
  	for (ID in ID.LIST ) {
 	eq=paste0(ID,"~CONDITION")
 	cpm=compare_means(as.formula(eq),data=abundab)
	 # pvalue and pajd variable for MAplot and Valcano plot
 	if (length(cpm$p)>0) {
 	pvalue[i]=cpm$p
 	padj[i]=cpm$p.adj
 	} else 
 	{
 	pvalue[i]=1
 	padj[i]=1
 	}
 		i=i+1
 }
 #pvalue=as.matrix(pvalue)
 #colnames(pvalue)="pvalue"
 #padj=as.matrix(padj)
 #colnames(padj)="padj"
 
 diff.tab=cbind(as.data.frame(stat),as.data.frame(baseMean),as.data.frame(fc),as.data.frame(log2FoldChange),as.data.frame(pvalue),as.data.frame(padj),as.data.frame(fc_cv),as.data.frame(log2FoldChange_cv))
 
 return(diff.tab)
 }
