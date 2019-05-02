#library(ggfortify)
#library(stringr)

library(pheatmap)
library(ape)
library(ggplot2)
library(phyloseq)
library(matrixStats)
#vignette(phyloseq-analysis)
#vignette("phyloseq-basics")
#theme_set(theme_bw())


################ PART A: INPUT FILE   ###########################



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

AbundPcaPlot = function(abund.tab,sp.meta, meanCountCutOff=1,groups="ALL",title = "Title"){
  

  otu.mean <- apply(X = abund.tab, MARGIN = 1, FUN = mean)
  otu.filter <- abund.tab[otu.mean >= meanCountCutOff,]
  otu=t(otu.filter)


  if (groups=="ALL") 
  {
  	condition=sp.meta[,2]
  	pca=prcomp(log(otu+1))
  	data_df <- data.frame(pc_1 = pca$x[,1], pc_2 = pca$x[,2], c = condition)
  }
  else
  {
    	temp=sp.meta[sp.meta$CONDITION %in% groups,]
	mylist=rownames(temp)
	mycondition=temp[,2]
	otu=otu[mylist,]
	pca=prcomp(log(otu+1))
  	data_df <- data.frame(pc_1 = pca$x[,1], pc_2 = pca$x[,2], c = mycondition)
	
  }
  p <-ggplot(data_df, aes(pc_1,pc_2)) + geom_point(aes(color=c,shape=c)) + labs(title = title)
  p
}



#' Function to create a heatmap  of abundance data (ROWs:OTU,GENE,KO, COLS: SAMPLE NAMES)  
#'
#' @param abund.tab abundance table as.matrix
#' @param sp.meta sample metadata as.matrix (ID,NICKNAME,CONDITION,SUBJECT,CATEGORY)
#' @param groups list of biological condition(Healthy vs Disease)
#' @return Heatmap plot
#' @examples
#' OtuPlotHeatmap(OTUDATA$otu.tab,OTUDATA$sp.meta, meanCountCutOff = 50, groups=c('protein','cholic'),fontsize_row=5,fontsize_col=5)

AbundPlotHeatmap <- function(abund.tab,sp.meta, meanCountCutOff = 1, groups="ALL",fontsize_row=5,fontsize_col=5) {

  counts <- abund.tab
  # Remove counts with too low count
  counts.mean <- apply(X=counts, MARGIN = 1, FUN = mean)
  expressed <- counts.mean >= meanCountCutOff
  counts <- counts[expressed,]

  log2.counts <- log2(counts + 1)

  if (groups=="ALL") 
  {
  	annot <- data.frame(Group = sp.meta$CONDITION)
  	rownames(annot) <- rownames(sp.meta)
  }
  else
  {
	temp=metadata[metadata$CONDITION %in% groups,]
	mylist=rownames(temp)
	annot <- data.frame(Group =temp[,2])
	rownames(annot) <- rownames(temp)
        log2.counts=log2.counts[,mylist]
	
  }
  pheatmap(log2.counts, scale="row", show_rownames = F, show_colnames = T, annotation_col = annot,fontsize_row=fontsize_row,fontsize_col=fontsize_col)

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
	{richness=richness_phyloseq(abund.tab) 
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








