
 
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
         
 # these varibales for MAplot and Valcano plot using 'ggmaplot'
 baseMean=(mean_a+mean_b)/2        
 fc=mean_a/(mean_b+0.000000001)
 log2FoldChange=log2(fc)
 
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
 
 diff.tab=cbind(as.data.frame(stat),as.data.frame(baseMean),as.data.frame(fc),as.data.frame(log2FoldChange),as.data.frame(pvalue),as.data.frame(padj))
 
 return(diff.tab)
 }
