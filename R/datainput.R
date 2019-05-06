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
 