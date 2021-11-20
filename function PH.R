library(fs)
library(vroom)
library(Biobase)
library(here)
library(magrittr)
library(data.table)
library(clusterProfiler)
library(edgeR)
library(openxlsx)
library(pheatmap)
library(gridExtra)
library(PoiClaClu)
library(RColorBrewer)
library(limma)
library(GO.db)
library(stringr)
library(dplyr)
library(ggplot2)
library(ComplexUpset)
#Functions for analysis of the liver vs primary hepatocyte study. Some things are slightly differnet in this setup, and hence functions were written specifically for this

#' count_matrix_loader
#'
#' @param file_type your raw data file type (presently optimized for txt)
#'
#' @return a count matrix containing raw count values for the experiment with group annotations.

count_matrix_assembly_PH <- function(file_type){
  count_file <- fs::dir_ls(here::here("data-raw/LiverVsPH/"),
                           regexp = file_type,
                           recurse = TRUE)
  count_matrix <- openxlsx::read.xlsx(xlsxFile = count_file,sheet = 1)
  count_matrix$Geneid<-stringr::str_remove_all(count_matrix$Geneid, "\\..*")
  rownames(count_matrix) <- count_matrix$Geneid
  count_matrix <- count_matrix %>%
    dplyr::select(-Geneid)

  return(count_matrix)
}

#' Load metadata and sort them according to count matrix input
#'
#' @param file_name the name of teh metadata file (default "metadata.csv")
#'
#' @return metadata file sorted according to count matrix order

 load_metadata_PH <- function(file_name) {
   data_file <- fs::dir_ls(here::here("data-raw/LiverVsPH/"),
                           regexp = file_name,
                           recurse = T)
   metadata <- openxlsx::read.xlsx(xlsxFile = data_file)
   return(metadata)
 }



#' RNAseq_processing
#'
#' @param count_matrix generated with the count matrix assembly function
#' @param metadata loaded in with the load_metadata function
#' @param design Generated with Generate_design_matrix function
#' @param ctrsts defined in the WATanalysis script
#'
#' @return a dgeResults list object

RNAseq_processing_PH <- function(count_matrix, metadata, design, ctrsts) {
  group <- as.matrix(metadata$Group)
  RNAseq <- edgeR::DGEList(counts = count_matrix, group = group)
  keep <- edgeR::filterByExpr(RNAseq, design = design, min.count = 20)
  RNAseq <- RNAseq[keep, , keep.lib.sizes = F]
  RNAseq <- edgeR::calcNormFactors(RNAseq)
  RNAseq <- edgeR::estimateDisp(RNAseq,design)
  efit <- edgeR::glmQLFit(RNAseq, design)
  dgeResults <- apply(ctrsts, 2, . %>%
                        edgeR::glmQLFTest(glmfit = efit, contrast = .) %>%
                        edgeR::topTags(n = Inf, p.value = 1) %>%
                        magrittr::extract2("table") %>%
                        data.table::as.data.table(keep.rownames = TRUE))
  return(dgeResults)
}

#' Generate design matrix
#'
#' @param metadata a metadata object generated through the load_metadata function
#'
#' @return a design matrix file


Generate_design_matrix_PH <- function(metadata){
  design <- stats::model.matrix( ~0+Group, metadata)
  colnames(design) <-
    stringr::str_remove_all(colnames(design), "\\(|\\)|Group|:")
  return(design)
}


#' Gene ontology enrichment analysis of genes generated from a results file
#'
#' @param result_list list of data.tables generated from edgeR. Must be data.table and contain a SYMBOL annotation column
#'
#' @return a list containing enrichresults for each element in the results file list

goAnalysis_PH <- function(result_list, ontology, direction){
  bg <- result_list[[1]]
  bg_list <- clusterProfiler::bitr(
    bg$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
  )

  goResult_list <- vector(mode = "list", length = length(result_list))
  if (direction == "Upregulated"){
    print("Running analysis on upregulated genes")
    for(i in 1:length(result_list)){
      sig_list<- result_list[[i]] %>%
        dplyr::filter(FDR<0.05) |>
        dplyr::filter(logFC>0)
      eg <- clusterProfiler::bitr(
        sig_list$SYMBOL,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
      )
      goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                             universe = bg_list$ENTREZID,
                                             OrgDb = org.Mm.eg.db,
                                             ont = ontology)
      goResult_list[[i]]<- goResults
    }


  }

  if(direction == "Downregulated"){
    print("Running analysis on downregulated genes")
    for(i in 1:length(result_list)){
      sig_list<- result_list[[i]] %>%
        dplyr::filter(FDR<0.05) |>
        dplyr::filter(logFC<0)
      eg <- clusterProfiler::bitr(
        sig_list$SYMBOL,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
      )
      goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                             universe = bg_list$ENTREZID,
                                             OrgDb = org.Mm.eg.db,
                                             ont = ontology)
      goResult_list[[i]]<- goResults
    }

  }
  if (direction == ""){
    for(i in 1:length(result_list)){
      sig_list<- result_list[[i]] %>%
        dplyr::filter(FDR<0.05)
      eg <- clusterProfiler::bitr(
        sig_list$SYMBOL,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
      )
      goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                             universe = bg_list$ENTREZID,
                                             OrgDb = org.Mm.eg.db,
                                             ont = ontology)
      goResult_list[[i]]<- goResults
    }
  }


  for (i in 1:length(goResult_list)){
    names(goResult_list)[i]<-names(result_list)[i]
  }
  return(goResult_list)

}



