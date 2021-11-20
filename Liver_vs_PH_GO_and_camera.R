edgeR_data <- list(Liver = NA,
                   Cell_susp = NA,
                   Prim_hep = NA,
                   Liv_vs_CS_WT = NA,
                   Liv_vs_PH_WT = NA,
                   CS_vs_PH_WT = NA,
                   Liv_vs_CS_KO = NA,
                   Liv_vs_PH_KO = NA,
                   CS_vs_PH_KO = NA)
  for (i in 1:9){
    edgeR_data[[i]]<- openxlsx::read.xlsx(here("data/LiverVsPH/edgeR.xlsx"),sheet = i)
  }
cpm_matrix <- openxlsx::read.xlsx(here::here("data/LiverVsPH/CPM_matrix.xlsx"))


metadata <- openxlsx::read.xlsx(here::here("data-raw/LiverVsPH/metadata.xlsx"))
metadata <- metadata %>%
  dplyr::filter(!Sample == "558L")
####GO CC on signficant genes####
go_sig_genes_CC <- goAnalysis_PH(edgeR_data, "CC","")

for (i in 1:length(go_sig_genes_CC)){
  go_sig_genes_CC[[i]] <- clusterProfiler::setReadable(go_sig_genes_CC[[i]], OrgDb = org.Mm.eg.db, keyType="ENTREZID")
}
saveRDS(go_sig_genes_CC, file = here::here("data/go_sig_CC.data.rds"))
test_CC <- readRDS(here::here("data/go_sig_CC.data.rds"))
#The following figure was used as figure 9H
Liver <- dotplot(test_CC[[1]], font.size = 14)+
  ggtitle("Liver - Genotype effect")

tiff("LiverGO.tif", units = "cm", width = 25, height = 15, res = 300)
Liver
dev.off()
#The following figure was used as figure 9I
PH <- dotplot(test_CC[[3]], font.size = 14)+
  ggtitle("Primary Hepatocyte - Genotype effect")

tiff("PHGO.tif", units = "cm", width = 25, height = 15, res = 300)
PH
dev.off()

#Run analysis for biological processes
go_sig_genes <- goAnalysis_PH(edgeR_data, "BP", "")

for (i in 1:length(go_sig_genes)){
  go_sig_genes[[i]] <- clusterProfiler::setReadable(go_sig_genes[[i]], OrgDb = org.Mm.eg.db, keyType="ENTREZID")
}
saveRDS(go_sig_genes, file = here::here("data/LiverVsPH/go_sig.data_BP.rds"))

#Generate heatmap for mitochondrial terms in BP analysis
test <- readRDS(here::here("data/LiverVsPH/go_sig.data_BP.rds"))
mito_genes <- test[[1]]@result

gene_list <- mito_genes %>%
  filter(Description == "cellular respiration"|Description =="NADH dehydrogenase complex assembly") %>%
  dplyr::select(geneID)

gene_list_unsplit <- c(unlist(str_split(gene_list[1,], "/")),unlist(str_split(gene_list[2,], "/")))

gene_list_unique <- unique(gene_list_unsplit)
cpm_key <-   clusterProfiler::bitr(
  gene_list_unique,
  fromType = "SYMBOL",
  toType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db"
)

cpm_annotated <- as.data.frame(cpm_matrix)
cpm_annotated <- cpm_annotated %>%
  dplyr::filter(ENSEMBL %in% cpm_key$ENSEMBL)
rownames(cpm_annotated)<-cpm_annotated$ENSEMBL
cpm_annotated <- cpm_annotated |>
  dplyr::select(- ENSEMBL, -SYMBOL)

columns <- colnames(cpm_annotated)
column_order <- c("544L", "548L", "554L", "556L", "564L", "540L", "542L", "546L", "552L", "562L", "566L", "544CS", "548CS", "554CS", "556CS", "558CS", "564CS", "540CS", "542CS", "546CS", "552CS", "562CS", "566CS", "544PH", "548PH", "554PH", "556PH", "558PH", "564PH", "540PH", "542PH", "546PH", "552PH", "562PH",  "566PH")

cpm_test <- cpm_annotated %>%
  dplyr::select(column_order)

conv <- clusterProfiler::bitr(rownames(cpm_test),
                                            fromType = "ENSEMBL",
                                            toType = "SYMBOL",
                                            OrgDb = "org.Mm.eg.db")
rownames(cpm_test) <- conv$SYMBOL
#generate and organize metadata for heatmap
meta_heat_map <- metadata %>%
 dplyr::arrange(match(Sample, column_order)) %>%
  dplyr::filter(!Sample == "558L") %>%
  dplyr::select(Sample, Group)
rownames(meta_heat_map)<-meta_heat_map$Sample
meta_heat_map <- meta_heat_map %>%
  dplyr::select(-Sample)
meta_heat_map <- meta_heat_map %>%
  dplyr::mutate(Group = case_when(
    Group == "WT_L"~"WT L",
    Group =="WT_CS"~"WT CS",
    Group == "WT_PH"~"WT PH",
    Group =="HNKO_L"~"HNKO L",
    Group =="HNKO_CS"~"HNKO CS",
    Group =="HNKO_PH"~"HNKO PH"
  ))

meta_heat_map$Group <- factor(meta_heat_map$Group, levels = c("WT L", "HNKO L", "WT CS", "HNKO CS", "WT PH", "HNKO PH"))
#This figure is displayed as figure 9J
OxPhos <- pheatmap(cpm_test,
         treeheight_col = 0,
         treeheight_row = 0,
         scale = "row",
         legend = T,
         na_col = "white",
         Colv = NA,
         na.rm = T,
         cluster_cols = F,
         fontsize_row = 8,
         fontsize_col = 11,
         cellwidth = 12,
         cellheight = 7,
         annotation_col = meta_heat_map
)

 tiff("Oxphos.tif", units = "cm", width = 25, height = 20, res = 300)
 OxPhos
 dev.off()


 #Code for generating upsetplot supporting fig. 5
 sig_genes <- list("Liver" = NA,
                   "Cell Suspension" = NA,
                   "Primary Hepatocytes" = NA)
 for (i in 1:3){
   sig_genes[[i]]<- edgeR_data[[i]]
   sig_genes[[i]]<- sig_genes[[i]] %>%
     dplyr::filter(FDR < 0.05)
   sig_genes[[i]]<-sig_genes[[i]]$SYMBOL
 }


 order_upset <- c("Liver", "Cell Suspension", "Primary Hepatocytes")

 upsetPlot <- UpSetR::upset(UpSetR::fromList(sig_genes),
                            sets = order_upset,
                            order.by = "freq",
                            keep.order = T,
                            text.scale = c(3, 2, 2, 1.2, 3, 3),
                            point.size = 4
 )

 tiff("UpsetPH.tif", units = "cm", width = 25, height = 20, res = 300)
 upsetPlot
 grid::grid.text("Main effect of Genotype", x=0.8, y =0.98, gp=grid::gpar(fontsize = 20))
 dev.off()
