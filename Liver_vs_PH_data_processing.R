####data import and trimming ####


counts <- count_matrix_assembly_PH("count_matrix.xlsx")

metadata <- load_metadata_PH("metadata.xlsx")


#Quality_control_plots(counts, metadata)
#QC reveals that 558L looks very odd. Remove and rerun. QCplots saved as QCplots_before_filtering
metadata <- metadata %>%
  dplyr::filter(!Sample == "558L")
counts <- counts %>%
  dplyr::select(- "558L")

Quality_control_plots(counts, metadata)
#MDS plots from post filtering was used for Figure 9G
design <- Generate_design_matrix_PH(metadata)

ctrsts <- limma::makeContrasts(
  Liver = HNKO_L - WT_L,
  Cell_susp = HNKO_CS - WT_CS,
  Prim_hep = HNKO_PH - WT_PH,
  Liv_vs_CS_WT = WT_L - WT_CS,
  Liv_vs_PH_WT = WT_L - WT_PH,
  CS_vs_PH_WT = WT_CS - WT_PH,
  Liv_vs_CS_KO = HNKO_L - HNKO_CS,
  Liv_vs_PH_KO = HNKO_L - HNKO_PH,
  CS_vs_PH_KO = HNKO_CS - HNKO_PH,
  levels = design)

all(metadata$Sample == colnames(counts))

#DGE analysis
group <- as.matrix(metadata$Group)
RNAseq <- edgeR::DGEList(counts = counts, group = group)
keep <- edgeR::filterByExpr(RNAseq, design = design, min.count = 20)
RNAseq <- RNAseq[keep, , keep.lib.sizes = F]
RNAseq <- edgeR::calcNormFactors(RNAseq)
counts_norm <- RNAseq$counts
cpm_matrix <- cpm(RNAseq, log = T)
key <- clusterProfiler::bitr(rownames(cpm_matrix), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
RNAseq <- edgeR::estimateDisp(RNAseq,design)
efit <- edgeR::glmQLFit(RNAseq, design)
dgeResults <- apply(ctrsts, 2, . %>%
                      edgeR::glmQLFTest(glmfit = efit, contrast = .) %>%
                      edgeR::topTags(n = Inf, p.value = 1) %>%
                      magrittr::extract2("table") %>%
                      data.table::as.data.table(keep.rownames = TRUE))
dgeResults_annotated <- dgeResults
for (i in 1:length(dgeResults_annotated)){
  data.table::setnames(dgeResults_annotated[[i]],names(dgeResults_annotated[[i]])[1], "ENSEMBL")
  ens2symbol <-
    clusterProfiler::bitr(dgeResults_annotated[[i]]$ENSEMBL,
                          fromType = 'ENSEMBL',
                          toType = 'SYMBOL',
                          OrgDb = "org.Mm.eg.db")
  ens2symbol %<>% data.table %>% data.table::setkey(ENSEMBL)
  dgeResults_annotated[[i]] <- dplyr::full_join(dgeResults_annotated[[i]], ens2symbol)
}

write.xlsx(dgeResults_annotated, file = here("data/LiverVsPH/edgeR.xlsx"), asTable = TRUE)


#Make and save a cpm_matrix with annotations
cpm_matrix_anno <- cpm_matrix
cpm_matrix_anno <- as.data.frame(cpm_matrix_anno)
cpm_matrix_anno <- cpm_matrix_anno %>%
  dplyr::mutate(ENSEMBL = rownames(cpm_matrix_anno))

ens2symbol <-
  clusterProfiler::bitr(cpm_matrix_anno$ENSEMBL,
                        fromType = 'ENSEMBL',
                        toType = 'SYMBOL',
                        OrgDb = "org.Mm.eg.db")
ens2symbol %<>% data.table %>% data.table::setkey(ENSEMBL)
cpm_matrix_anno <- dplyr::full_join(cpm_matrix_anno, ens2symbol)
write.xlsx(cpm_matrix_anno, file = here("data/LiverVsPH/CPM_matrix.xlsx"), asTable = TRUE)

#Code for generating figure 9G
group <- as.matrix(metadata$Group)

RNAseq <- edgeR::DGEList(counts = counts, group = group)

mdsData <- plotMDS(RNAseq, plot = FALSE)
mdsData <-
  mdsData$eigen.vectors %>% as.data.table() %>%
  dplyr::mutate(ID = rownames(RNAseq$samples)) %>%
  dplyr::mutate(Group = metadata$Group) %>%
  dplyr::select(ID, Group, V1, V2, V3)


setnames(mdsData,
         c("V1", "V2", "V3", "ID", "Group"),
         c("dim1", "dim2", "dim3", "ID", "Group"))

pBase <-
  ggplot(mdsData, aes(x = dim1, y = dim2, colour = Group)) +
  geom_point(size = 10) +
  #geom_label(show.legend = FALSE, size = 5) +
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5))+
  ggtitle("MDS Plot")

tiff("MDSplot.tif", units = "cm", width = 20, height = 20, res = 300)
pBase
dev.off()
