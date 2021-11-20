####data import and trimming ####

devtools::load_all()
counts <- count_matrix_assembly("counts.csv.gz")

metadata <- load_metadata("022_Metadata.xlsx")
metadata <- metadata %>% dplyr::mutate(Group = paste(Genotype, Age, sep = "_"))

#Quality_control_plots(counts, metadata)


metadata <- metadata %>%
    dplyr::filter(!Basespace.ID == "022_17")
counts <- counts %>%
    dplyr::select(- "022_17")

design <- Generate_design_matrix(metadata)

ctrsts <- limma::makeContrasts(
    Genotype_3d = KO_3d - WT_3d,
    Genotype_6d = KO_6d - WT_6d,
    Genotype_12d = KO_12d - WT_12d,
    Genotype_21d = KO_21d - WT_21d,
    levels = design)

counts <- counts %>% dplyr::select(metadata$Basespace.ID)
all(metadata$Basespace.ID == colnames(counts))

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

#write.xlsx(dgeResults_annotated, file = here("data/edgeR.xlsx"), asTable = TRUE)


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
#write.xlsx(cpm_matrix_anno, file = here("data/CPM_matrix.xlsx"), asTable = TRUE)
