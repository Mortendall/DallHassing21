devtools::load_all()

#load proteomics data
Proteomics_dataset <- openxlsx::read.xlsx(here::here("data-raw/liverProteomics/Proteomics dataset.xlsx"))
setup_proteomics <- openxlsx::read.xlsx("data-raw/liverProteomics/Setup.xlsx")
setup_proteomics <- setup_proteomics |>  dplyr::filter(setup_proteomics$sample %in% colnames(Proteomics_dataset))
setup_proteomics <- setup_proteomics %>%
    tidyr::unite(group, c("Genotype", "Time"), remove = F)
data.table::setDT(setup_proteomics)
setup_proteomics <- setup_proteomics %>%
    dplyr::mutate(ID = as.character(ID))
data.table::setkey(setup_proteomics, ID)
setup_proteomics <- setup_proteomics %>%
    dplyr::mutate(sample = ID)%>%
    dplyr::select(-ID)
#normalize
normalized_proteomics <- limma::normalizeBetweenArrays(log(as.matrix(Proteomics_dataset[,-c(1:2)])), method = "quantile")
rownames(normalized_proteomics) <- Proteomics_dataset$Protein_ID


#remove rows with more than two missing samples pr group
missingSamples_proteomics <- data.table::data.table(is.na(normalized_proteomics), keep.rownames = TRUE) %>%
    data.table::melt(measure.vars = colnames(normalized_proteomics), variable.name = "sample")
missingSamples_proteomics <- merge(setup_proteomics, missingSamples_proteomics, by = "sample")
data.table::setnames(missingSamples_proteomics, "rn", "Accession")
missingSamples_proteomics <- missingSamples_proteomics %>%
    dplyr::group_by(Accession, group) %>%
    dplyr::mutate(nMissing = sum(value))%>%
    reshape2::dcast(Accession ~ group, value.var = "nMissing", fun.aggregate = mean)

setup_proteomics <- setup_proteomics %>%
    dplyr::group_by(group)

cutoff <- 2

tooManyMissing_proteomics <- missingSamples_proteomics %>%
    dplyr::filter(KO_12 > cutoff|
               KO_21 > cutoff |
               KO_3 > cutoff |
               KO_6 > cutoff |
               WT_3 > cutoff|
               WT_6 > cutoff |
               WT_12 > cutoff |
               WT_21 > cutoff)

normalized_proteomics_res <- normalized_proteomics[!(rownames(normalized_proteomics) %in% tooManyMissing_proteomics$Accession), ]

#ensure that setup and colnames are matched
all(colnames(normalized_proteomics_res) == setup_proteomics$sample)

#make mdsPlots to look for outliers
mdsData_proteomics <- limma::plotMDS(normalized_proteomics_res,  plot = FALSE)
mdsData_proteomics <- as.data.frame(mdsData_proteomics)
mdsData_proteomics <- cbind.data.frame(setup_proteomics, mdsData_proteomics)

mdsData_proteomics <-
    mdsData_proteomics %>% data.table::as.data.table() %>%
    dplyr::select(sample, group, eigen.vectors.1, eigen.vectors.2, eigen.vectors.3)

data.table::setnames(mdsData_proteomics,
         c("eigen.vectors.1", "eigen.vectors.2", "eigen.vectors.3", "sample", "group"),
         c("V1", "V2", "V3", "ID", "Group")
         )

ggplot2::ggplot(mdsData_proteomics, aes(x = V1, y = V2, colour = Group)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")+
    ggplot2::geom_point(size = 5) +
    ggplot2::theme_bw()

ggplot2::ggplot(mdsData_proteomics, aes(x = V1, y = V3, colour = Group)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")+
    ggplot2::geom_point(size = 5) +
    ggplot2::theme_bw()

ggplot2::ggplot(mdsData_proteomics, aes(x = V2, y = V3, colour = Group)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")+
    ggplot2::geom_point(size = 5) +
    ggplot2::theme_bw()

#create design matrix and contrast matrix
design_proteomics <- stats::model.matrix(~ 0 + group, setup_proteomics)

colnames(design_proteomics) <- stringr::str_remove_all(colnames(design_proteomics), "group")
fit_proteomics <- limma::lmFit(normalized_proteomics_res, design = design_proteomics, method = "robust")
cont.matrix_proteomics <- limma::makeContrasts(
    Genotype_3 = KO_3 - WT_3,
    Genotype_6 = KO_6 - WT_6,
    Genotype_12 = KO_12 - WT_12,
    Genotype_21 = KO_21 - WT_21,
    HNKO_age_21_3 = KO_21 - KO_3,
    WT_age_21_3 = WT_21 - WT_3,
    Interaction = (KO_21 - WT_21) - (KO_3 - WT_3),
    levels = design_proteomics)

#make models and create data table for export
fit2_proteomics <- limma::contrasts.fit(fit_proteomics, cont.matrix_proteomics)
fit2_proteomics <- limma::eBayes(fit2_proteomics, trend = TRUE, robust = TRUE)
resultTables_proteomics <- list(
    Genotype_3 = limma::topTable(fit2_proteomics, coef = "Genotype_3", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    Genotype_6 = limma::topTable(fit2_proteomics, coef = "Genotype_6", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    Genotype_12 = limma::topTable(fit2_proteomics, coef = "Genotype_12", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    Genotype_21 = limma::topTable(fit2_proteomics, coef = "Genotype_21", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    HNKO_age_21_3 = limma::topTable(fit2_proteomics, coef = "HNKO_age_21_3", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    WT_age_21_3 = limma::topTable(fit2_proteomics, coef = "WT_age_21_3", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    Interaction = limma::topTable(fit2_proteomics, coef = "Interaction", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE)
)

lapply(resultTables_proteomics, data.table::setnames, "rn", "Accession")
conv_prot <- Proteomics_dataset[, 1:2]
colnames(conv_prot)[2] <- "Accession"
conv_prot <- data.table::as.data.table(conv_prot)
data.table::setkey(conv_prot, Accession)

for (i in resultTables_proteomics){
    i[, Gene:=conv_prot[Accession, Gene_names]]
}
openxlsx::write.xlsx(resultTables_proteomics, file = here::here("data/limma_results_proteomics_time_course_liver.xlsx"))

#####GO analysis and upsetplot#####
#these figures were used in figure 2.
#layout is created in file GOAnalysis_Figures.R
resultTables_proteomics <- list("HNKO 3d" = NA,
                                "HNKO 6d" = NA,
                                "HNKO 12d" = NA,
                                "HNKO 21d" = NA
)
for (i in 1:4){
    resultTables_proteomics[[i]] <- openxlsx::read.xlsx(here::here("data/limma_results_proteomics_time_course_liver.xlsx"),i)
}

proteomics_sig <- resultTables_proteomics
for (i in 1:4){
    proteomics_sig[[i]]<-proteomics_sig[[i]] %>%
        dplyr::filter(adj.P.Val < 0.05)
    proteomics_sig[[i]]<-proteomics_sig[[i]]$Gene
}

order_upset <- c("HNKO 21d", "HNKO 12d", "HNKO 6d","HNKO 3d")

upset_data <- data.frame("Gene"= proteomics_sig[[1]],
                         "Group"= names(proteomics_sig[1]))
for (i in 2:length(proteomics_sig)){
    upset_data <- dplyr::add_row(upset_data,"Gene"=proteomics_sig[[i]],
                                 "Group"=names(proteomics_sig)[i])
}
upset_data <- dplyr::distinct(upset_data)
upset_data <- upset_data |>
    dplyr::mutate("TRUE"=TRUE)
upset_wide <- tidyr::pivot_wider(upset_data, names_from = Group, values_from = "TRUE",values_fill = FALSE) |>
    dplyr::filter(!is.na(Gene))

rownames(upset_wide)<-upset_wide$Gene

upset_wide <- dplyr::select(upset_wide,-Gene)
upsetPlot <- ComplexUpset::upset(upset_wide, order_upset, name = "", sort_sets = F, themes = ComplexUpset::upset_modify_themes(list(
    'intersections_matrix'=theme(text = element_text(size = 16)))))
#ggplotify to use the object in patchwork

upsetPlot <- ggplotify::as.ggplot(upsetPlot)
upsetPlot <- upsetPlot + ggplot2::ggtitle("Proteins with effect of genotype")+ggplot2::theme(plot.title = element_text(size = 18,
                                                                                                                       hjust = 0.5,
                                                                                                                       vjust = 0.95))
#save RDS object to load into RNAseq analysis
saveRDS(upsetPlot, here::here("data/proteinUpset.RDS"))


#GO for overlap genes
overlap_genes <- as.data.frame(proteomics_sig[[1]]) %>%
    dplyr::filter(proteomics_sig[[1]]%in% proteomics_sig[[2]] &
                      proteomics_sig[[1]]%in% proteomics_sig[[3]]  &
                      proteomics_sig[[1]]%in% proteomics_sig[[4]]
    )
overlap_genes <- overlap_genes %>%
    dplyr::filter(!is.na(overlap_genes))
colnames(overlap_genes) <- "Symbol"

overlap_entrez <- bitr(overlap_genes$Symbol,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = "org.Mm.eg.db",
                       drop = T)


overlap_bg = bitr(resultTables_proteomics[[1]]$Gene,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = "org.Mm.eg.db",
                  drop = T)

goResults_overlap <- enrichGO(gene = overlap_entrez$ENTREZID,
                              universe = overlap_bg$ENTREZID,
                              OrgDb = org.Mm.eg.db,
                              ont = "MF")
protGO <- enrichplot::dotplot(goResults_overlap)+ggplot2::ggtitle("Overlap Proteins GO-terms")+ggplot2::theme(plot.title = element_text(size = 18, hjust = 0.5))
saveRDS(protGO, here::here("data/protGO.RDS"))
tiff("GOProt_MF.tif", unit = "cm", height = 10, width = 25, res = 300)
protGO
dev.off()
