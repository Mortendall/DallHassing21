devtools::load_all()

#load proteomics data
Proteomics_dataset <- openxlsx::read.xlsx(here::here("data-raw/plasmaProteomics/Plasma Proteomics dataset sorted.xlsx"))
setup_proteomics <- openxlsx::read.xlsx("data-raw/plasmaProteomics/Setup_plasma.xlsx")
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
#openxlsx::write.xlsx(resultTables_proteomics, file = here::here("data/limma_results_proteomics_time_course.xlsx"))

######Upsetplot#####
#used to generate Figure 2G
significant_proteins_upset <- list("HNKO 3d" = resultTables_proteomics[[1]],
                                   "HNKO 6d" = resultTables_proteomics[[2]],
                                   "HNKO 12d" = resultTables_proteomics[[3]],
                                   "HNKO 21d" = resultTables_proteomics[[4]])
for (i in 1:4){
    significant_proteins_upset [[i]]  <- significant_proteins_upset [[i]] %>%
        dplyr::filter(adj.P.Val < 0.05)
    significant_proteins_upset [[i]]  <- significant_proteins_upset [[i]]$Gene
}

order_upset <- c("HNKO 21d", "HNKO 12d", "HNKO 6d","HNKO 3d")
upset <-UpSetR::upset(UpSetR::fromList(significant_proteins_upset),
                      sets = order_upset,
                      order.by = "freq",
                      keep.order = T,
                      text.scale = 2
)

grid::grid.text("Plasma Proteomics", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 30))


tiff(here::here("data/plasma_upset.tif"), unit = "cm", height = 14, width = 25, res = 300)
upset
grid::grid.text("Plasma Proteomics", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 30))
dev.off()
