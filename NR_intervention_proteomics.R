devtools::load_all()

#load proteomics data
Proteomics_dataset <- data.table::fread(here::here("data-raw/NR_intervention/proteomics/expressions.csv"), header = T)
setup_proteomics <- data.table::fread("data-raw/NR_intervention/proteomics/setup.csv")
data.table::setnames(setup_proteomics, c("ID", "Genotype", "Treatment"))
setup_proteomics[, group:=paste(Genotype, Treatment, sep = "_")]
setup_proteomics <- setup_proteomics |>  dplyr::filter(setup_proteomics$ID %in% colnames(Proteomics_dataset))
data.table::setDT(setup_proteomics)
setup_proteomics <- setup_proteomics %>%
    dplyr::mutate(ID = as.character(ID))
data.table::setkey(setup_proteomics, ID)
setup_proteomics <- setup_proteomics %>%
    dplyr::mutate(sample = ID)%>%
    dplyr::select(-ID)
#normalize
normalized_proteomics <- limma::normalizeBetweenArrays(log(as.matrix(Proteomics_dataset[,-c(1:2)])), method = "quantile")
rownames(normalized_proteomics) <- Proteomics_dataset$Accession

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
    dplyr::filter(KO_Control > cutoff |
                      KO_NR > cutoff |
                      WT_Control > cutoff |
                      WT_NR > cutoff)

normalized_proteomics_res <- normalized_proteomics[!(rownames(normalized_proteomics) %in% tooManyMissing_proteomics$Accession), ]

#ensure that setup and colnames are matched
all(colnames(normalized_proteomics_res) == setup_proteomics$sample)

#make mdsPlots to look for outliers
mdsData_proteomics <- limma::plotMDS(normalized_proteomics_res,  plot = FALSE)
mdsData_proteomics <- as.data.frame(mdsData_proteomics$eigen.vectors)
mdsData_proteomics <- cbind.data.frame(setup_proteomics, mdsData_proteomics)

mdsData_proteomics <-
    mdsData_proteomics %>% data.table::as.data.table() %>%
    dplyr::select(sample, group, V1, V2, V3)

data.table::setnames(mdsData_proteomics,
                     c("sample", "group"),
                     c("ID", "Group")
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

#mouse 330 was steatotic so it is removed
setup_proteomics <- setup_proteomics |>
    dplyr::filter(sample != "330")
#330 is the steatotic mouse
normalized_proteomics_res <- as.data.frame(normalized_proteomics_res)%>%
    dplyr::select(!"330")

design <- stats::model.matrix(~ 0 + group, setup_proteomics)

colnames(design) <- stringr::str_remove_all(colnames(design), "group")
fit <- limma::lmFit(normalized_proteomics_res, design = design, method = "robust")

cont.matrix <- limma::makeContrasts(
    Treatment_in_WT = WT_NR - WT_Control,
    Treatment_in_KO = KO_NR - KO_Control,
    KO_in_Control = KO_Control - WT_Control,
    KO_in_NR = KO_NR - WT_NR,
    Interaction = (KO_NR - KO_Control) - (WT_NR - WT_Control),
    levels = design
)
fit2 <- limma::contrasts.fit(fit, cont.matrix)
fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

resultTables <- list(
    Treatment_in_WT = limma::topTable(fit2, coef = "Treatment_in_WT", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    Treatment_in_KO = limma::topTable(fit2, coef = "Treatment_in_KO", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    KO_in_Control = limma::topTable(fit2, coef = "KO_in_Control", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    KO_in_NR = limma::topTable(fit2, coef = "KO_in_NR", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    Interaction = limma::topTable(fit2, coef = "Interaction", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE)
)

lapply(resultTables, setnames, "rn", "Accession")
conv <- Proteomics_dataset[, 1:2]
data.table::setkey(conv, Accession)

for (i in resultTables){
    i[, Gene:=conv[Accession, Gene]]
}
#openxlsx::write.xlsx(resultTables, file = here::here("data/limma_results_NR_int.xlsx"))

#####GO analysis#####
limma_results <- vector(mode = "list", length = 4)

temporary_result <- data.frame(NA, NA, NA, NA)

for (i in 1:4){
    temporary_result <- openxlsx::read.xlsx(here::here("data/limma_results_NR_int.xlsx"),i)
    limma_results[[i]]<- temporary_result
}

names(limma_results)[[1]]<- "Treatment in WT"
names(limma_results)[[2]]<- "Treatment in KO"
names(limma_results)[[3]]<- "Genotype in control"
names(limma_results)[[4]]<- "Genotype in NR"

#Generate GO-data and save it

GO_results_MF <- vector(mode = "list", length = 4)
temporary_result <- data.frame(NA, NA, NA, NA,NA,NA,NA,NA,NA)

names(GO_results_MF)[[1]]<- "Treatment in WT"
names(GO_results_MF)[[2]]<- "Treatment in KO"
names(GO_results_MF)[[3]]<- "Genotype in control"
names(GO_results_MF)[[4]]<- "Genotype in NR"

bg <- clusterProfiler::bitr(limma_results$`Treatment in WT`$Gene,
                                      fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = "org.Mm.eg.db",
                                      drop = T)
limma_results_sig <- limma_results
for (i in 1:4){
    limma_results_sig[[i]]<-limma_results_sig[[i]] %>%
        dplyr::filter(adj.P.Val < 0.05)
    limma_results_sig[[i]]<-limma_results_sig[[i]]$Gene
}

for (i in 1:4){
    eg <- clusterProfiler::bitr(
        limma_results_sig[[i]],
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )
    goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                           universe = bg$ENTREZID,
                                           OrgDb = org.Mm.eg.db,
                                           ont = "MF")
    goResults <- clusterProfiler::setReadable(goResults, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    GO_results_MF[[i]]<- goResults@result
}
#openxlsx::write.xlsx(GO_results_MF, here::here("data/go_data_NR_MF.xlsx"))

GO_results_CC <- vector(mode = "list", length = 4)
temporary_result <- data.frame(NA, NA, NA, NA,NA,NA,NA,NA,NA)

names(GO_results_CC)[[1]]<- "Treatment in WT"
names(GO_results_CC)[[2]]<- "Treatment in KO"
names(GO_results_CC)[[3]]<- "Genotype in control"
names(GO_results_CC)[[4]]<- "Genotype in NR"

for (i in 1:4){
    eg <- clusterProfiler::bitr(
        limma_results_sig[[i]],
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )
    goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                           universe = bg$ENTREZID,
                                           OrgDb = org.Mm.eg.db,
                                           ont = "CC")
    goResults <- clusterProfiler::setReadable(goResults, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    GO_results_CC[[i]]<- goResults@result
}
#openxlsx::write.xlsx(GO_results_CC, here::here("data/go_data_NR_CC.xlsx"))

#####Effect of NR in HNKO#####
#This code was used to generate figure 7A+B
NR_effect <- limma_results[[2]]
HNKO_effect <- limma_results[[3]]

NR_effect_sig <- NR_effect %>%
    dplyr::filter(adj.P.Val < 0.05)

HNKO_effect_sig <- HNKO_effect %>%
    dplyr::filter(adj.P.Val < 0.05)


NR_entrez <- clusterProfiler::bitr(NR_effect_sig$Gene,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = "org.Mm.eg.db",
                  drop = T)
NR_entrez_bg <- clusterProfiler::bitr(NR_effect$Gene,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = "org.Mm.eg.db",
                     drop = T)

go_Results_NR <- clusterProfiler::enrichGO(gene = NR_entrez$ENTREZID,
                          universe = NR_entrez_bg$ENTREZID,
                          OrgDb = org.Mm.eg.db,
                          ont = "MF")

go_Results_NR_CC <- clusterProfiler::enrichGO(gene = NR_entrez$ENTREZID,
                             universe = NR_entrez_bg$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "CC")
GO_NR <- enrichplot::dotplot(go_Results_NR)+ggtitle("HNKO Con vs NR - MF" )
GO_NR_CC <- enrichplot::dotplot(go_Results_NR_CC)+ggtitle("HNKO Con vs NR - CC")

# tiff(here::here("data/GO_NR_MF.tif"), unit = "cm", height = 20, width = 30, res = 300)
# GO_NR/GO_NR_CC
# dev.off()

#these figures were used in supporting figure 3.
con_entrez <- clusterProfiler::bitr(HNKO_effect_sig$Gene,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Mm.eg.db",
                   drop = T)


go_Results_con <- clusterProfiler::enrichGO(gene = con_entrez$ENTREZID,
                           universe = NR_entrez_bg$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "MF")
go_Results_con_cc <- clusterProfiler::enrichGO(gene = con_entrez$ENTREZID,
                              universe = NR_entrez_bg$ENTREZID,
                              OrgDb = org.Mm.eg.db,
                              ont = "CC")
con_plot <- enrichplot::dotplot(go_Results_con)+ggtitle("Control WT vs HNKO \n Molecular Function" )
con_plot_cc <- enrichplot::dotplot(go_Results_con_cc)+ggtitle("Control WT vs HNKO \n Cellular Component" )

# tiff(here::here("data/GO_con_MF_CC.tif"), unit = "cm", height = 20, width = 30, res = 300)
# con_plot/con_plot_cc
# dev.off()

#####Upsetplot#####
#this was used for figure 6I
limma_results <- vector(mode = "list", length = 4)

temporary_result <- data.frame(NA, NA, NA, NA)

for (i in 1:4){
    temporary_result <- openxlsx::read.xlsx(here::here("data/limma_results_NR_int.xlsx"), i)
    limma_results[[i]]<- temporary_result
}

names(limma_results)[[1]]<- "Treatment in WT"
names(limma_results)[[2]]<- "Treatment in KO"
names(limma_results)[[3]]<- "Genotype in control"
names(limma_results)[[4]]<- "Genotype in NR"
limma_results_sig <- limma_results
for (i in 1:4){
    limma_results_sig[[i]]<-limma_results_sig[[i]] %>%
        dplyr::filter(adj.P.Val < 0.05)
    limma_results_sig[[i]]<-limma_results_sig[[i]]$Gene
}

order_upset <- c("Genotype in NR", "Genotype in control", "Treatment in KO","Treatment in WT")
upset_NR <-UpSetR::upset(UpSetR::fromList(limma_results_sig),
                         sets = order_upset,
                         order.by = "freq",
                         keep.order = T,
                         text.scale = 3.5,
)

grid::grid.text("Significantly altered proteins", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 20))

# tiff(here::here("data/UpsetProtein_NR.tif"), unit = "cm", height = 25, width = 35, res = 300)
# upset_NR
# grid::grid.text("Significantly altered proteins", x=0.7, y = 0.95, gp=grid::gpar(fontsize = 30))
# dev.off()

#####Igraph of oxidoreducaseactivity proteins#####
#This script was used to generate Fig 7C

go_results <- list("Treatment in WT" = NA,
                   "Treatment in KO" = NA,
                   "Genotype in control" = NA,
                   "Genotype in NR" = NA)
for (i in 1:4){
    go_results[[i]]<- openxlsx::read.xlsx(here::here("data/go_data_NR_MF.xlsx"),i)
}

oxphos_proteins <- go_results$`Treatment in KO`$geneID[[1]]
oxphos_proteins <- unlist(str_split(oxphos_proteins, "/"))


g1 <- igraph::graph(oxphos_proteins, isolates = oxphos_proteins)
plot(g1, edge.arrow.size = 0.5)
#load the manual annotation of OXRED proteins
group_anno <- openxlsx::read.xlsx(here::here("data-raw/NR_intervention/Proteomics/annotation_oxRed.xlsx"), colNames = T)

all(group_anno$SYMBOL == oxphos_proteins)

#try to link together the terms that are related by group

group_anno <- group_anno %>%
    dplyr::mutate(color = case_when(
        Group == "NAD" ~ "grey",
        Group == "NADP" ~ "orange",
        Group == "FAD" ~ "red",
        Group == "FMN(H2)" ~"turquoise",
        Group == "Oxidase" ~ "green1",
        Group == "OxPhos" ~ "pink",
        Group == "NAD/OxPhos" ~ "cyan",
        Group == "ROS" ~ "yellow",
        Group == "Other" ~ "green2"
    ))
Group <- group_anno$SYMBOL
color_vector <- group_anno$color
color_vector[73:81]<- "white"
g3 <- igraph::graph.data.frame(group_anno,
                               directed = F)
V(g3)$label.cex = 0.65
 # tiff("data/Network_plot.tif", unit = "cm", height = 40, width = 40, res = 300)
 # set.seed(42)
 # plot(g3, vertex.color = color_vector, vertex.label.cex = 1.5, vertex.size = 15)
 # legend("topleft",legend = unique(group_anno$Group),
 #        pt.bg  = unique(group_anno$color),
 #        pch    = 21,
 #        cex    = 1.5,
 #        bty    = "n",
 #        title  = "Groups")
 # title("Proteins from oxidoreductase activity \n classified by co-factor", cex.main = 2)
 # dev.off()

#####Heatmaps#####
#Heatmap prep
go_results <- list("Treatment in WT" = NA,
                   "Treatment in KO" = NA,
                   "Genotype in control" = NA,
                   "Genotype in NR" = NA)
for (i in 1:4){
    go_results[[i]]<- openxlsx::read.xlsx(here::here("data/go_data_NR_MF.xlsx"),i)
}
#Script was used to generate Supporting Fig. 4A
oxphos_proteins <- go_results$`Treatment in KO`$geneID[[1]]
oxphos_proteins <- unlist(str_split(oxphos_proteins, "/"))


expressions <- data.table::fread(here::here("data-raw/NR_intervention/Proteomics/expressions.csv"), header = TRUE)
setup <- data.table::fread(here::here("data-raw/NR_intervention/Proteomics/setup.csv"))
data.table::setnames(setup, c("sample", "Genotype", "Treatment"))
setup[, group:=paste(Genotype, Treatment, sep = "_")]
setup <- setup[sample != "330"]
#330 is the steatotic mouse
expressions <- expressions %>%
    dplyr::select(!"330")

res <- limma::normalizeBetweenArrays(log(as.matrix(expressions[,-c(1:2)])), method = "quantile")
res <- as.data.frame(res)
res <- res %>%
    dplyr::mutate(Gene = expressions$Gene)


res_ox <- res %>%
    dplyr::filter(Gene %in% oxphos_proteins) %>%
    dplyr::distinct(Gene, .keep_all = T)
rownames(res_ox) <- res_ox$Gene
res_ox <- res_ox %>%
    dplyr::select(-Gene)

setup_ordered <- setup
setup_ordered <- setup_ordered %>%
    dplyr::mutate(
        group = dplyr::case_when(
            group == "WT_Control" ~ "WT Control",
            group == "KO_Control" ~ "HNKO Control",
            group == "WT_NR" ~ "WT NR",
            group == "KO_NR" ~ "HNKO NR"
        )
    )
order <- c("WT Control", "HNKO Control", "WT NR", "HNKO NR")

setup_ordered <- setup_ordered %>%
    dplyr::arrange(Treatment,desc(Genotype))

setup_ordered$sample <- as.character(setup_ordered$sample)
res_ox <- res_ox %>%
    dplyr::select(setup_ordered$sample)

key <- as.data.frame(setup_ordered)

key <- key %>%
    dplyr::select(group)
rownames(key) <- setup_ordered$sample
key$group <- factor(key$group, c("WT Control", "HNKO Control", "WT NR", "HNKO NR"))



OxRed <- pheatmap::pheatmap(res_ox,
                            treeheight_col = 0,
                            treeheight_row = 0,
                            scale = "row",
                            legend = T,
                            na_col = "white",
                            Colv = NA,
                            na.rm = T,
                            cluster_cols = F,
                            fontsize_row = 5,
                            fontsize_col = 8,
                            cellwidth = 7,
                            cellheight = 5,
                            annotation_col = key,
                            show_colnames = F,
                            show_rownames = T,
                            main = "Oxidation-reduction process"
)

# tiff(here::here("data/Heatmap_OxRed.tif"), unit = "cm", height = 10, width = 15, res = 300)
# OxRed
# dev.off()

#####generate heatmap of NAD associated proteins####
#Used to generate 7D
group_anno <- openxlsx::read.xlsx(here::here("data-raw/NR_intervention/Proteomics/annotation_oxRed.xlsx"), colNames = T)
group_anno <- group_anno %>%
    dplyr::filter(Group== "NAD" | Group == 'NAD/OxPhos' | Group == "NADP")
NAD_proteins <- group_anno$SYMBOL
res_ox <- res %>%
    dplyr::filter(Gene %in% NAD_proteins) %>%
    dplyr::distinct(Gene, .keep_all = T) %>%
    dplyr::arrange(Gene)
rownames(res_ox) <- res_ox$Gene
res_ox <- res_ox %>%
    dplyr::select(-Gene)

setup_ordered <- setup
setup_ordered <- setup_ordered %>%
    dplyr::mutate(
        group = dplyr::case_when(
            group == "WT_Control" ~ "WT Control",
            group == "KO_Control" ~ "HNKO Control",
            group == "WT_NR" ~ "WT NR",
            group == "KO_NR" ~ "HNKO NR"
        )
    )
order <- c("WT Control", "HNKO Control", "WT NR", "HNKO NR")


setup_ordered <- setup_ordered %>%
    dplyr::arrange(Treatment,desc(Genotype))

setup_ordered$sample <- as.character(setup_ordered$sample)
res_ox <- res_ox %>%
    dplyr::select(setup_ordered$sample)


setup_ordered <- setup_ordered %>%
    dplyr::filter(!sample == 330)
key <- as.data.frame(setup_ordered)

key <- key %>%
    dplyr::select(group)
rownames(key) <- setup_ordered$sample
key$group <- factor(key$group, c("WT Control", "HNKO Control", "WT NR", "HNKO NR"))

NADAsso <- pheatmap::pheatmap(res_ox,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              scale = "row",
                              legend = T,
                              na_col = "white",
                              Colv = NA,
                              na.rm = T,
                              cluster_cols = F,
                              cluster_rows = F,
                              fontsize_row = 8,
                              fontsize_col = 8,
                              cellwidth = 10,
                              cellheight = 8,
                              annotation_col = key,
                              show_colnames = F,
                              show_rownames = T,
                              main = "NAD/NADP associated proteins"
)
# tiff(here::here("data/NADAssociatedProtHeatmap.tif"), unit = "cm", height = 15, width = 20, res = 300)
# NADAsso
# dev.off()

#####heatmap for OxPhos classified proteins#####
#USed for figure 8G
group_anno <- openxlsx::read.xlsx(here::here("data-raw/NR_intervention/Proteomics/annotation_oxRed.xlsx"), colNames = T)
group_anno <- group_anno %>%
    dplyr::filter(Group== "OxPhos" | Group == 'NAD/OxPhos')
ox_proteins <- group_anno$SYMBOL
res_ox <- res %>%
    dplyr::filter(Gene %in% ox_proteins) %>%
    dplyr::distinct(Gene, .keep_all = T) %>%
    dplyr::arrange(Gene)
rownames(res_ox) <- res_ox$Gene
res_ox <- res_ox %>%
    dplyr::select(-Gene)

setup_ordered <- setup
setup_ordered <- setup_ordered %>%
    dplyr::mutate(
        group = dplyr::case_when(
            group == "WT_Control" ~ "WT Control",
            group == "KO_Control" ~ "HNKO Control",
            group == "WT_NR" ~ "WT NR",
            group == "KO_NR" ~ "HNKO NR"
        )
    )
order <- c("WT Control", "HNKO Control", "WT NR", "HNKO NR")


setup_ordered <- setup_ordered %>%
    dplyr::arrange(Treatment,desc(Genotype))

class(colnames(res_ox))
setup_ordered$sample <- as.character(setup_ordered$sample)
res_ox <- res_ox %>%
    dplyr::select(setup_ordered$sample)


setup_ordered <- setup_ordered %>%
    dplyr::filter(!sample == 330)
key <- as.data.frame(setup_ordered)

key <- key %>%
    dplyr::select(group)
rownames(key) <- setup_ordered$sample
key$group <- factor(key$group, c("WT Control", "HNKO Control", "WT NR", "HNKO NR"))



NADOxPhos <- pheatmap::pheatmap(res_ox,
                                treeheight_col = 0,
                                treeheight_row = 0,
                                scale = "row",
                                legend = T,
                                na_col = "white",
                                Colv = NA,
                                na.rm = T,
                                cluster_cols = F,
                                cluster_rows = F,
                                fontsize_row = 8,
                                fontsize_col = 8,
                                cellwidth = 10,
                                cellheight = 8,
                                annotation_col = key,
                                show_colnames = F,
                                show_rownames = T,
                                main = "OxPhos Proteins"
)
 # tiff(here::here("data/NADOxPhos.tif"), unit = "cm", height = 15, width = 20, res = 300)
 # NADOxPhos
 # dev.off()

#This was used for supporting figure 4D-E
NAD_genes <- c("Nampt", "Nmnat1", "Nmnat2", "Nmnat3", "Nadsyn1", "Nnt", "Nnmt", "Tdo2", "Afmid", "Kmo", "Qprt", "Kynu", "Ido2", "Adk", "Cd73", "Naprt", "Nmrk1", "Nmrk2", "Slc25a1", "Ent1", "Ent2", "Ent4", "Slc12a8", "Haao", "Aspdh")
NAD_genes <- sort(NAD_genes)

res_NAD <- res %>%
    dplyr::filter(Gene %in% NAD_genes) %>%
    dplyr::distinct(Gene, .keep_all = T)

rownames(res_NAD) <- res_NAD$Gene
res_NAD <- res_NAD %>%
    dplyr::arrange(Gene) %>%
    dplyr::select(-Gene)

setup_ordered <- setup
order <- c("WT Control", "HNKO Control", "WT NR", "HNKO NR")

setup_ordered <- setup_ordered %>%
    dplyr::arrange(Treatment,desc(Genotype))


setup_ordered$sample <- as.character(setup_ordered$sample)
res_NAD <- res_NAD %>%
    dplyr::select(setup_ordered$sample)

key <- as.data.frame(setup_ordered)

key <- key %>%
    dplyr::select(group)
rownames(key) <- setup_ordered$sample

key$group <- factor(key$group, levels = c("WT_Control","KO_Control", "WT_NR", "KO_NR"))

NAD_syn <- pheatmap::pheatmap(res_NAD,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              scale = "row",
                              legend = T,
                              na_col = "white",
                              Colv = NA,
                              na.rm = T,
                              cluster_cols = F,
                              fontsize_row = 9,
                              fontsize_col = 8,
                              cellwidth = 8,
                              cellheight = 10,
                              annotation_col = key,
                              show_colnames = F,
                              show_rownames = T,
                              main = "NAD-synthesis",
                              cluster_rows = F
)

# tiff(here::here("data/Heatmap_NADSYN.tif"), unit = "cm", height = 10, width = 15, res = 300)
# NAD_syn
#
# dev.off()


#Do the same for NAD consumers
NAD_consumers <- c("Sirt1", "Sirt2", "Sirt3", "Sirt4", "Sirt5", "Sirt6", "Sirt7", "Parp1", "Parp4", "Parp14", "Parp3", "Parp12", "Parp12", "Parp9", "Parp10", "Cd38", "Nadk", "Sarm1")
NAD_consumers <- sort(NAD_consumers)
res_con <- res %>%
    dplyr::filter(Gene %in% NAD_consumers) %>%
    dplyr::distinct(Gene, .keep_all = T)

rownames(res_con) <- res_con$Gene
res_con <- res_con %>%
    dplyr::arrange(Gene) %>%
    dplyr::select(-Gene)


res_con <- res_con %>%
    dplyr::select(setup_ordered$sample)


NAD_con <- pheatmap::pheatmap(res_con,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              scale = "row",
                              legend = T,
                              na_col = "white",
                              Colv = NA,
                              na.rm = T,
                              cluster_cols = F,
                              fontsize_row = 9,
                              fontsize_col = 8,
                              cellwidth = 8,
                              cellheight = 10,
                              annotation_col = key,
                              show_colnames = F,
                              show_rownames = T,
                              main = "NAD-consumers",
                              cluster_rows = F
)

# tiff(here::here("data/Heatmap_NADCON.tif"), unit = "cm", height = 10, width = 15, res = 300)
# NAD_con
#
# dev.off()
