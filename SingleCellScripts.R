
#The data pre-processing was performed by Oliver Knights Møller and Jonatan Thompson from the CBMR single cell transcriptomics platform and is reproduced here for transparancy

SCOP_2021_0132.210202_liver <- readRDS(here::here("data-raw/NR_intervention/SingleNucleusRNAseq/SCOP_2021_0132.210202_liver.raw.rds"))
### Visualizing HTO counts and filtering
VlnPlot(SCOP_2021_0132.210202_liver, features = c('nCount_HTO'), log = T)
SCOP_2021_0132.210202_liver <- subset(SCOP_2021_0132.210202_liver, subset = nCount_HTO > 100 & nCount_HTO < 1250)
VlnPlot(SCOP_2021_0132.210202_liver, features = c('nCount_HTO'), log = T)

### HTO normalization and demultiplexing
SCOP_2021_0132.210202_liver <- NormalizeData(SCOP_2021_0132.210202_liver, assay = 'HTO', normalization.method = 'CLR', verbose = F)
SCOP_2021_0132.210202_liver <- HTODemux(SCOP_2021_0132.210202_liver, assay = 'HTO', positive.quantile = 0.99, verbose = F)
table(SCOP_2021_0132.210202_liver$HTO_classification.global)

### Remove negatives and proccess and visualize samples
SCOP_2021_0132.210202_liver <- subset(SCOP_2021_0132.210202_liver, idents = 'Negative', invert = TRUE)
# Calculate a HTO UMAP embedding
DefaultAssay(SCOP_2021_0132.210202_liver) <- 'HTO'
SCOP_2021_0132.210202_liver <- ScaleData(SCOP_2021_0132.210202_liver,
                                         features = rownames(SCOP_2021_0132.210202_liver),
                                         verbose = F)
SCOP_2021_0132.210202_liver <- RunPCA(SCOP_2021_0132.210202_liver,
                                      features = rownames(SCOP_2021_0132.210202_liver),
                                      reduction.name = 'hto.pca',
                                      approx = F,
                                      verbose = F)
SCOP_2021_0132.210202_liver <- RunUMAP(SCOP_2021_0132.210202_liver,
                                       reduction = 'hto.pca',
                                       dims = 1:max(SCOP_2021_0132.210202_liver$nFeature_HTO),
                                       reduction.name = 'hto.umap',
                                       reduction.key = 'HUMAP_',
                                       verbose = F)



DimPlot(SCOP_2021_0132.210202_liver)

### Subset singlets and visualize RNA features
# Subset singlets
Idents(SCOP_2021_0132.210202_liver) <- 'HTO_classification.global'
SCOP_2021_0132.210202_liver <- subset(SCOP_2021_0132.210202_liver, idents = 'Singlet')

# Now we filter out outliers based on feature count
VlnPlot(SCOP_2021_0132.210202_liver, features = c('nFeature_RNA', 'nCount_RNA'), ncol = 2)
SCOP_2021_0132.210202_liver <- subset(SCOP_2021_0132.210202_liver, subset = nFeature_RNA > 0 & nFeature_RNA < 5500)

VlnPlot(SCOP_2021_0132.210202_liver, features = c('nFeature_RNA', 'nCount_RNA'), ncol = 2)

### Proccess and generate elbowplot to inspect sufficient number of dimensions to use
DefaultAssay(SCOP_2021_0132.210202_liver) <- 'RNA'
#MD: Jonathan Thompson recommended SCtransform rather than NormalizeData. Try that one instead.
SCOP_2021_0132.210202_liver <- NormalizeData(SCOP_2021_0132.210202_liver, verbose = F)
SCOP_2021_0132.210202_liver <- ScaleData(SCOP_2021_0132.210202_liver, verbose = F)
SCOP_2021_0132.210202_liver <- FindVariableFeatures(SCOP_2021_0132.210202_liver, verbose = T)
SCOP_2021_0132.210202_liver <- RunPCA(SCOP_2021_0132.210202_liver, verbose = F)

ElbowPlot(SCOP_2021_0132.210202_liver, ndims = 50)

### Run UMAP and visualize
SCOP_2021_0132.210202_liver <- RunUMAP(SCOP_2021_0132.210202_liver,
                                       dims = 1:10,
                                       verbose = F)
Idents(SCOP_2021_0132.210202_liver) <- 'HTO_classification'

DimPlot(SCOP_2021_0132.210202_liver)

#####End of single cell omics platform analysis#####
## From here, all subsequent analysis has been performed by the first authors, while the single cell platform has provided protocols and feedhack
#Run UMAP and clustering for singlets
VlnPlot(object = SCOP_2021_0132.210202_liver,
        features =c("nCount_RNA",
                    "nFeature_RNA"),
        pt.size = 0.1,
        group.by = "hash.ID")
#scaling and transformation of data using glmGamPoi method
SCOP <- SCTransform(SCOP_2021_0132.210202_liver)
npcs <- 35 # max n PCs to check
randomSeed <- 12345
set.seed(randomSeed)
#runPCA
SCOP <- RunPCA(SCOP,
               npcs=npcs,
               seed.use=randomSeed,
               verbose=F)
#runUMAP
SCOP <- RunUMAP(SCOP,
                dims = 1:npcs,
                seed.use = randomSeed+1)
SCOP <- FindNeighbors(SCOP,
                      reduction = "pca",
                      dims = 1:npcs)

SCOP <- FindClusters(SCOP,
                     random.seed = randomSeed,
                     resolution = 1.2,
                     verbose = FALSE)

#saveRDS(SCOP, here::here("data/single_cell/SCOP_SCtransform.rds"))
######Find markers for the clusters####
SCOP_markers <- FindAllMarkers(SCOP, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
#openxlsx::write.xlsx(SCOP_markers, here::here("data/single_cell/markers.xlsx"))
SCOP_markers <- SCOP_markers %>%
    dplyr::filter(p_val_adj<0.05) %>%
    dplyr::group_by(cluster)
SCOP_markers %>% dplyr::top_n(n = 5, wt = avg_log2FC) %>% View()

DimPlot(SCOP,
        group.by="hash.ID",
        repel = T,
        pt.size = 0.1, label.size = 5,
        label = TRUE) + NoLegend()

#add genotype to metadata
Idents(SCOP)<-"HTO_classification"
SCOP <- RenameIdents(SCOP,
                     "322-KO"="HNKO",
                     "342-WT"="WT",
                     "318-WT"="WT",
                     "328-KO"="HNKO",
                     "338-WT"="WT",
                     "344-KO"="HNKO")
SCOP$genotype <- Idents(SCOP)

######The following code was used to create supporting figure 6A####
#Plot shows a nicer classification by individual
Idents(SCOP)<-"HTO_classification"
SCOP <- RenameIdents(SCOP,
                     "322-KO"="HNKO 1",
                     "342-WT"="WT 1",
                     "318-WT"="WT 2",
                     "328-KO"="HNKO 3",
                     "338-WT"="WT 3",
                     "344-KO"="HNKO 2")
SCOP$new_name <- Idents(SCOP)
Dimplot_cell <- DimPlot(
    SCOP,
    repel = T,
    pt.size = 0.3, label.size = 5,
    label = TRUE) +
    NoLegend() +
    ggplot2::ggtitle("Liver cell populations grouped by individual") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 18),
                   strip.text.x = ggplot2::element_text(size = 18))


# tiff(here::here("Data/single_cell/Dimplot_individual.tif"), units = "cm", width = 30, height = 20, res = 300)
# Dimplot_cell
# dev.off()

#####Add celltypes to metadata####
cell_types <- c("HNKO 1",
                "Hepatocytes 1",
                "Periportal 1",
                "HNKO 2",
                "Hepatocytes 2",
                "Hepatocytes 3",
                "HNKO 3",
                "Periportal 2",
                "Endothelial",
                "Central 1",
                "Periportal 3",
                "Hepatocytes 4",
                "Stellate",
                "CD45+ 1",
                "Central HNKO 1",
                "Kupfer",
                "Central HNKO 2",
                "Cholangiocyte 1",
                "CD45+ 2",
                "Cholangiocyte 2",
                "CD45+ 3")

Idents(SCOP) <- "seurat_clusters"

names(cell_types)<-levels(SCOP)
SCOP <- RenameIdents(SCOP, cell_types)
SCOP$cell_types <- Idents(SCOP)
#factor order of levels
SCOP$cell_types<- factor(SCOP$cell_types, levels = c(
    "Hepatocytes 1",
    "Hepatocytes 2",
    "Hepatocytes 3",
    "Hepatocytes 4",
    "Central 1",
    "Periportal 1",
    "Periportal 2",
    "Periportal 3",
    "HNKO 1",
    "HNKO 2",
    "HNKO 3",
    "Central HNKO 1",
    "Central HNKO 2",
    "Endothelial",
    "Cholangiocyte 1",
    "Cholangiocyte 2",
    "Stellate",
    "Kupfer",
    "CD45+ 1",
    "CD45+ 2",
    "CD45+ 3"))

plot_markers <- c("Acaa1b", "Scd1","Rbp4", "Slc1a2" , "Glul","Cyp2f2","Cyp17a1", "Stab2", "Ptprb", "Ctnnd2","Spp1", "Dcn", "Rbms3",  "Hdac9","Inpp4b", "Timd4", "Ptprc")
Seurat::Idents(SCOP)<-"cell_types"
#####The following code was used for supporting fig 6B####
#Figure shows dotplot to show validity of markers
Dotplot <- Seurat::DotPlot(SCOP,  features = plot_markers)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14,angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 14),
                   axis.title.x = ggplot2::element_text(size = 0),
                   axis.title.y = ggplot2::element_text(size = 0))
# tiff(here::here("Data/single_cell/Dotplot_all_cells.tif"), units = "cm", width = 30, height = 20, res = 300)
# Dotplot
# dev.off()

#####Dimplot of liver cell populations split by genotype####
#The following code was made to generate figure 10A
Idents(SCOP)<-"cell_types"

SCOP$genotype <- factor(SCOP$genotype, c("WT", "HNKO"))
Dimplot_cell <- DimPlot(
    SCOP,
    split.by = "genotype",
    repel = T,
    pt.size = 0.1, label.size = 5,
    label = TRUE) +
    NoLegend() +
    ggplot2::ggtitle("Liver cell populations") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 18),
                   strip.text.x = ggplot2::element_text(size = 18))
# tiff(here::here("Data/single_cell/Dimplot.tif"), units = "cm", width = 30, height = 20, res = 300)
# Dimplot_cell
# dev.off()

#####subset hepatocyte populations####
hepatocytes <- c("0","1", "2", "3", "4", "5", "6", "7", "9", "11", "14", "16")
SCOP_hepatocytes <- subset(SCOP, subset = seurat_clusters %in% hepatocytes)
SCOP_hepatocytes <- RunPCA(SCOP_hepatocytes,
                           npcs=npcs,
                           seed.use=randomSeed,
                           verbose=F)
#runUMAP
SCOP_hepatocytes  <- RunUMAP(SCOP_hepatocytes,
                             dims = 1:npcs,
                             seed.use = randomSeed+1)
SCOP_hepatocytes  <- FindNeighbors(SCOP_hepatocytes,
                                   reduction = "pca",
                                   dims = 1:npcs)

SCOP_hepatocytes  <- FindClusters(SCOP_hepatocytes ,
                                  random.seed = randomSeed,
                                  resolution = 1.2,
                                  verbose = FALSE)
hepatocytes_features <- FindAllMarkers(SCOP_hepatocytes, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
hepatocytes_features <- hepatocytes_features  %>%
    dplyr::filter(p_val_adj<0.05) %>%
    dplyr::group_by(cluster)
hepatocytes_features %>% dplyr::top_n(n = 5, wt = avg_log2FC) %>% View()
Idents(SCOP_hepatocytes)<-"HTO_classification"
SCOP_hepatocytes <- RenameIdents(SCOP_hepatocytes,
                                 "322-KO"="HNKO",
                                 "342-WT"="WT",
                                 "318-WT"="WT",
                                 "328-KO"="HNKO",
                                 "338-WT"="WT",
                                 "344-KO"="HNKO")
SCOP_hepatocytes$genotype <- Idents(SCOP_hepatocytes)
Idents(SCOP_hepatocytes)<-"HTO_classification"
SCOP_hepatocytes <- RenameIdents(SCOP_hepatocytes,
                                 "322-KO"="HNKO 1",
                                 "342-WT"="WT 1",
                                 "318-WT"="WT 2",
                                 "328-KO"="HNKO 3",
                                 "338-WT"="WT 3",
                                 "344-KO"="HNKO 2")
SCOP_hepatocytes$new_name <- Idents(SCOP_hepatocytes)
#####Dimplot by donor mouse####
#The following code was used to generate supporting fig 7A

Hep_by_ind <- DimPlot(
    SCOP_hepatocytes,
    repel = T,
    pt.size = 0.1, label.size = 5,
    label = TRUE) +
    NoLegend() +
    ggplot2::ggtitle("Hepatocyte Populations grouped by individual") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 18),
                   strip.text.x = ggplot2::element_text(size = 18))

  # tiff(here::here("Data/single_cell/Hep_by_ind.tif"), units = "cm", width = 30, height = 20, res = 300)
  # Hep_by_ind
  # dev.off()

  #####Clustering of hepatocytes by marker genes and necrosis scores####

#####annotate based on top non-pseudogene marker####

marker_genes <- c("Sds_1",
                  "Cyp3a41a",
                  "Ntrk2",
                  "Sord",
                  "Rtn4",
                  "Zbtb16",
                  "ND4",
                  "Sds_2",
                  "Slco1a1",
                  "Slc1a2_1",
                  "Slc1a2_2",
                  "Chrm3",
                  "Slc1a2_3")
Idents(SCOP_hepatocytes)<-"seurat_clusters"
names(marker_genes)<-levels(SCOP_hepatocytes)
SCOP_hepatocytes <- RenameIdents(SCOP_hepatocytes, marker_genes)
SCOP_hepatocytes$marker_genes<- Idents(SCOP_hepatocytes)
Idents(SCOP_hepatocytes)<-"marker_genes"
SCOP_hepatocytes$genotype <- factor(SCOP_hepatocytes$genotype, c("WT", "HNKO"))
#####hepatocyt plot colored by fibrosis scores####
Idents(SCOP_hepatocytes)<-"HTO_classification"
SCOP_hepatocytes <- RenameIdents(SCOP_hepatocytes,
                                 "322-KO"="Necrosis Score 2",
                                 "342-WT"="Necrosis Score 0",
                                 "318-WT"="Necrosis Score 0",
                                 "328-KO"="Necrosis Score 0",
                                 "338-WT"="Necrosis Score 0",
                                 "344-KO"="Necrosis Score 2")
SCOP_hepatocytes$necrosis <- Idents(SCOP_hepatocytes)
SCOP_hepatocytes$necrosis <- factor(SCOP_hepatocytes$necrosis, c("Necrosis Score 0", "Necrosis Score 2"))
SCOP_hepatocytes$genotype <- factor(SCOP_hepatocytes$genotype, c("WT", "HNKO"))

##The following code was used to generate figure 10B+C
plot1 <- DimPlot(
    SCOP_hepatocytes,
    split.by = "genotype",
    group.by = "necrosis",
    pt.size = 0.1
) + ggplot2::ggtitle("Hepatocytes split by Genotype") + ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
    strip.text.x = ggplot2::element_text(size = 12)
) #+ NoLegend()

plot2 <- DimPlot(
    SCOP_hepatocytes,
    split.by = "genotype",
    group.by = "marker_genes",
    label = T,
    repel = T,
    pt.size = 0.05,
    label.size = 4
) + NoLegend() + ggplot2::ggtitle("Cluster Genes for Hepatocyte Subpopulations") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
                   strip.text.x = ggplot2::element_text(size = 12))
plot1
plot2
plot_grid <- gridExtra::grid.arrange(plot2, plot1, ncol = 2)
plot_grid <- ggplotify::as.ggplot(plot_grid)

# tiff(here::here("Data/single_cell/hepatocyte_sub.tif"), units = "cm", width = 40, height = 20, res = 300)
# plot_grid
# dev.off()
#####Save hepatocyte analysis and markers####
# saveRDS(SCOP_hepatocytes, here::here("data/single_cell/scop_hepatocytes.rds"))
# openxlsx::write.xlsx(hepatocytes_features, here::here("data/single_cell/hepatocyte_markers.xlsx"))

#####Analysis of the ND4 cluster####
ND4_cluster <- subset(SCOP_hepatocytes@meta.data, marker_genes == "ND4")
FeaturePlot(SCOP_hepatocytes, "ND4", split.by = "genotype")
ND_4_genes <- hepatocytes_features %>% dplyr::filter(cluster ==6)
ND_4_genes <- slice_head(ND_4_genes, n = 20)
ND_4_genes <- ND_4_genes %>% dplyr::arrange(gene)
DotPlot(SCOP_hepatocytes, group.by = "new_name", features = ND_4_genes$gene)
SCOP_hepatocytes$new_name <- factor(SCOP_hepatocytes$new_name, c("HNKO 1", "HNKO 2", "HNKO 3", "WT 1", "WT 2", "WT 3"))
Idents(SCOP_hepatocytes)<-"new_name"
#This code was used to generate supporting fig 7C
ND4_dot <- DotPlot(SCOP_hepatocytes, group.by = "new_name", features = ND_4_genes$gene)+
    ggplot2::ggtitle("Top20 genes from the ND4-cluster")+
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(
            size = 11,
            angle = 45,
            hjust = 1
        ),
        axis.text.y = ggplot2::element_text(size = 14),
        axis.title.x = ggplot2::element_text(size = 0),
        axis.title.y = ggplot2::element_text(size = 0),
        plot.title = ggplot2::element_text(hjust = 0.5)
    )

# tiff(here::here("data/single_cell/ND4dot.tif"), units = "cm", width = 20, height = 15, res = 300)
# ND4_dot
# dev.off()

#####Run GO analysis on ND4 cluster####
ND_4 <- hepatocytes_features %>% dplyr::filter(cluster ==6)

universe <- rownames(SCOP_hepatocytes)
bg = clusterProfiler::bitr(universe,
                           fromType = "SYMBOL",
                           toType = "ENTREZID",
                           OrgDb = "org.Mm.eg.db",
                           drop = T)
eg<- clusterProfiler::bitr(ND_4$gene,
                           fromType = "SYMBOL",
                           toType = "ENTREZID",
                           OrgDb = "org.Mm.eg.db",
                           drop = T)

goResults_ND_4 <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                            universe = bg$ENTREZID,
                                            OrgDb = org.Mm.eg.db,
                                            ont = "MF")
#The following code was used to generate supporting figure 10C
ND4_GO <- clusterProfiler::dotplot(goResults_ND_4)+
    ggplot2::theme(title = element_text(size = 16))+
    ggplot2::ggtitle("Markers from the ND4 cluster")

# tiff(here::here("data/single_cell/ND4GO_MF.tif"), units = "cm", width = 20, height = 15, res = 300)
# ND4_GO
# dev.off()
