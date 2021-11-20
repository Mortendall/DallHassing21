devtools::load_all()
edgeR_data <- list("HNKO 3d" = NA,
                   "HNKO 6d" = NA,
                   "HNKO 12d" = NA,
                   "HNKO 21d" = NA)

for (i in 1:4){
    edgeR_data[[i]]<- openxlsx::read.xlsx(here::here("data/edgeR.xlsx"),i)
}

edgeR_sig <- edgeR_data
for (i in 1:4){
    edgeR_sig[[i]]<-edgeR_sig[[i]] %>%
        dplyr::filter(FDR < 0.05)
    edgeR_sig[[i]]<-edgeR_sig[[i]]$SYMBOL
}
cpm_matrix <- openxlsx::read.xlsx(here("data/CPM_matrix.xlsx"))
metadata <- load_metadata("022_Metadata.xlsx")
metadata <- metadata %>% dplyr::mutate(Group = paste(Genotype, Age, sep = "_"))

#####Upsetplot#####
order_upset <- c("HNKO 21d", "HNKO 12d", "HNKO 6d","HNKO 3d")
upset_data <- data.frame("Gene"= edgeR_sig[[1]],
                         "Group"= names(edgeR_sig[1]))
for (i in 2:length(edgeR_sig)){
    upset_data <- dplyr::add_row(upset_data,"Gene"=edgeR_sig[[i]],
                         "Group"=names(edgeR_sig)[i])
}
upset_data <- dplyr::distinct(upset_data)
upset_data <- upset_data |>
    dplyr::mutate("TRUE"=TRUE)
upset_wide <- tidyr::pivot_wider(upset_data, names_from = Group, values_from = "TRUE",values_fill = FALSE) |>
    dplyr::filter(!is.na(Gene))

rownames(upset_wide)<-upset_wide$Gene

upset_wide <- dplyr::select(upset_wide,-Gene)

upsetRNA <- ComplexUpset::upset(upset_wide, order_upset, name = "", sort_sets = F, themes = ComplexUpset::upset_modify_themes(list(
    'intersections_matrix'=theme(text = element_text(size = 16)))))
#ggplotify to use the object in patchwork
upsetRNA <- ggplotify::as.ggplot(upsetRNA)
upsetRNA <- upsetRNA + ggplot2::ggtitle("Genes with effect of genotype")+ggplot2::theme(plot.title = element_text(size = 18,
                                                                                                                  hjust = 0.5,
                                                                                                                  vjust = 0.95))
#####code for GO analysis####
overlap_genes <- as.data.frame(edgeR_sig[[1]]) %>%
    dplyr::filter(edgeR_sig[[1]] %in% edgeR_sig[[2]] &
                      edgeR_sig[[1]] %in% edgeR_sig[[3]] &
                      edgeR_sig[[1]] %in% edgeR_sig[[4]]
    )
overlap_genes <- overlap_genes %>%
    dplyr::filter(!is.na(overlap_genes))
colnames(overlap_genes) <- "Symbol"

overlap_entrez <- clusterProfiler::bitr(overlap_genes$Symbol,
                                        fromType = "SYMBOL",
                                        toType = "ENTREZID",
                                        OrgDb = "org.Mm.eg.db",
                                        drop = T)


overlap_bg = bitr(edgeR_data[[1]]$SYMBOL,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = "org.Mm.eg.db",
                  drop = T)

goResults_overlap <- enrichGO(gene = overlap_entrez$ENTREZID,
                              universe = overlap_bg$ENTREZID,
                              OrgDb = org.Mm.eg.db,
                              ont = "MF")
rnaGO <- enrichplot::dotplot(goResults_overlap)+ggplot2::ggtitle("Overlap Genes GO-terms")+ggplot2::theme(plot.title = element_text(size = 18,
                                                                                                                                    hjust = 0.5))

# tiff("UpsetRNA.tif", unit = "cm", height = 25, width = 50, res = 600)
 # upsetRNA
 # grid::grid.text("Genes with effect of genotype", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 30))
 # dev.off()

#####Create layout with proteomics figures#####
# load the proteomics plots in. Created in the proteomics TC project
upsetProt <- readRDS(here::here("data/proteinUpset.RDS"))
goProt <- readRDS(here::here("data/protGO.RDS"))

#createlayout for fig. 3
 tiff(here::here("data/figures/upset_GO.tif"), width = 55, height = 50, units = "cm", res = 300)
 (upsetRNA+rnaGO)/(upsetProt+goProt)
 dev.off()

#####Visualization of Epcam, Ncam and Krt7#####
#Supporting Fig. 2
cpm_candidates <- cpm_matrix %>%
    dplyr::filter(SYMBOL %in% c("Ncam1", "Epcam", "Krt7", "Spp1", "Sox9", "Alb" ))
cpm_cand_long <-cpm_candidates %>%
    dplyr::select(-ENSEMBL) %>%
    tidyr::pivot_longer(cols = -SYMBOL, names_to = "ID", values_to = "CPM")
cpm_cand_long$ID <- as.character(cpm_cand_long$ID)
setup <- metadata %>%
    dplyr::mutate(Genotype = case_when(
        Genotype == "KO"~"HNKO",
        Genotype == "WT"~"WT",
        TRUE ~ as.character(Genotype))) %>%
    dplyr::mutate(Group = paste(Genotype, Age, sep = " ")) %>%
    dplyr::select(Basespace.ID, Group, Genotype)
cpm_joined <- dplyr::left_join(cpm_cand_long, setup, by = c("ID"="Basespace.ID"))
cpm_joined$Group <- factor(cpm_joined$Group, levels = c("WT 3d", "WT 6d", "WT 12d", "WT 21d", "HNKO 3d", "HNKO 6d", "HNKO 12d", "HNKO 21d"))

sum_stats <- cpm_joined%>%
    dplyr::select(SYMBOL, CPM, Group, Genotype) %>%
    group_by(SYMBOL, Group) %>%
    rstatix::get_summary_stats(type = "mean_sd")
cpm_joined$Genotype <- factor(cpm_joined$Genotype, c("WT", "HNKO"))
Ncam <- ggplot2::ggplot()+
    geom_point(data = subset(cpm_joined, SYMBOL == "Ncam1"), aes(x = Group, y = CPM, color = Genotype))+
    geom_errorbar(data = subset(sum_stats, SYMBOL == "Ncam1"), aes(x = Group, y = mean, ymin = mean-sd, ymax = mean+sd), stat = "identity")+
    ggtitle("Ncam")+
    xlab("")+
    ylab("log2CPM")+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5,
                                    size = 18),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12))+
              geom_text(aes(label = "**", x = "HNKO 3d", y = 2.7), size = 5)+
              geom_text(aes(label = "**", x = "HNKO 6d", y = 2.7), size = 5)+
              geom_text(aes(label = "*", x = "HNKO 12d", y = 2.7), size = 5)+
              geom_text(aes(label = "**", x = "HNKO 21d", y = 2.7), size = 5)



Epcam <- ggplot2::ggplot()+
    geom_point(data = subset(cpm_joined, SYMBOL == "Epcam"), aes(x = Group, y = CPM, color = Genotype))+
    geom_errorbar(data = subset(sum_stats, SYMBOL == "Epcam"), aes(x = Group, y = mean, ymin = mean-sd, ymax = mean+sd), stat = "identity")+
    ggtitle("Epcam")+
    xlab("")+
    ylab("log2CPM")+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5,
                                    size = 18),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12)
    )+
    geom_text(aes(label = "*", x = "HNKO 3d", y = 3), size = 5)+
    geom_text(aes(label = "**", x = "HNKO 6d", y = 3), size = 5)+
    geom_text(aes(label = "**", x = "HNKO 12d", y = 3.5), size = 5)+
    geom_text(aes(label = "**", x = "HNKO 21d", y = 3), size = 5)

Krt7 <- ggplot2::ggplot()+
    geom_point(data = subset(cpm_joined, SYMBOL == "Krt7"), aes(x = Group, y = CPM, color = Genotype))+
    geom_errorbar(data = subset(sum_stats, SYMBOL == "Krt7"), aes(x = Group, y = mean, ymin = mean-sd, ymax = mean+sd), stat = "identity")+
    ggtitle("Krt7")+
    xlab("")+
    ylab("log2CPM")+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 18),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12)
            )+
    geom_text(aes(label = "*", x = "HNKO 3d", y = 2.7), size = 5)+
    geom_text(aes(label = "**", x = "HNKO 6d", y = 2.5), size = 5)+
    geom_text(aes(label = "*", x = "HNKO 21d", y = 1), size = 5)

Sox9 <- ggplot2::ggplot()+
    geom_point(data = subset(cpm_joined, SYMBOL == "Sox9"), aes(x = Group, y = CPM, color = Genotype))+
    geom_errorbar(data = subset(sum_stats, SYMBOL == "Sox9"), aes(x = Group, y = mean, ymin = mean-sd, ymax = mean+sd), stat = "identity")+
    ggtitle("Sox9")+
    xlab("")+
    ylab("log2CPM")+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5,
                                    size = 18),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12)
    )+
    geom_text(aes(label = "*", x = "HNKO 3d", y = 3), size = 5)+
    geom_text(aes(label = "**", x = "HNKO 6d", y = 3.5), size = 5)
Spp1 <- ggplot2::ggplot()+
    geom_point(data = subset(cpm_joined, SYMBOL == "Spp1"), aes(x = Group, y = CPM, color = Genotype))+
    geom_errorbar(data = subset(sum_stats, SYMBOL == "Spp1"), aes(x = Group, y = mean, ymin = mean-sd, ymax = mean+sd), stat = "identity")+
    ggtitle("Spp1")+
    xlab("")+
    ylab("log2CPM")+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5,
                                    size = 18),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12)
    )+
    geom_text(aes(label = "**", x = "HNKO 3d", y = 9.5), size = 5)+
    geom_text(aes(label = "**", x = "HNKO 6d", y = 9.6), size = 5)+
    geom_text(aes(label = "**", x = "HNKO 12d", y = 9.5), size = 5)+
    geom_text(aes(label = "**", x = "HNKO 21d", y = 9.5), size = 5)

Alb <- ggplot2::ggplot()+
    geom_point(data = subset(cpm_joined, SYMBOL == "Alb"), aes(x = Group, y = CPM, color = Genotype))+
    geom_errorbar(data = subset(sum_stats, SYMBOL == "Alb"), aes(x = Group, y = mean, ymin = mean-sd, ymax = mean+sd), stat = "identity")+
    ggtitle("Alb")+
    xlab("")+
    ylab("log2CPM")+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5,
                                    size = 18),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12)
    )+
    geom_text(aes(label = "**", x = "HNKO 6d", y = 16.4), size = 5)+
    geom_text(aes(label = "p=0.09", x = "HNKO 3d", y = 16.6), size = 5)


CPM_exp <- (Ncam+Epcam+Krt7)/(Spp1+Sox9+Alb)

    # tiff(here::here("data/CPM_HPC.tif"), unit = "cm", height = 20, width = 50, res = 600)
    # CPM_exp
    # dev.off()



