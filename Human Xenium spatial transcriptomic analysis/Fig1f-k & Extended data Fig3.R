library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(harmony)
library(ggrastr)


####Xenium data loading and sample assignment####
xenium <- LoadXenium(path, fov = "fov", molecule.coordinates = TRUE, segmentations = "cell")
coords <- GetTissueCoordinates(xenium, which = "centroids")
xenium$sample_id <- case_when(
  coords$x > 5000 & coords$y > 15000 ~ "C-SCLC-1",
  coords$x >= 5200 & coords$y < 15000 & coords$y > 10300 ~ "C-SCLC-2",
  coords$x > 5500 & coords$y < 10300 & coords$y > 5400 ~ "C-SCLC-3",
  coords$x > 6000 & coords$y < 5400 ~ "C-SCLC-4",
  coords$x < 5000 & coords$y > 12500 ~ "C-SCLC-5",
  coords$x < 5200 & coords$y < 12500 & coords$y >= 8400 ~ "C-SCLC-6",
  coords$x < 5500 & coords$y < 8400 & coords$y > 4000 ~ "C-SCLC-7",
  coords$x < 6000 & coords$y < 4000 ~ "T-SCLC-1"
)


# manual corrections for boundary cells
xenium$sample_id[is.na(xenium$sample_id) &
                   coords$x >= 5000 & coords$x <= 5500 &
                   coords$y > 10000 & coords$y < 10500] <- "C-SCLC-6"
xenium$sample_id[coords$x >= 9750 & coords$x <= 10000 &
                   coords$y >= 5000 & coords$y <= 5500] <- "C-SCLC-3"
xenium$sample_id[coords$x >= 3000 & coords$x <= 3500 &
                   coords$y >= 8250 & coords$y <= 8500] <- "C-SCLC-6"
xenium$sample_id[coords$x >= 1000 & coords$x <= 2000 &
                   coords$y >= 8200 & coords$y <= 8500] <- "C-SCLC-6"
xenium$sample_id <- factor(xenium$sample_id)


####Xenium quality control and sketching####
xenium.obj <- subset(xenium, subset = nCount_Xenium > 0 & nFeature_Xenium > 10)
xenium.obj <- NormalizeData(xenium.obj)
xenium.obj <- FindVariableFeatures(xenium.obj)
xenium.obj <- SketchData(xenium.obj, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(xenium.obj) <- "sketch"
xenium.obj <- FindVariableFeatures(xenium.obj)
xenium.obj <- ScaleData(xenium.obj)
xenium.obj <- RunPCA(xenium.obj, npcs = 30)
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30, reduction.name = "umap.sketch.before")


####Xenium sample integration with Harmony####
xenium.obj <- JoinLayers(xenium.obj)
xenium.obj <- RunHarmony(xenium.obj, reduction = "pca",
                         reduction.save = "harmony.sketch", group.by.vars = "sample_id")
xenium.obj <- FindNeighbors(xenium.obj, reduction = "harmony.sketch", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 2)
xenium.obj <- RunUMAP(xenium.obj, reduction = "harmony.sketch", dims = 1:30,
                      reduction.name = "umap.sketch.after", return.model = TRUE)


# Xenium full data projection
xenium.obj <- ProjectIntegration(xenium.obj, sketched.assay = "sketch",
                                 assay = "Xenium", reduction = "harmony.sketch")
xenium.obj <- ProjectData(xenium.obj, assay = "Xenium",
                          full.reduction = "harmony.sketch.full", sketched.assay = "sketch",
                          sketched.reduction = "harmony.sketch", umap.model = "umap.sketch.after",
                          dims = 1:30, refdata = list(cluster_full = "seurat_clusters"))
DefaultAssay(xenium.obj) <- "Xenium"
xenium.obj$cluster_full <- factor(xenium.obj$cluster_full, levels = 0:35)
Idents(xenium.obj) <- "cluster_full"
xenium.obj <- JoinLayers(xenium.obj)


####Xenium cell type annotation####
xenium.obj$celltype <- case_when(
  xenium.obj$cluster_full %in% c(0,1,2,4,5,6,7,9,14,17,18,20,25,28,29,30,31,32,34) ~ "Epithelial",
  xenium.obj$cluster_full %in% c(8,11,19,26,27) ~ "Fibroblast",
  xenium.obj$cluster_full %in% c(3,23,24) ~ "Myeloid",
  xenium.obj %in% c(12,13) ~ "Endothelial",
  xenium.obj$cluster_full == 21 ~ "B",
  xenium.obj$cluster_full == 15 ~ "Plasma",
  .default = "Unknown"
)
xenium.obj$celltype <- factor(xenium.obj$celltype,
                              levels = c("Epithelial", "Myeloid", "T", "B", "Plasma", "Fibroblast", "Endothelial", "Unknown"))


#epithelial subtype
epi.cells <- WhichCells(xenium.obj, expression = celltype == "Epithelial")
expr <- FetchData(xenium.obj,
                  vars = c("INSM1", "UCHL1", "SYP", "MUC1", "NKX2-1"),
                  cells = epi.cells)
is_sclc <- (expr$INSM1 > 0 | expr$UCHL1 > 0 | expr$SYP > 0) & expr$MUC1 == 0
is_luad <- (expr[["NKX2-1"]] > 0 | expr$MUC1 > 0) & expr$INSM1 == 0 & expr$UCHL1 == 0 & expr$SYP == 0
is_hybrid <- (expr$INSM1 > 0 | expr$UCHL1 > 0 | expr$SYP > 0) & expr$MUC1 > 0

xenium.obj$celltype2 <- as.character(xenium.obj$celltype)
xenium.obj$celltype2[epi.cells[is_sclc]] <- "SCLC"
xenium.obj$celltype2[epi.cells[is_luad]] <- "LUAD"
xenium.obj$celltype2[epi.cells[is_hybrid]] <- "SCLC/LUAD"
xenium.obj$celltype2[setdiff(epi.cells, epi.cells[is_sclc | is_luad | is_hybrid])] <- "Other epithelial"
xenium.obj$celltype2 <- factor(xenium.obj$celltype2,
                               levels = c("LUAD", "SCLC/LUAD", "SCLC", "Other epithelial",
                                          "Myeloid", "T", "B", "Plasma", "Fibroblast", "Endothelial", "Unknown"))


####Xenium spatial visualization of cell types####
celltype_cols <- c(
  "Epithelial"  = "#E8593C",
  "Fibroblast"  = "#3B8BD4",
  "Myeloid"     = "#2ECC71",
  "T"           = "#FF4DB8",
  "Endothelial" = "#9B58B6",
  "B"           = "#F39C12",
  "Plasma"      = "#FFD700",
  "Unknown"     = "#808080"
)

p <- ImageDimPlot(xenium.obj, group.by = "celltype", cols = celltype_cols,
                  nmols = 500, size = 0.25, dark.background = FALSE) +
  labs(fill = "Celltype") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    panel.spacing = unit(0, "lines")
  )

p <- rasterise(p, dpi = 300)
pdf("Figure1f_ImageDimPlot_celltype.pdf", height = 10, width = 30)
p
dev.off()


####Xenium UMAP visualization by sample & celltype####
sample.cols <- setNames(brewer.pal(8, "Paired"), levels(xenium.obj$sample_id))
umap_theme <- theme_classic() +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16)
  )
p <- DimPlot(xenium.obj, reduction = "full.umap.sketch.after",
             group.by = "sample_id", cols = sample.cols, alpha = 0.15, raster = TRUE) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "n = 1,099,001", color = "Sample") +
  umap_theme
pdf("FigureS3c_UMAP_sample.pdf", height = 8, width = 10)
p
dev.off()


p <- DimPlot(xenium.obj, reduction = "full.umap.sketch.after",
             group.by = "celltype", cols = celltype_cols, alpha = 0.15, raster = TRUE) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "n = 1,099,001", color = "Celltype") +
  umap_theme
pdf("FigureS3d_UMAP_celltype.pdf", height = 8, width = 10)
p
dev.off()


####Xenium celltype marker dotplot####
markers <- c("EPCAM", "NKX2-1", "MUC1", "INSM1", "UCHL1", "SYP",
             "ASCL1", "NEUROD1", "MYC", "KRAS", "EGFR",
             "CD14", "CD68", "FCGR2A",
             "CD3E", "CD4", "CD8A",
             "MS4A1", "CD19",
             "XBP1", "MZB1",
             "COL5A1", "POSTN",
             "PECAM1", "CD34")
p <- DotPlot(xenium.obj, features = rev(markers), group.by = "celltype",
             scale = TRUE, dot.scale = 8) +
  coord_flip() +
  scale_colour_gradient(low = "white", high = "#08519C") +
  xlab("") + ylab("Cell types") +
  theme(
    panel.border = element_rect(color = "black", size = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
pdf("FigureS3e_xenium.obj_celltype.dotplot.pdf", width = 6, height = 7)
p
dev.off()


####Xenium epithelial subtype UMAP####
p <- DimPlot(epi, reduction = "umap.epi", group.by = "celltype2",
             cols = c("LUAD" = "#E8593C", "SCLC" = "#3B8BD4",
                      "SCLC/LUAD" = "#A855F7", "Other epithelial" = "#A6761D"),
             alpha = 0.15, raster = TRUE) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "", color = "Cell type") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16)
  )
pdf("Figure1g_UMAP_epicell_subtype.pdf", height = 8, width = 10)
p
dev.off()


####Xenium epithelial subtype marker dotplot####
marker_state <- c("EPCAM", "NCAM1", "INSM1", "UCHL1", "SYP", "ASCL1", "NEUROD1",
                  "POU2F3", "YAP1", "NKX2-1", "MUC1", "MYC", "MYCL", "MYCN",
                  "RB1", "KRAS", "EGFR")
Idents(xenium.obj) <- "celltype2"

p <- DotPlot(xenium.obj, features = rev(marker_state), group.by = "celltype2",
             idents = c("LUAD","SCLC/LUAD","SCLC","Other epithelial"), scale = TRUE, dot.scale = 8) + 
  coord_flip() +
  scale_colour_gradient(low = "white", high = "#08519C")+
  xlab("") +
  ylab("") +
  theme(
    panel.border = element_rect(color = "black", size = 0.8)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
pdf("Figure1h_epithelial_subtype_marker_dotplot.pdf", width = 5, height = 8)
p
dev.off()


####Xenium epithelial subtype proportion barplot####
prop_tumor <- xenium.obj@meta.data %>%
  filter(celltype2 %in% c("LUAD", "SCLC/LUAD", "SCLC", "Other epithelial")) %>%
  count(sample_id, celltype2) %>%
  group_by(sample_id) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(celltype2 = factor(celltype2,
                            levels = c("Other epithelial", "LUAD", "SCLC", "SCLC/LUAD")))

pdf("Figure1i_barplot_epithelial_subtype_proportion.pdf", width = 6, height = 4)
ggplot(prop_tumor, aes(x = sample_id, y = pct, fill = celltype2)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c(
    "LUAD" = "#E8A598", "SCLC" = "#7BB5D6",
    "SCLC/LUAD" = "#C4A8D8", "Other epithelial" = "#EDE0B0")) +
  labs(x = "Sample", y = "Percent (%)", fill = "State") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 14)
  )
dev.off()


####Xenium SCLC/LUAD hybrid cell marker expression####
hybrid_cells_subset <- subset(xenium.obj2, subset = celltype2 == "SCLC/LUAD")
gene <- c("MYC","RB1","ASCL1","NEUROD1","POU2F3","YAP1","NKX2-1","MUC1")

pdf("Figure1j_SCLC&LUAD_gene_exp_in_sample_dotplot.pdf", height = 4, width = 6)
DotPlot(hybrid_cells_subset, features = rev(gene), group.by = "sample_id", dot.scale = 8) +
  coord_flip() +
  scale_colour_gradient(low = "white", high = "#08519C") +
  xlab("") + ylab("") +
  theme(
    panel.border = element_rect(color = "black", size = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
dev.off()


####C-SCLC-7 epithelial subtype annotation####
obj_s <- subset(xenium.obj, sample_id == "C-SCLC-7")
DefaultAssay(obj_s) <- "Xenium"
obj_s$celltype3 <- as.character(obj_s$celltype2)
obj_s$celltype3[!obj_s$celltype3 %in% c("SCLC", "SCLC/LUAD", "LUAD", "Other epithelial")] <- "Non-epithelial"
obj_s$celltype3 <- factor(obj_s$celltype3,
                          levels = c("SCLC", "SCLC/LUAD", "LUAD", "Other epithelial", "Non-epithelial"))
celltype3_cols <- c(
  "LUAD"             = "#E8A598",
  "SCLC"             = "#7BB5D6",
  "SCLC/LUAD"        = "#C4A8D8",
  "Other epithelial" = "#EDE0B0",
  "Non-epithelial"   = "#D3D3D3"
)


####C-SCLC-7 spatial marker gene expression####
genes <- c("MYC", "RB1", "ASCL1", "NEUROD1")
p <- ImageFeaturePlot(obj_s, fov = "fov", features = genes,
                      size = 0.6, combine = FALSE, dark.background = FALSE)
p <- lapply(p, function(x) {
  x + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    panel.spacing = unit(0, "lines")
  )
})
p2 <- wrap_plots(p, ncol = 4)
png("FigureS3f_C-SCLC-7_spatial_marker.png", height = 2500, width = 6000, res = 300)
p2
dev.off()


####C-SCLC-7 big field spatial visualization####
cropped.coords <- Crop(obj_s[["fov"]], x = c(1500, 3500), y = c(6000, 8000), coords = "tissue")
obj_s[["bigfield"]] <- cropped.coords
DefaultBoundary(obj_s[["bigfield"]]) <- "segmentation"

p <- ImageDimPlot(obj_s, fov = "bigfield", alpha = 0.8, size = 1,
                  molecules = c("MYC", "RB1", "NEUROD1"),
                  mols.cols = c("MYC" = "black", "RB1" = "#ED2123", "NEUROD1" = "#4068B2"),
                  nmols = 20000, mols.size = 0.02,
                  group.by = "celltype3", cols = celltype3_cols, dark.background = FALSE) +
  labs(fill = "Cell type") +
  scale_x_continuous(breaks = seq(0, 10000, by = 100)) +
  scale_y_continuous(breaks = seq(0, 10000, by = 100)) +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    panel.spacing = unit(0, "lines")
  )

pdf("Figure1k_C-SCLC-7_bigfield_molecule.pdf", height = 18, width = 20)
p
dev.off()


####C-SCLC-7 small field spatial visualization####
cropped.coords <- Crop(obj_s[["fov"]], x = c(1900, 2200), y = c(7100, 7400), coords = "tissue")
obj_s[["field1"]] <- cropped.coords
DefaultBoundary(obj_s[["field1"]]) <- "segmentation"

p <- ImageDimPlot(obj_s, fov = "field1", size = 0.6,
                  group.by = "celltype3", cols = celltype3_cols,
                  border.size = 0.1, dark.background = FALSE) +
  labs(fill = "Cell types") +
  theme(
    panel.border = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    panel.spacing = unit(0, "lines")
  )

pdf("Figure1k_C-SCLC-7_smallfield_imagedimplot.pdf", height = 6, width = 8)
p
dev.off()


p <- ImageDimPlot(obj_s, fov = "field1", alpha = 0.5, size = 1,
                  molecules = c("MYC", "RB1", "NEUROD1"),
                  mols.cols = c("MYC" = "black", "RB1" = "#ED2123", "NEUROD1" = "#4068B2"),
                  nmols = 30000, mols.size = 1,
                  group.by = "celltype3", cols = celltype3_cols, dark.background = FALSE) +
  labs(fill = "Cell types") +
  theme(
    panel.border = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    panel.spacing = unit(0, "lines")
  )

pdf("Figure1k_C-SCLC-7_smallfield_molecule.pdf", height = 6, width = 8)
p
dev.off()
