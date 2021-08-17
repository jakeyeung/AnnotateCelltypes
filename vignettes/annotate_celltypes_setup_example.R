# Jake Yeung
# Date of Creation: 2021-08-17
# File: ~/projects/AnnotateCelltypes/vignettes/annotate_celltypes_example.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)

library(Seurat)
library(umap)
library(irlba)

library(AnnotateCelltypes)

# Functions ---------------------------------------------------------------


outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/tmp/cluster_free_celltyping_outputs_withprob.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)

# Load PBMC available from 10x  -------------------------------------------

# pbmc.data <- Read10X(data.dir = "/home/jovyan/data/pbmc3k_filtered_gene_bc_matrices/hg19")
pbmc.data <- Read10X(data.dir = "/home/jyeung/data/public_data/filtered_gene_bc_matrices/hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# Subset fractioon of the data ----------------------------------------------------

jfrac <- 0.8
ncells <- length(colnames(pbmc@assays$RNA@counts))
set.seed(0)
sampled.cells <- sample(x = colnames(pbmc@assays$RNA@counts), size = round(ncells * jfrac), replace = F)
holdout.cells <- colnames(pbmc)[!colnames(pbmc) %in% sampled.cells]


pbmc.sub <- subset(pbmc, cells = sampled.cells)

# Make celltype reference using fraction of the data -----------------------------

pbmc.sub <- NormalizeData(pbmc.sub, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc.sub <- FindVariableFeatures(pbmc.sub, selection.method = "vst", nfeatures = 2000)  # can play with this but for reference doesn't matter

# probably not necessary
all.genes <- rownames(pbmc.sub)
pbmc.sub <- ScaleData(pbmc.sub, features = all.genes)
pbmc.sub <- RunPCA(pbmc.sub, features = VariableFeatures(object = pbmc.sub), npcs = 50)
pbmc.sub <- RunUMAP(pbmc.sub, dims = 1:50)

DimPlot(pbmc.sub, reduction = "pca")
DimPlot(pbmc.sub, reduction = "umap")

# get celltypes by clustering the reference data
pbmc.sub <- FindNeighbors(pbmc.sub, dims = 1:10)
pbmc.sub <- FindClusters(pbmc.sub, resolution = 0.5)

pbmc.markers <- FindAllMarkers(pbmc.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# pbmc.markers %>%
#   group_by(cluster) %>%
#   top_n(n = 2, wt = avg_log2FC)


# FeaturePlot(pbmc.sub, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
FeaturePlot(pbmc.sub, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "TUBB1"))

top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  # top_n(n = 10, wt = avg_log2FC)
  top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc.sub, features = top10$gene) + NoLegend()


new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc.sub)


pbmc.sub <- RenameIdents(pbmc.sub, new.cluster.ids)
DimPlot(pbmc.sub, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# Create reference data  --------------------------------------------------

counts.raw <- pbmc.sub@assays$RNA@counts
genes.keep <- VariableFeatures(pbmc.sub)

counts.raw.filt <- counts.raw[genes.keep, ]

mats.keep.pseudobulk <- lapply(new.cluster.ids, function(clst){
  mat.keep <- subset(x = pbmc.sub, idents = clst)@assays$RNA@counts[genes.keep, ]
  # return pseudobulk
  rowSums(mat.keep)
}) %>%
  bind_cols() %>%
  as.data.frame()

colnames(mats.keep.pseudobulk) <- new.cluster.ids
rownames(mats.keep.pseudobulk) <- genes.keep

#' zero counts are problematic when calculating likelihoods, add pseudocount to compensate
mats.keep.pseudobulk.pcount <- mats.keep.pseudobulk + 1
mats.keep.pseudobulk.pcount.frac <- sweep(mats.keep.pseudobulk.pcount, MARGIN = 2, STATS = colSums(mats.keep.pseudobulk.pcount), FUN = "/")



# Test on one cell --------------------------------------------------------

jcell <- counts.raw.filt[, 1]
jname <- colnames(counts.raw.filt)[1]
jident <- Idents(pbmc.sub)[jname]

#' note the input here is unnormalized counts. No need for log-transform or scaling rows or columns.
#' should be robust to downsampling too.
LL <- CalculateMultinomLL(jcell, mats.keep.pseudobulk.pcount.frac)

dat.LL <- data.frame(ctype = colnames(LL), LL = LL[1, ], stringsAsFactors = FALSE) %>%
  ungroup() %>%
  mutate(prob = softmax(LL))

ggplot(dat.LL, aes(x = ctype, y = LL)) +
  geom_col() +
  theme_bw() +
  ylab("Log likelihood (larger is better)") +
  ggtitle(paste("Probabilistic inference of cell", jname)) +
  xlab("") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.LL, aes(x = ctype, y = prob)) +
  geom_col() +
  theme_bw() +
  ylab("Probability of Celltype") +
  ggtitle(paste("Probabilistic inference of cell", jname)) +
  xlab("") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



print(LL)
P <- softmax(LL)
# convert to probability (sum exponential of logLikelihood as denominator)
print(P)

ctype.pred <- colnames(P)[which.max(P)]
ctype.real <- jident

print(paste("Ctype pred:", ctype.pred, ". Ctype label:", ctype.real))

# Validate the reference can at least recapitulate training cells  --------

LL.all <- apply(counts.raw.filt, 2, function(jcell){
  LL <- CalculateMultinomLL(jcell, mats.keep.pseudobulk.pcount.frac)
})
rownames(LL.all) <- colnames(mats.keep.pseudobulk.pcount.frac)

# ASIDE begin: this is actually an interpretable "latent space"

pca.LL <- prcomp(LL.all, center = TRUE, scale. = TRUE)
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(pca.LL$rotation, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.idents <- data.frame(cell = names(Idents(pbmc.sub)), ctype = Idents(pbmc.sub), stringsAsFactors = FALSE)
dat.umap.long.annot <- left_join(dat.umap.long, dat.idents)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = ctype)) +
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = cbPalette) +
  ggtitle("Cell type predictions from cluster-free model")

# ASIDE end: this is actually an interpretable "latent space"

# let's see how we did

head(LL.all)

LL.max <- apply(LL.all, 2, function(jcol){
  which.max(jcol)
})

ctype.pred.vec <- factor(rownames(LL.all)[LL.max], levels = rownames(LL.all))
ctype.label.vec <- Idents(pbmc.sub)[colnames(LL.all)]

dat.pred <- data.frame(pred = ctype.pred.vec, label = ctype.label.vec, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(pair = interaction(pred, label, sep = "|"),
         is.correct = pred == label)

# count correct/incorrect by celltype
dat.pred.sum <- dat.pred %>%
  group_by(label) %>%
  summarise(ncorrect = length(which(is.correct)),
            nwrong = length(which(!is.correct)),
            accuracy = ncorrect / (ncorrect + nwrong))

#' how well you predict probably depends on how "fine" your cell types are
#' also my "ground truth" come from some clustering analysis, rather than real FACS data
print(dat.pred.sum)


# Test reference on hold-out cells ----------------------------------------

#' Now that we have a reference set, we can have data
#' streaming in one cell at a time and do inference.
#' This means no log transforms, no column or row scalings.

pbmc.holdout <- subset(pbmc, cells = holdout.cells)

counts.raw.ho <- pbmc.holdout@assays$RNA@counts
counts.raw.filt.ho <- counts.raw.ho[genes.keep, ]

#' Loop through the inference machine one cell at a time
LL.all.ho <- apply(counts.raw.filt.ho, 2, function(jcell){
  LL <- CalculateMultinomLL(jcell, mats.keep.pseudobulk.pcount.frac)
})
rownames(LL.all.ho) <- colnames(mats.keep.pseudobulk.pcount.frac)

LL.max.ho <- apply(LL.all.ho, 2, function(jcol){
  which.max(jcol)
})
ctype.pred.vec.ho <- factor(rownames(LL.all.ho)[LL.max.ho], levels = rownames(LL.all.ho))

ctypes.hash <- hash::hash(c(colnames(LL.all.ho), colnames(LL.all)),
                          c(as.character(ctype.pred.vec.ho), as.character(ctype.pred.vec)))

# Do UMAP of full dataset, color by cluster-free celltype predictions --------

#' there are probably ways to do this dim-reduction without normalizing and scaling data...
#'

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)  # can play with this but for reference doesn't matter

# probably not necessary
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc.sub), npcs = 50)
pbmc <- RunUMAP(pbmc, dims = 1:50)

cells.all <- names(Idents(pbmc))
preds.all <- sapply(cells.all, function(x) ctypes.hash[[x]])
Idents(pbmc) <- preds.all

DimPlot(pbmc, reduction = "pca")

#' The cluster-free annotator seems to be doing well
#' The cell type of every cell in this plot was colored one by one, without knowledge of the clustering information
#' I show clustering just to visualize the outputs that these cell-by-cell annotations make sense
DimPlot(pbmc, reduction = "umap") + ggtitle("UMAP of full dataset", "Cell type inferred one at a time using a reference dataset")



# Write prob matrix to ouptput --------------------------------------------

outdir <- "/home/jyeung/projects/AnnotateCelltypes/inst/extdata"
# save output
outrefmat <- file.path(outdir, "prob_mat_PBMC_reference.txt")
outrawmat <- file.path(outdir, "raw_count_PBMC_test_set.txt")
# fwrite(x = mats.keep.pseudobulk.pcount.frac, file = outrefmat)
write.table(x = mats.keep.pseudobulk.pcount.frac, file = outrefmat, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# save raw count matrix
write.table(x = counts.raw.filt.ho, file = outrawmat, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


dev.off()

# rmarkdown::render("~/scripts/ml_playground/cluster_free_celltyping.R")



