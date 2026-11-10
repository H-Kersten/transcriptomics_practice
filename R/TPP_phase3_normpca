# ==============================================================================
# 03_normalise_pca.R  --  Phase 2, step 2
#
# Normalise, pick HVGs, PCA, and an UNINTEGRATED UMAP.
#
# The unintegrated embedding is the point of this script. Harmony cannot be
# judged without a baseline: you need to see how much of the structure is
# donor, how much is round, how much is sort, and how much is the round-3
# quality collapse, BEFORE anything is corrected.
#
# Doublets are still retained. They are removed in step 4, after we have looked
# at where they sit.
#
# Bioconductor packages are loaded BEFORE Seurat so Seurat's accessors win the
# masking contest. This is what broke 02_merge.R on the first run.
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(Seurat)
  library(SeuratObject)
})

set.seed(1234)

# ---- config ------------------------------------------------------------------

in_obj    <- "data/processed/merged_raw.rds"
out_obj   <- "data/processed/merged_pca.rds"
fig_dir   <- "results/figures/integration"
tab_dir   <- "results/tables"

n_hvg     <- 2000
n_pcs     <- 50
umap_dims <- 1:30
drop_mt_from_hvg <- TRUE   # keeps the damaged round-3 samples from driving a PC

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()

# ---- load --------------------------------------------------------------------

seu <- readRDS(in_obj)
assay <- SeuratObject::DefaultAssay(seu)
message(sprintf("Loaded %d genes x %d nuclei, assay '%s'",
                nrow(seu), ncol(seu), assay))

# ---- normalise ---------------------------------------------------------------
# LogNormalize, not SCTransform: at 108k nuclei SCT's residual matrix will not
# fit in 16 GB.

seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize",
                             scale.factor = 1e4, verbose = FALSE)

# ---- variable features -------------------------------------------------------

seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst",
                                    nfeatures = n_hvg, verbose = FALSE)

hvg <- SeuratObject::VariableFeatures(seu)

if (drop_mt_from_hvg) {
  mt_hvg <- grep("^MT-", hvg, value = TRUE)
  if (length(mt_hvg)) {
    message("Removing ", length(mt_hvg), " MT- genes from the HVG set: ",
            paste(mt_hvg, collapse = ", "))
    hvg <- setdiff(hvg, mt_hvg)
    SeuratObject::VariableFeatures(seu) <- hvg
  }
}
message("HVGs used: ", length(hvg))

# ---- scale + PCA -------------------------------------------------------------

seu <- Seurat::ScaleData(seu, features = hvg, verbose = FALSE)
seu <- Seurat::RunPCA(seu, features = hvg, npcs = n_pcs, verbose = FALSE)

# scale.data is ~1.7 GB dense and nothing downstream needs it; recompute on a
# gene subset later if you want heatmaps
seu[[assay]]$scale.data <- NULL
gc(verbose = FALSE)

png(file.path(fig_dir, "elbow.png"), width = 1000, height = 700, res = 130)
print(Seurat::ElbowPlot(seu, ndims = n_pcs))
dev.off()

# ---- are the top PCs technical? ----------------------------------------------
# Correlate every PC with the covariates we know are broken. If PC1 tracks
# contamination or depth rather than cell type, Harmony on donor will not save
# you and the round-3 samples have to go.

emb <- SeuratObject::Embeddings(seu, "pca")
md  <- seu@meta.data

covars <- c(paste0("nCount_", assay), paste0("nFeature_", assay),
            "percent_mt", "contamination", "doublet_score")
covars <- intersect(covars, colnames(md))

cor_tab <- sapply(covars, function(v) {
  x <- md[[v]]
  if (grepl("^nCount|^nFeature", v)) x <- log10(x + 1)
  apply(emb, 2, function(pc) cor(pc, x, method = "spearman"))
})
cor_tab <- as.data.frame(cor_tab)
cor_tab$PC <- rownames(cor_tab)
cor_tab <- cor_tab[, c("PC", covars)]

write.csv(cor_tab, file.path(tab_dir, "pc_covariate_cor.csv"), row.names = FALSE)

message("\nSpearman correlation of the first 10 PCs with technical covariates:")
print(head(cor_tab, 10), digits = 2, row.names = FALSE)

# variance explained by round / sort / donor, per PC (eta-squared via one-way F)
grp_tab <- do.call(rbind, lapply(intersect(c("brnum", "round", "sort", "sample_id"),
                                           colnames(md)), function(g) {
  r2 <- apply(emb[, 1:10, drop = FALSE], 2, function(pc) {
    summary(lm(pc ~ md[[g]]))$r.squared
  })
  data.frame(variable = g, PC = names(r2), r_squared = as.numeric(r2))
}))
write.csv(grp_tab, file.path(tab_dir, "pc_group_variance.csv"), row.names = FALSE)

# ---- unintegrated UMAP -------------------------------------------------------

message("\nRunning UMAP on uncorrected PCs (this is the slow part) ...")
seu <- Seurat::RunUMAP(seu, dims = umap_dims, reduction = "pca",
                       reduction.name = "umap_unint",
                       reduction.key  = "umapunint_",
                       verbose = FALSE)

# round-3 and 17c as explicit labels, so they are visible on every plot
seu$qc_group <- "ok"
seu$qc_group[as.character(seu$round) == "3"] <- "round3"
seu$qc_group[seu$sample_id == "17c"]         <- "17c"
seu$qc_group <- factor(seu$qc_group, levels = c("ok", "round3", "17c"))

seu$log10_umi <- log10(md[[paste0("nCount_", assay)]] + 1)

plot_by <- function(g, ...) {
  Seurat::DimPlot(seu, reduction = "umap_unint", group.by = g,
                  raster = TRUE, raster.dpi = c(1200, 1200),
                  pt.size = 0.4, ...) +
    ggtitle(g) + theme(legend.text = element_text(size = 7))
}

p <- (plot_by("qc_group") | plot_by("round")) /
     (plot_by("sort")     | plot_by("brnum")) /
     (plot_by("doublet_class") |
        Seurat::FeaturePlot(seu, reduction = "umap_unint",
                            features = "log10_umi", raster = TRUE,
                            raster.dpi = c(1200, 1200), pt.size = 0.4))

png(file.path(fig_dir, "umap_unintegrated.png"),
    width = 2000, height = 2600, res = 150)
print(p)
dev.off()

png(file.path(fig_dir, "umap_unintegrated_contam.png"),
    width = 1100, height = 900, res = 140)
print(Seurat::FeaturePlot(seu, reduction = "umap_unint",
                          features = "contamination", raster = TRUE,
                          raster.dpi = c(1200, 1200), pt.size = 0.4))
dev.off()

# ---- save --------------------------------------------------------------------

saveRDS(seu, out_obj)
message("\nWrote ", out_obj)
message("Figures in ", fig_dir)
message(sprintf("Elapsed: %.1f min",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

------###----- post-hoc checks
print(seu[["pca"]], dims = 1:4, nfeatures = 15)

emb <- SeuratObject::Embeddings(seu, "pca")
md  <- seu@meta.data
d   <- log10(md$nCount_decontX + 1)

cor(resid(lm(emb[, 1] ~ d)),
    resid(lm(md$contamination ~ d)), method = "spearman")

tapply(emb[, 1], md$sort,     median)
tapply(emb[, 1], md$qc_group, median)

hvg <- SeuratObject::VariableFeatures(seu)
cts <- SeuratObject::LayerData(seu, layer = "counts")
nz  <- tabulate(cts@i + 1L, nbins = nrow(cts))
names(nz) <- rownames(cts)

summary(nz[hvg])
sum(nz[hvg] < 500)
c("SNAP25","SYT1","RBFOX3","PLP1","AQP4","GFAP","CSF1R","MBP") %in% hvg

pc1 <- emb[, 1]
print(quantile(pc1, c(0, .001, .01, .25, .5, .75, .99, .999, 1)), digits = 3)
c(hi = sum(pc1 > 10), lo = sum(pc1 < -10))

ext <- pc1 > quantile(pc1, 0.999)
sort(table(droplevels(md$sample_id[ext])), decreasing = TRUE)
table(md$qc_group[ext])
table(md$doublet_class[ext])

sum(pc1 == stats::median(pc1))
summary(nz[hvg] / ncol(seu))
round(nz[c("SNAP25","SYT1","RBFOX3","PLP1","AQP4","MBP")] / ncol(seu), 3)

n_det <- Matrix::colSums(cts[hvg, , drop = FALSE] > 0)
summary(n_det)
sum(n_det == 0)
tapply(n_det, md$sample_id, median)

frac <- nz / ncol(seu)
keep <- names(frac)[frac >= 0.01]
length(keep)

seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst",
                                    nfeatures = 3000, verbose = FALSE)
hvg2 <- intersect(SeuratObject::VariableFeatures(seu), keep)
hvg2 <- setdiff(hvg2, grep("^MT-", hvg2, value = TRUE))
hvg2 <- head(hvg2, 2000)
SeuratObject::VariableFeatures(seu) <- hvg2

summary(nz[hvg2] / ncol(seu))
c("SNAP25","SYT1","RBFOX3","PLP1","AQP4","GFAP") %in% hvg2


### more pca, somehow different idk
seu <- Seurat::ScaleData(seu, features = hvg2, verbose = FALSE)
seu <- Seurat::RunPCA(seu, features = hvg2, npcs = 50, verbose = FALSE)

emb <- SeuratObject::Embeddings(seu, "pca")
print(seu[["pca"]], dims = 1:3, nfeatures = 12)

tapply(emb[, 1], md$sort, median)
tapply(emb[, 2], md$sort, median)
print(quantile(emb[, 1], c(0, .01, .25, .5, .75, .99, 1)), digits = 3)
