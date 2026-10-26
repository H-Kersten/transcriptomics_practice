# Phase 1 — per-sample QC, ambient correction, doublet removal
# Reads raw 10x matrices, writes one filtered .rds per sample.
#
# Thresholds settled on sample 1c (see README):
#   lower  = 300    top of the ambient plateau; testing the plateau costs
#                   power via the BH correction and runtime, and buys nothing
#   niters = 30000  enough resolution at this test count (Limited/FALSE == 0)
#   mito   = 10%    fixed, not adaptive — an adaptive mito threshold is lenient
#                   on exactly the samples that are most contaminated.
#                   Corroborated: isOutlier() resolved to 10.66% on 1c
#   no upper count filter — the high-count mode is neurons, not artefacts

library(Seurat)
library(Matrix)
library(DropletUtils)
library(celda)
library(scDblFinder)
library(SingleCellExperiment)
library(scuttle)
library(BiocParallel)

set.seed(1234)

# ---- config ------------------------------------------------------------

dev_mode    <- TRUE    # TRUE = first 2 samples only, to check the loop runs
use_workers <- 0       # 0 = serial. 3 gave ~3x on 16 GB with nothing else open

ed_lower  <- 300
ed_niters <- 30000
ed_fdr    <- 0.001
mito_max  <- 10
bg_max    <- 50000     # empty droplets sampled for the decontX background

raw_dir   <- "data/raw/samples"
out_dir   <- "data/processed/filtered"
qc_dir    <- "results/tables/qc"
fig_dir   <- "results/figures/qc"

for (d in c(out_dir, qc_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

BPP <- if (use_workers > 0) SnowParam(workers = use_workers) else SerialParam()

# ---- sample selection --------------------------------------------------
# 20 of 26 libraries: one per donor x sort. Keeps all 10 donors and both
# sort levels. Donor-level n stays 10, which is what the pseudobulk DE
# in Phase 3 depends on — never drop whole donors to save time.

meta <- read.csv("data/processed/sample_metadata.csv", stringsAsFactors = FALSE)

pick <- do.call(rbind, lapply(split(meta, meta$brnum), function(d) {
  do.call(rbind, lapply(split(d, d$sort), function(s) s[1, , drop = FALSE]))
}))

samples <- pick$sample_id
if (dev_mode) samples <- head(samples, 2)

message("selected ", length(samples), " samples")
print(table(pick$sort, pick$sex))

# ---- helpers -----------------------------------------------------------

read_sample <- function(sid, root = raw_dir) {
  m <- Read10X(file.path(root, sid))
  colnames(m) <- paste0(sid, "_", colnames(m))
  m
}

process_sample <- function(sid) {

  out_file <- file.path(out_dir, paste0(sid, ".rds"))
  qc_file  <- file.path(qc_dir,  paste0(sid, "_qc.csv"))
  if (file.exists(out_file)) { message("skip ", sid); return(invisible(NULL)) }

  message("\n=== ", sid, " === ", format(Sys.time(), "%H:%M:%S"))
  t0 <- Sys.time()
  set.seed(1234)

  raw <- read_sample(sid)

  ## 1. empty droplets -----------------------------------------------------
  e <- emptyDrops(raw, lower = ed_lower, niters = ed_niters, BPPARAM = BPP)
  is_cell <- !is.na(e$FDR) & e$FDR <= ed_fdr

  lim <- sum(!is_cell & e$Limited, na.rm = TRUE)
  if (lim > 0) warning(sid, ": ", lim, " barcodes limited by niters — raise it")

  cells <- raw[, is_cell]

  ## 2. QC metrics and adaptive thresholds --------------------------------
  sce  <- SingleCellExperiment(assays = list(counts = cells))
  mito <- grep("^MT-", rownames(sce), value = TRUE)
  if (length(mito) != 13) warning(sid, ": found ", length(mito), " MT- genes, expected 13")

  qc <- perCellQCMetrics(sce, subsets = list(mito = mito))

  lo_counts <- isOutlier(qc$sum,      log = TRUE, type = "lower", nmads = 3)
  lo_genes  <- isOutlier(qc$detected, log = TRUE, type = "lower", nmads = 3)
  hi_mito   <- qc$subsets_mito_percent > mito_max
  discard   <- lo_counts | lo_genes | hi_mito

  # thresholds the MAD rule actually resolved to, for the report
  thr_counts <- attr(lo_counts, "thresholds")[["lower"]]
  thr_genes  <- attr(lo_genes,  "thresholds")[["lower"]]

  sce <- sce[, !discard]

  ## 3. doublets -----------------------------------------------------------
  # per sample, never pooled: artificial doublets must be built from nuclei
  # that could actually have shared a droplet
  sce <- scDblFinder(sce, BPPARAM = SerialParam())

  ## 4. ambient RNA --------------------------------------------------------
  # background taken from the droplets emptyDrops rejected, subsampled —
  # the profile is an average, so 50k estimates it as well as 1.3M
  bg_idx <- which(!is_cell & Matrix::colSums(raw) >= 10)
  if (length(bg_idx) > bg_max) bg_idx <- sample(bg_idx, bg_max)
  sce <- decontX(sce, background = raw[, bg_idx])

  ## 5. assemble and save --------------------------------------------------
  # both matrices kept: raw filtered counts as RNA, corrected as decontX.
  # Which one feeds the Phase 3 pseudobulk is a decision for Phase 3.
  obj <- CreateSeuratObject(counts = counts(sce), project = sid)
  obj[["decontX"]] <- CreateAssay5Object(counts = round(decontXcounts(sce)))

  obj$sample_id     <- sid
  obj$doublet_class <- sce$scDblFinder.class
  obj$doublet_score <- sce$scDblFinder.score
  obj$contamination <- sce$decontX_contamination
  obj$percent_mt    <- qc$subsets_mito_percent[!discard]

  mrow <- meta[meta$sample_id == sid, ][1, ]
  for (v in c("brnum", "round", "sort", "age", "sex", "pmi")) obj[[v]] <- mrow[[v]]

  saveRDS(obj, out_file)

  ## 6. per-sample QC record ----------------------------------------------
  # written per sample, not accumulated in memory, so a crash at sample 17
  # doesn't lose the first 16
  qc_row <- data.frame(
    sample_id      = sid,
    brnum          = mrow$brnum,
    sort           = mrow$sort,
    n_barcodes     = ncol(raw),
    n_emptydrops   = sum(is_cell),
    n_limited      = lim,
    n_low_counts   = sum(lo_counts),
    n_low_genes    = sum(lo_genes),
    n_high_mito    = sum(hi_mito),
    n_discarded    = sum(discard),
    n_after_qc     = ncol(sce),
    n_doublets     = sum(sce$scDblFinder.class == "doublet"),
    thr_counts     = thr_counts,
    thr_genes      = thr_genes,
    median_umi     = median(Matrix::colSums(counts(sce))),
    median_genes   = median(qc$detected[!discard]),
    median_mito    = median(qc$subsets_mito_percent[!discard]),
    median_contam  = median(sce$decontX_contamination),
    minutes        = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  )
  write.csv(qc_row, qc_file, row.names = FALSE)
  print(qc_row)

  ## 7. QC figure ----------------------------------------------------------
  png(file.path(fig_dir, paste0(sid, "_qc.png")), width = 1200, height = 400)
  par(mfrow = c(1, 3))
  hist(log10(qc$sum), breaks = 100, main = paste(sid, "— log10 UMIs"), xlab = "")
  plot(log10(qc$sum), log10(qc$detected), pch = 16, cex = 0.3,
       main = "counts vs genes", xlab = "log10 UMI", ylab = "log10 genes")
  plot(sce$decontX_contamination, obj$percent_mt, pch = 16, cex = 0.3,
       main = "contamination vs mito", xlab = "decontX", ylab = "% mito")
  dev.off()

  rm(raw, cells, sce, obj, e); gc()
  invisible(NULL)
}

# ---- run ---------------------------------------------------------------

for (sid in samples) {
  tryCatch(process_sample(sid),
           error = function(err) message("FAILED ", sid, ": ", conditionMessage(err)))
}

# ---- collect -----------------------------------------------------------

qc_files <- list.files(qc_dir, pattern = "_qc\\.csv$", full.names = TRUE)
if (length(qc_files) > 0) {
  qc_all <- do.call(rbind, lapply(qc_files, read.csv))
  write.csv(qc_all, "results/tables/qc_summary.csv", row.names = FALSE)
  print(qc_all)
}
