# =============================================================================
# Phase 1 — per-sample QC, ambient RNA correction, doublet detection
# =============================================================================
#
# Human hippocampus snRNA-seq (GEO GSE225158-style layout).
# Reads raw 10x matrices, writes one filtered Seurat object per sample.
#
# In   data/raw/samples/<sample_id>/        raw 10x matrix (barcodes/features/mtx)
#      data/processed/sample_metadata.csv   one row per library
#
# Out  data/processed/filtered/<id>.rds     Seurat object, RNA + decontX assays
#      results/tables/qc/<id>_qc.csv        one-row QC record per sample
#      results/tables/qc_summary.csv        all samples combined, with flags
#      results/figures/qc/<id>_qc.png       per-sample diagnostic panel
#      results/figures/qc/overview.png      contamination vs depth, all samples
#
# Usage   set dev_mode <- TRUE for a 2-sample trial run, then FALSE.
#         Restartable: samples with an existing .rds are skipped, and a failure
#         on one sample is logged without stopping the loop.
#
# Runtime ~15-50 min per sample serial, dominated by emptyDrops. 3 Snow workers
#         gave roughly 3x on 16 GB with the browser closed.
#
# -----------------------------------------------------------------------------
# Parameter choices (measured on sample 1c, see README for the full comparison)
# -----------------------------------------------------------------------------
# lower = 300     Top of the ambient plateau. Barcodes at or below this are
#                 pooled to estimate the ambient profile and never tested.
#                 Measured alternatives on 1c:
#                   lower=200, niters=1e5  -> 4881 nuclei, 123 min
#                   lower=250, niters=5e4  -> 4608 nuclei,  49 min
#                   lower=300, niters=3e4  -> 4659 nuclei,  37 min
#                 Testing the plateau inflates the BH correction without
#                 recovering cells, so a lower threshold can yield fewer nuclei.
#
# niters = 30000  Monte Carlo p-values have a floor of 1/(niters+1). At
#                 niters=1e4 every barcode hit that floor and failed BH
#                 correction across 72k tests, so emptyDrops returned only the
#                 barcodes its own `retain` rule kept. Verify per sample that
#                 table(Sig, Limited) has an empty FALSE/TRUE cell.
#
# mito = 10%      Fixed rather than adaptive: an adaptive threshold is most
#                 lenient on the most contaminated samples. Corroborated by
#                 isOutlier(), which resolved to 10.66% on 1c.
#
#                 APPLIED AFTER decontX, on corrected counts. In snRNA-seq the
#                 cytoplasm is stripped during isolation, so mitochondrial
#                 reads are mostly ambient carryover rather than cell damage.
#                 On 25c the filter removed 2870/5046 nuclei when applied to
#                 raw counts and 146/5046 after correction, with median mito
#                 falling from 10.7% to 4.1%.
#
# No upper count filter. The high-count mode in these libraries is neurons
# (SNAP25/SYT1/RBFOX3 enriched 17-40x over the low mode), and the counts-vs-
# genes relation is a single clean curve with no off-curve population. Cutting
# the high tail would delete the neurons. Multiplets are scDblFinder's job.
# =============================================================================

library(Seurat)
library(Matrix)
library(DropletUtils)
library(celda)
library(scDblFinder)
library(SingleCellExperiment)
library(scuttle)
library(BiocParallel)

set.seed(1234)

# ---- config -----------------------------------------------------------------

dev_mode    <- FALSE   # TRUE = first 2 samples only
use_workers <- 3       # 0 = serial

ed_lower  <- 300
ed_niters <- 30000
ed_fdr    <- 0.001     # strict: >70k tests per sample, so 0.05 would admit
                       # thousands of ambient barcodes as a phantom population
mito_max  <- 10
bg_max    <- 50000     # empty droplets sampled for the decontX background

raw_dir <- "data/raw/samples"
out_dir <- "data/processed/filtered"
qc_dir  <- "results/tables/qc"
fig_dir <- "results/figures/qc"

for (d in c(out_dir, qc_dir, fig_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

BPP <- if (use_workers > 0) SnowParam(workers = use_workers) else SerialParam()

# ---- sample selection -------------------------------------------------------
# 20 of 26 libraries: one per donor x sort. Keeps all 10 donors and both sort
# levels. Donor-level n is what the Phase 3 pseudobulk depends on, so libraries
# are dropped from over-represented donors and whole donors never are.

meta <- read.csv("data/processed/sample_metadata.csv", stringsAsFactors = FALSE)

pick <- do.call(rbind, lapply(split(meta, meta$brnum), function(donor) {
  do.call(rbind, lapply(split(donor, donor$sort), function(s) s[1, , drop = FALSE]))
}))

samples <- pick$sample_id
if (dev_mode) samples <- head(samples, 2)

stopifnot(all(samples %in% list.files(raw_dir)))
message("selected ", length(samples), " samples")
print(table(pick$sort, pick$sex))

# ---- helpers ----------------------------------------------------------------

read_sample <- function(sid, root = raw_dir) {
  m <- Read10X(file.path(root, sid))
  colnames(m) <- paste0(sid, "_", colnames(m))
  m
}

process_sample <- function(sid) {

  out_file <- file.path(out_dir, paste0(sid, ".rds"))
  qc_file  <- file.path(qc_dir,  paste0(sid, "_qc.csv"))
  if (file.exists(out_file)) {
    message("skip ", sid)
    return(invisible(NULL))
  }

  message("\n=== ", sid, " === ", format(Sys.time(), "%H:%M:%S"))
  t0 <- Sys.time()
  set.seed(1234)

  raw <- read_sample(sid)

  ## 1. empty droplets --------------------------------------------------------
  # Tests each barcode's expression *profile* against the ambient profile,
  # rather than thresholding on total counts. A small nucleus with few UMIs is
  # retained if its composition is non-ambient.

  e       <- emptyDrops(raw, lower = ed_lower, niters = ed_niters, BPPARAM = BPP)
  is_cell <- !is.na(e$FDR) & e$FDR <= ed_fdr

  lim <- sum(!is_cell & e$Limited, na.rm = TRUE)
  if (lim > 0) {
    warning(sid, ": ", lim, " barcodes limited by niters - raise it")
  }

  cells <- raw[, is_cell]

  ## 2. QC metrics and adaptive thresholds ------------------------------------
  # MAD-based so each sample gets its own cutoff: a shallow library should not
  # be judged against a deep one. sum and detected are not redundant - high
  # counts across few genes means degraded RNA.
  #
  # Note these thresholds typically resolve below the lower=300 floor that
  # emptyDrops already imposes, so they rarely fire here. Retained because they
  # cost nothing and would matter on data entering the pipeline another way.

  sce  <- SingleCellExperiment(assays = list(counts = cells))
  mito <- grep("^MT-", rownames(sce), value = TRUE)
  if (length(mito) != 13) {
    warning(sid, ": found ", length(mito), " MT- genes, expected 13")
  }

  qc <- perCellQCMetrics(sce, subsets = list(mito = mito))

  lo_counts <- isOutlier(qc$sum,      log = TRUE, type = "lower", nmads = 3)
  lo_genes  <- isOutlier(qc$detected, log = TRUE, type = "lower", nmads = 3)
  discard   <- lo_counts | lo_genes

  thr_counts <- attr(lo_counts, "thresholds")[["lower"]]
  thr_genes  <- attr(lo_genes,  "thresholds")[["lower"]]

  sce     <- sce[, !discard]
  qc_keep <- qc[!discard, ]   # kept in step with sce through every subsetting

  ## 3. doublets --------------------------------------------------------------
  # Per sample, never pooled: artificial doublets must be built from nuclei
  # that could actually have shared a droplet. Run on integer counts, since
  # the simulation sums pairs of real cells.
  # Flagged, not removed - doublets should sit between clusters on the Phase 2
  # UMAP, and if they instead form their own blob the caller has found
  # something else.

  sce <- scDblFinder(sce, BPPARAM = SerialParam())

  ## 4. ambient RNA -----------------------------------------------------------
  # Background taken from the droplets emptyDrops rejected rather than inferred
  # from the cells. Subsampled because the profile is an average over barcodes,
  # so 50k estimates it as well as 1.3M at a fraction of the memory.

  bg_idx <- which(!is_cell & Matrix::colSums(raw) >= 10)
  if (length(bg_idx) > bg_max) bg_idx <- sample(bg_idx, bg_max)
  bg <- SingleCellExperiment(assays = list(counts = raw[, bg_idx]))

  sce <- decontX(sce, background = bg)

  ## 4b. mito filter, on corrected counts -------------------------------------
  # Whatever mitochondrial signal survives decontX is genuine damage rather
  # than soup. See the header for the measured effect of the ordering.

  dec      <- decontXcounts(sce)
  mito_dec <- 100 * Matrix::colSums(dec[mito, , drop = FALSE]) / Matrix::colSums(dec)
  hi_mito  <- mito_dec > mito_max

  message(sid, ": mito filter removes ", sum(hi_mito), " of ", length(hi_mito),
          " (median contamination ", round(median(sce$decontX_contamination), 3), ")")

  sce      <- sce[, !hi_mito]
  qc_keep  <- qc_keep[!hi_mito, ]
  mito_dec <- mito_dec[!hi_mito]

  ## 5. assemble and save -----------------------------------------------------
  # Both matrices kept: raw filtered counts as RNA, corrected as decontX.
  # Which feeds the Phase 3 pseudobulk is a decision deferred to Phase 3.

  obj <- CreateSeuratObject(counts = counts(sce), project = sid)
  obj[["decontX"]] <- CreateAssay5Object(counts = round(decontXcounts(sce)))

  obj$sample_id     <- sid
  obj$doublet_class <- sce$scDblFinder.class
  obj$doublet_score <- sce$scDblFinder.score
  obj$contamination <- sce$decontX_contamination
  obj$percent_mt    <- mito_dec

  # TODO: add the condition/diagnosis column once confirmed in the metadata -
  # Phase 3 needs it in the design matrix.
  mrow <- meta[meta$sample_id == sid, ]
  if (nrow(mrow) != 1) {
    stop(sid, ": metadata join matched ", nrow(mrow), " rows")
  }
  mrow <- mrow[1, ]
  for (v in c("brnum", "round", "sort", "age", "sex", "pmi")) obj[[v]] <- mrow[[v]]

  saveRDS(obj, out_file)

  ## 6. per-sample QC record --------------------------------------------------
  # Written per sample rather than accumulated in memory, so a crash at sample
  # 17 does not lose the first 16.

  qc_row <- data.frame(
    sample_id       = sid,
    brnum           = mrow$brnum,
    sort            = mrow$sort,
    n_barcodes      = ncol(raw),
    n_emptydrops    = sum(is_cell),
    n_limited       = lim,
    n_low_counts    = sum(lo_counts),
    n_low_genes     = sum(lo_genes),
    n_high_mito     = sum(hi_mito),
    n_discarded     = sum(discard) + sum(hi_mito),
    n_after_qc      = ncol(sce),
    n_doublets      = sum(sce$scDblFinder.class == "doublet"),
    thr_counts      = thr_counts,
    thr_genes       = thr_genes,
    median_umi      = median(Matrix::colSums(counts(sce))),
    median_genes    = median(qc_keep$detected),
    median_mito_raw = median(qc_keep$subsets_mito_percent),
    median_mito_dec = median(mito_dec),
    median_contam   = median(sce$decontX_contamination),
    minutes         = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  )
  write.csv(qc_row, qc_file, row.names = FALSE)
  print(qc_row)

  ## 7. per-sample QC figure --------------------------------------------------

  png(file.path(fig_dir, paste0(sid, "_qc.png")), width = 1200, height = 400)
  par(mfrow = c(1, 3))
  hist(log10(qc_keep$sum), breaks = 100,
       main = paste(sid, "- log10 UMIs"), xlab = "")
  plot(log10(qc_keep$sum), log10(qc_keep$detected), pch = 16, cex = 0.3,
       main = "counts vs genes", xlab = "log10 UMI", ylab = "log10 genes")
  plot(sce$decontX_contamination, qc_keep$subsets_mito_percent, pch = 16, cex = 0.3,
       main = "contamination vs raw mito", xlab = "decontX", ylab = "% mito (raw)")
  dev.off()

  rm(raw, cells, sce, obj, e, dec, bg)
  gc()
  invisible(NULL)
}

# ---- run --------------------------------------------------------------------

for (sid in samples) {
  tryCatch(
    process_sample(sid),
    error = function(err) message("FAILED ", sid, ": ", conditionMessage(err))
  )
}

# ---- collect ----------------------------------------------------------------

qc_files <- list.files(qc_dir, pattern = "_qc\\.csv$", full.names = TRUE)

if (length(qc_files) > 0) {

  qc_all <- do.call(rbind, lapply(qc_files, read.csv))

  # Flags for Phase 2, not exclusions. Whether a flagged sample is dropped is
  # decided on the integrated UMAP: a sample that mixes normally with the rest
  # is kept regardless of its contamination estimate.
  #
  # Contamination is partly mechanical - a droplet picks up a roughly fixed
  # amount of ambient RNA, so it is a larger fraction of a small glial nucleus
  # than of a 25k-UMI neuron. Read median_contam against median_umi, not alone.
  qc_all$flag_contam <- qc_all$median_contam > 0.40
  qc_all$flag_sort   <- qc_all$sort == "PI+NeuN+" & qc_all$median_umi < 6000

  write.csv(qc_all, "results/tables/qc_summary.csv", row.names = FALSE)
  print(qc_all)

  message("total nuclei retained: ", sum(qc_all$n_after_qc))
  message("contamination vs depth, r = ",
          round(cor(log10(qc_all$median_umi), qc_all$median_contam), 3))

  png(file.path(fig_dir, "overview.png"), width = 800, height = 600)
  plot(log10(qc_all$median_umi), qc_all$median_contam, pch = 16,
       xlab = "log10 median UMI", ylab = "median decontX contamination",
       main = "ambient contamination vs library depth")
  text(log10(qc_all$median_umi), qc_all$median_contam, qc_all$sample_id,
       pos = 4, cex = 0.7)
  dev.off()
}

writeLines(capture.output(sessionInfo()), "results/tables/phase1_sessionInfo.txt")
