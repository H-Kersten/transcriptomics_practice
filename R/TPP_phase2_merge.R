# ==============================================================================
# 02_merge.R  --  Phase 2, step 1
#
# Merge the 20 per-sample Seurat objects from Phase 1 into one Seurat v5 object.
#
#   - decontX-corrected counts only (raw RNA counts stay on disk per sample)
#   - doublets RETAINED and labelled; they are filtered in step 3, after we have
#     checked on the UMAP that they sit between clusters
#   - genes detected in fewer than `min_cells_per_gene` nuclei are dropped
#   - built as a single `counts` layer, NOT via merge(), so memory stays
#     predictable and v5 does not split into 20 layers
#
# Peak memory is around 6 GB. Run this in a fresh R session with nothing else
# loaded.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# ---- config ------------------------------------------------------------------

in_dir             <- "data/processed/filtered"
out_obj            <- "data/processed/merged_raw.rds"
out_tab            <- "results/tables/merge_summary.csv"
assay_in           <- "decontX"   # assay to pull counts from
assay_out          <- "decontX"   # assay name in the merged object
min_cells_per_gene <- 3

# ---- setup -------------------------------------------------------------------

stopifnot(dir.exists(in_dir))
dir.create(dirname(out_obj), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_tab), recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()

files <- sort(list.files(in_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0L) stop("No .rds files found in ", in_dir)
message(sprintf("Found %d sample objects in %s", length(files), in_dir))

# counts matrix from either a v3 Assay or a v5 Assay5
get_counts <- function(obj, assay) {
  a <- obj[[assay]]
  if (inherits(a, "Assay5")) {
    SeuratObject::LayerData(a, layer = "counts")
  } else {
    SeuratObject::GetAssayData(a, slot = "counts")
  }
}

# ---- load samples ------------------------------------------------------------

mats      <- vector("list", length(files))
mds       <- vector("list", length(files))
ref_genes <- NULL

for (i in seq_along(files)) {
  f  <- files[i]
  id <- sub("\\.rds$", "", basename(f))

  obj <- readRDS(f)

  if (!assay_in %in% Assays(obj)) {
    stop(sprintf("%s: no '%s' assay (found: %s)",
                 id, assay_in, paste(Assays(obj), collapse = ", ")))
  }

  m <- get_counts(obj, assay_in)

  # enforce one shared gene universe, in identical row order, before cbind
  if (is.null(ref_genes)) {
    ref_genes <- rownames(m)
  } else {
    if (!setequal(rownames(m), ref_genes)) {
      stop(sprintf("%s: gene set differs from the first sample", id))
    }
    m <- m[ref_genes, , drop = FALSE]
  }

  # barcodes repeat across 10x runs -- prefix with the sample id
  colnames(m) <- paste(id, colnames(m), sep = "_")

  md <- obj@meta.data
  rownames(md) <- colnames(m)
  if (!"sample_id" %in% colnames(md)) md$sample_id <- id

  mats[[i]] <- m
  mds[[i]]  <- md

  message(sprintf("  [%2d/%2d] %-8s %7d nuclei", i, length(files), id, ncol(m)))

  rm(obj, m, md)
  gc(verbose = FALSE)
}

names(mats) <- names(mds) <- sub("\\.rds$", "", basename(files))

# ---- harmonise metadata ------------------------------------------------------

common  <- Reduce(intersect, lapply(mds, colnames))
all_col <- Reduce(union,     lapply(mds, colnames))
dropped <- setdiff(all_col, common)

if (length(dropped)) {
  message("Metadata columns not present in every sample, dropped: ",
          paste(dropped, collapse = ", "))
}
message("Metadata carried through: ", paste(common, collapse = ", "))

md_all <- do.call(rbind, lapply(mds, function(d) d[, common, drop = FALSE]))
rm(mds); gc(verbose = FALSE)

# ---- cbind -------------------------------------------------------------------

message("Binding ", length(mats), " matrices ...")
counts <- do.call(cbind, mats)
rm(mats); gc(verbose = FALSE)

if (!inherits(counts, "dgCMatrix")) {
  counts <- methods::as(counts, "CsparseMatrix")
}

stopifnot(identical(colnames(counts), rownames(md_all)))
message(sprintf("Merged matrix: %d genes x %d nuclei (%.2f GB)",
                nrow(counts), ncol(counts),
                as.numeric(object.size(counts)) / 1024^3))

# ---- gene filter -------------------------------------------------------------
# count non-zeros per row straight from the sparse structure -- avoids
# materialising a 36k x 95k logical matrix

n_cells <- tabulate(counts@i + 1L, nbins = nrow(counts))
keep    <- n_cells >= min_cells_per_gene

message(sprintf("Genes detected in >= %d nuclei: %d of %d (dropping %d)",
                min_cells_per_gene, sum(keep), length(keep), sum(!keep)))

counts <- counts[keep, , drop = FALSE]
gc(verbose = FALSE)

# ---- build object ------------------------------------------------------------

seu <- CreateSeuratObject(
  counts       = counts,
  meta.data    = md_all,
  assay        = assay_out,
  project      = "hippo_snrna",
  min.cells    = 0,
  min.features = 0
)
rm(counts); gc(verbose = FALSE)

DefaultAssay(seu) <- assay_out
Idents(seu)       <- "sample_id"

# make the grouping variables explicit factors so downstream tables behave
for (v in intersect(c("sample_id", "brnum", "round", "sort", "sex",
                      "doublet_class"), colnames(seu@meta.data))) {
  seu@meta.data[[v]] <- factor(seu@meta.data[[v]])
}

message(sprintf("Seurat object: %d genes x %d nuclei (%.2f GB)",
                nrow(seu), ncol(seu),
                as.numeric(object.size(seu)) / 1024^3))

# ---- per-sample summary ------------------------------------------------------

md       <- seu@meta.data
count_col <- paste0("nCount_",   assay_out)
feat_col  <- paste0("nFeature_", assay_out)

summ <- do.call(rbind, lapply(split(md, md$sample_id), function(d) {
  data.frame(
    sample_id      = d$sample_id[1],
    brnum          = if ("brnum" %in% names(d)) as.character(d$brnum[1]) else NA,
    round          = if ("round" %in% names(d)) as.character(d$round[1]) else NA,
    sort           = if ("sort"  %in% names(d)) as.character(d$sort[1])  else NA,
    n_nuclei       = nrow(d),
    median_umi     = median(d[[count_col]]),
    median_genes   = median(d[[feat_col]]),
    doublet_rate   = if ("doublet_class" %in% names(d))
                       mean(d$doublet_class == "doublet") else NA_real_,
    median_contam  = if ("contamination" %in% names(d))
                       median(d$contamination) else NA_real_,
    median_mt      = if ("percent_mt" %in% names(d))
                       median(d$percent_mt) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
summ <- summ[order(summ$sort, summ$sample_id), ]
rownames(summ) <- NULL

print(summ, digits = 3)
write.csv(summ, out_tab, row.names = FALSE)

message(sprintf("\nTotal: %d nuclei, %d flagged as doublets (%.1f%%)",
                ncol(seu),
                sum(md$doublet_class == "doublet"),
                100 * mean(md$doublet_class == "doublet")))

# ---- save --------------------------------------------------------------------

saveRDS(seu, out_obj)
message("Wrote ", out_obj)
message("Wrote ", out_tab)
message(sprintf("Elapsed: %.1f min",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
