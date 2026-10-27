# =============================================================================
# Phase 0 - Environment, sample metadata, and data organisation
# Project : snRNA-seq pipeline on GSE264624 (LIBD human anterior hippocampus)
# Author  : Hendrik
# =============================================================================
#
# Runs top to bottom in a fresh R session. Package installation is guarded, so
# re-running does not reinstall anything already present.
#
# Troubleshooting history is preserved in the APPENDIX at the bottom of this
# file - commented out, with notes on what failed and why. Nothing there is
# meant to be executed.
#
# Recorded environment (see also results/sessionInfo_phase0.txt):
#   R 4.4.1 + Rtools44, Windows 11, 16 GB RAM
#   Bioconductor 3.20
#   Seurat 5.6.0 / SeuratObject 5.4.0
#   humanHippocampus2024 installed from GitHub, not Bioconductor (see Appendix C)
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Project folders
# -----------------------------------------------------------------------------

dirs <- c("data/raw", "data/processed", "data/processed/filtered",
          "R", "analysis", "results", "results/figures", "results/tables")
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))


# -----------------------------------------------------------------------------
# 2. renv
# -----------------------------------------------------------------------------
# Initialised AFTER packages were already installed in the user library, which
# caused the hydrate problems in Appendix B. On a clean machine, run this first.

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (!file.exists("renv.lock")) renv::init()


# -----------------------------------------------------------------------------
# 3. Packages
# -----------------------------------------------------------------------------
# One-time setup. The `need()` guard means re-running this block is a no-op.

need <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()                  # recorded: 3.20

cran_pkgs <- c("Seurat", "Matrix", "harmony", "remotes",
               "ggplot2", "patchwork", "dplyr", "tidyr", "readr")

bioc_pkgs <- c("DropletUtils", "scDblFinder", "celda", "DESeq2",
               "scuttle", "SingleCellExperiment", "fgsea",
               "ExperimentHub", "SummarizedExperiment")

if (length(need(cran_pkgs))) install.packages(need(cran_pkgs))
if (length(need(bioc_pkgs))) BiocManager::install(need(bioc_pkgs))

# GitHub-only packages.
#   presto               - fast Wilcoxon marker tests, used by Seurat in Phase 2
#   humanHippocampus2024 - not on Bioconductor 3.20 (see Appendix C)
if (length(need("presto")))
  remotes::install_github("immunogenomics/presto")
if (length(need("humanHippocampus2024")))
  remotes::install_github("christinehou11/humanHippocampus2024")

# Verification - all must be TRUE before continuing.
pkgs <- c(cran_pkgs, bioc_pkgs, "presto", "humanHippocampus2024")
stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))

renv::snapshot()


# -----------------------------------------------------------------------------
# 4. Sample metadata
# -----------------------------------------------------------------------------
# GEO itself carries no demographics - every GSM record lists only
# "tissue: brain (hippocampus)". The sample annotation lives in the authors'
# processed object instead, distributed via ExperimentHub.
#   EH9605 = spatial (SpatialExperiment)
#   EH9606 = snRNA-seq (SingleCellExperiment)   <- the one we want

library(ExperimentHub)
library(SummarizedExperiment)

ehub <- ExperimentHub()
query(ehub, "humanHippocampus2024")     # snapshotDate recorded: 2024-10-24

sce_ref <- ehub[["EH9606"]]             # ~2 GB on first call, cached after
names(colData(sce_ref))

meta <- unique(as.data.frame(colData(sce_ref)[
  , c("Sample", "brnum", "round", "sort", "age", "sex", "pmi")
]))

# colData returns Sample as a factor; convert before any joining or comparison.
meta$Sample    <- as.character(meta$Sample)

# Files use "1c", the authors' object uses "1c-scp". Anchored with $ so only a
# trailing suffix is stripped.
meta$sample_id <- sub("-scp$", "", meta$Sample)

stopifnot(nrow(meta) == 26)


# -----------------------------------------------------------------------------
# 5. Design checks
# -----------------------------------------------------------------------------
# These determine what can be modelled in Phase 3. Findings are noted inline.

table(meta$brnum)                             # 10 donors; 3 gave 4 libraries, 7 gave 2
table(meta$brnum, meta$sort)                  # every donor has both sorts, balanced
table(meta$sex,   meta$round)                 # sex spread across rounds, not confounded
table(meta$brnum, meta$round)                 # round nearly NESTED within donor
table(unique(meta[, c("brnum", "sex")])$sex)  # donor-level n: 4 F, 6 M
summary(meta$age)                             # 33.4 - 61.5, too narrow for a contrast
summary(meta$pmi)                             # 13.5 - 33.5 h, QC covariate only

# Consequences carried into later phases:
#   - Unit of replication is the DONOR (n = 10), not the library (n = 26).
#   - DE design is ~ sex alone. A round term cannot be fitted on 10 observations,
#     and round is near-nested in donor anyway. Costs power, not validity.
#   - `sort` (PI+ vs PI+NeuN+) is a DESIGN variable, not a batch effect.
#     Never pass it to Harmony. Balanced, so no covariate needed at pseudobulk.
#   - Harmony integrates over brnum and round only.


# -----------------------------------------------------------------------------
# 6. Extract and reorganise the GEO archive
# -----------------------------------------------------------------------------
# GEO ships 78 files flat with GSM prefixes. Read10X() needs one directory per
# sample containing exactly barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz.

src       <- "data/raw/extracted"
dest_root <- "data/raw/samples"

if (!dir.exists(src)) {
  untar("data/raw/GSE264624_RAW.tar", exdir = src)
}

f <- list.files(src)
stopifnot(length(f) == 78)

# Suffix stripped from the end rather than splitting on "-", so a sample ID
# containing its own hyphen would not silently break.
info <- data.frame(
  file = f,
  gsm  = sub("^(GSM[0-9]+)_.*$", "\\1", f),
  id   = sub("^GSM[0-9]+_(.*)-(barcodes\\.tsv|features\\.tsv|matrix\\.mtx)\\.gz$", "\\1", f),
  type = sub("^.*-(barcodes\\.tsv|features\\.tsv|matrix\\.mtx)\\.gz$", "\\1", f),
  stringsAsFactors = FALSE
)

stopifnot(length(unique(info$id)) == 26)
stopifnot(all(table(info$id) == 3))

new_name <- c("barcodes.tsv" = "barcodes.tsv.gz",
              "features.tsv" = "features.tsv.gz",
              "matrix.mtx"   = "matrix.mtx.gz")

# Copy rather than move, so data/raw/extracted survives as a fallback.
for (i in seq_len(nrow(info))) {
  d <- file.path(dest_root, info$id[i])
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.copy(from      = file.path(src, info$file[i]),
            to        = file.path(d, new_name[info$type[i]]),
            overwrite = TRUE)
}

# Every sample folder must contain exactly 3 files.
stopifnot(all(sapply(list.dirs(dest_root, recursive = FALSE),
                     function(d) length(list.files(d))) == 3))


# -----------------------------------------------------------------------------
# 7. Reconcile IDs, then write the metadata table
# -----------------------------------------------------------------------------
# Both must be empty: every folder has metadata, every metadata row has files.

stopifnot(length(setdiff(info$id, meta$sample_id)) == 0)
stopifnot(length(setdiff(meta$sample_id, info$id)) == 0)

write.csv(meta, "data/processed/sample_metadata.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 8. Helper functions
# -----------------------------------------------------------------------------
# Move to R/functions.R and source() once a third script needs them.

library(Seurat)
library(Matrix)

read_sample <- function(sid, root = "data/raw/samples") {
  m <- Read10X(file.path(root, sid))
  # 10x barcodes come from a fixed whitelist, so the same sequence appears in
  # every sample. Prefixing at read time prevents collisions on merge.
  colnames(m) <- paste0(sid, "_", colnames(m))
  m
}

summarise_sample <- function(m, sid) {
  tot <- Matrix::colSums(m)
  data.frame(
    sample_id   = sid,
    n_genes     = nrow(m),
    n_barcodes  = ncol(m),
    total_umi   = sum(tot),
    bc_over_500 = sum(tot > 500),
    max_umi     = max(tot),
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# 9. First look at one sample
# -----------------------------------------------------------------------------
# GEO deposited RAW, unfiltered matrices: 36,601 genes x ~1.6M barcodes.
# Almost all barcodes are empty droplets - empty droplet removal is ours to do.

m   <- read_sample("1c")
dim(m)                                  # 36601 x 1616924
m[1:5, 1:5]
sum(m) / 1e6                            # 95.6M total UMIs

tot <- Matrix::colSums(m)
summary(tot)                            # median 1 UMI - the matrix is mostly empty
sum(tot == 0)                           # 382,913 barcodes with nothing
sum(tot > 500)                          # 6,712 plausible nuclei, rough cut

# Barcode rank plot: plateau of real nuclei, sharp cliff, then the ambient
# plateau at 200-300 UMIs. That ambient level sits close to a weak nucleus at
# ~1,000 UMIs, which is the quantitative case for decontX in Phase 1.
r <- sort(tot, decreasing = TRUE)

png("results/figures/barcode_rank_1c.png", width = 900, height = 700, res = 120)
plot(seq_along(r), r, log = "xy", type = "l",
     xlab = "barcode rank", ylab = "total UMI",
     main = "Sample 1c - barcode ranks")
dev.off()

rm(m, tot, r); gc()


# -----------------------------------------------------------------------------
# 10. Multi-sample read test and memory budget
# -----------------------------------------------------------------------------
# Four samples only. Confirms consistent gene counts and measures the memory
# footprint before committing to all 26.

test_ids <- meta$sample_id[1:4]
mats  <- list()
stats <- list()

for (sid in test_ids) {
  message("reading ", sid, " ...")
  mats[[sid]]  <- read_sample(sid)
  stats[[sid]] <- summarise_sample(mats[[sid]], sid)
}

stats <- do.call(rbind, stats)
stats
write.csv(stats, "results/tables/raw_sample_stats_test4.csv", row.names = FALSE)

format(object.size(mats), units = "GB")  # 2.3 GB for 4 -> ~15 GB for 26

# CONCLUSION: all 26 raw matrices will not fit in 16 GB. Phase 1 must process
# one sample at a time - read, filter, write .rds, then rm() and gc(). The loop
# should skip samples whose output already exists so it is restartable.

rm(mats); gc()


# -----------------------------------------------------------------------------
# 11. Session record
# -----------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "results/sessionInfo_phase0.txt")
renv::snapshot()


# =============================================================================
# APPENDIX - troubleshooting log
# =============================================================================
# Preserved for traceability. None of this should be run. Each block records a
# dead end hit during setup, in the order encountered.
#
#
# A. renv was never installed
# ---------------------------
# The "Use renv with this project" box was missed when creating the project, so
# renv::snapshot() failed with "there is no package called 'renv'". Packages had
# already gone into the user library rather than a project-local one.
#
#
# B. Stale 00LOCK directory blocked the renv install
# --------------------------------------------------
# install.packages("renv") failed with:
#   ERROR: failed to lock directory '.../win-library/4.4' for modifying
#   Try removing '.../win-library/4.4/00LOCK-renv'
#
# R creates 00LOCK-<pkg> during installation and removes it on success. An
# earlier interrupted install had left it behind. First attempt did not clear it:
#
#   lib <- "C:/Users/Hendrik/AppData/Local/R/win-library/4.4"
#   unlink(file.path(lib, "00LOCK-renv"), recursive = TRUE)        # no effect
#
# Succeeded with force = TRUE:
#
#   unlink(file.path(lib, "00LOCK-renv"), recursive = TRUE, force = TRUE)
#
# Alternative that sidesteps the mechanism entirely:
#   install.packages("renv", INSTALL_opts = "--no-lock")
#
#
# C. renv::init() left the project library incomplete
# ---------------------------------------------------
# After init() restarted the session, .libPaths() pointed at the project library
# and Seurat was gone. renv::hydrate() pulled the CRAN packages - but downloaded
# current versions from CRAN rather than linking the cache, so it compiled for
# ~30 minutes. Bioconductor packages were NOT covered and had to be reinstalled
# by hand:
#
#   renv::hydrate()
#   BiocManager::install(c("scDblFinder", "celda", "DESeq2"))
#   install.packages("harmony")
#
# During hydrate, R offered to update base/recommended packages
# (Matrix, MASS, mgcv, survival, ...). DECLINED - answer "n". Updating Matrix
# breaks packages compiled against its C API (irlba, used by Seurat's RunPCA)
# with errors like "function 'as_cholmod_sparse' not provided by package 'Matrix'",
# which surface much later and look nothing like a packaging problem.
#
# Lesson: initialise renv BEFORE installing anything. Section 2 above reflects
# that ordering even though it is not what happened here.
#
#
# D. humanHippocampus2024 is not available for Bioconductor 3.20
# --------------------------------------------------------------
#   BiocManager::install("humanHippocampus2024")
#   Warning: package 'humanHippocampus2024' is not available for
#            Bioconductor version '3.20'
#
# The package entered Bioconductor at release 3.21, which requires R 4.5. This
# machine runs R 4.4.1 and is therefore pinned to Bioc 3.20. Reinstalling cannot
# fix it.
#
# Resolved by installing the package from GitHub instead. The wrapper is pure R
# with nothing compiled, and the ExperimentHub records happen to predate the
# 3.20 snapshot (snapshotDate 2024-10-24), so EH9605/EH9606 resolve normally:
#
#   remotes::install_github("christinehou11/humanHippocampus2024")
#
# This is luck, not design. An R upgrade to 4.5+ will become necessary for other
# packages eventually - deferred because the project was mid-sprint.
#
#
# E. First renv::snapshot() captured almost nothing
# -------------------------------------------------
# snapshot() scans R/ for library() calls and records only what the code uses.
# Running it before any scripts existed produced a near-empty lockfile. Not a
# bug. The lockfile grows as scripts are written; re-snapshot at the end of each
# phase. renv::snapshot(type = "all") captures the whole library regardless.
# =============================================================================
