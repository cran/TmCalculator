## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo      = TRUE,
  message   = FALSE,
  warning   = FALSE,
  fig.align = "center",
  fig.width = 7,
  fig.height = 6
)

## ----packages-----------------------------------------------------------------
# Bioconductor core
#BiocManager::install(c(
#  "BSgenomeForge",   # forge custom BSgenome packages from NCBI
#  "BSgenome",
#  "GenomicRanges"
#))

library(TmCalculator)
library(BSgenome)
library(GenomicRanges)

## ----forge, eval=FALSE--------------------------------------------------------
# library(BSgenomeForge)
# 
# # Download the FASTA and build the data package (takes ~1–2 min)
# forgeBSgenomeDataPkgFromNCBI(
#   assembly_accession = "GCF_000005845.2",
#   pkg_maintainer     = "Junhui Li <ljh.biostat@gmail.com>",
#   destdir            = "."        # writes BSgenome.Ecoli.NCBI.ASM584v2/ here
# )
# 
# # Install the freshly forged package
# install.packages(
#   "./BSgenome.Ecoli.NCBI.ASM584v2",
#   repos = NULL,
#   type  = "source"
# )

## ----setup-ecoli-bsgenome, eval=FALSE, echo=FALSE, message=FALSE, warning=FALSE----
# ecoli_pkg   <- "BSgenome.Ecoli.NCBI.ASM584v2"
# genome_obj  <- "Ecoli"   # BSgenomeObjname in DESCRIPTION; not the package name
# 
# .ecoli_genome_ready <- function() {
#   if (!requireNamespace(ecoli_pkg, quietly = TRUE)) return(FALSE)
#   exists(genome_obj, envir = asNamespace(ecoli_pkg), inherits = FALSE)
# }
# 
# if (!requireNamespace(ecoli_pkg, quietly = TRUE)) {
#   if (!requireNamespace("remotes", quietly = TRUE)) {
#     utils::install.packages("remotes", repos = "https://cloud.r-project.org")
#   }
#   remotes::install_github(
#     "JunhuiLi1017/BSgenome.Ecoli.NCBI.ASM584v2",
#     upgrade = "never",
#     quiet = TRUE
#   )
# }
# 
# if (!.ecoli_genome_ready()) {
#   if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
#   }
#   if (!requireNamespace("BSgenomeForge", quietly = TRUE)) {
#     BiocManager::install("BSgenomeForge", ask = FALSE, update = FALSE)
#   }
#   pkgdir <- BSgenomeForge::forgeBSgenomeDataPkgFromNCBI(
#     assembly_accession = "GCF_000005845.2",
#     pkg_maintainer     = "Junhui Li <ljh.biostat@gmail.com>",
#     destdir            = tempdir()
#   )
#   utils::install.packages(pkgdir, repos = NULL, type = "source", quiet = TRUE)
#   if (ecoli_pkg %in% loadedNamespaces()) {
#     unloadNamespace(ecoli_pkg)
#   }
# }
# 
# if (!.ecoli_genome_ready()) {
#   stop(
#     "Could not load genome object '", genome_obj, "' from package '", ecoli_pkg, "'.\n",
#     "Run the manual 'forge' chunk above, or forgeBSgenomeDataPkgFromNCBI() locally.",
#     call. = FALSE
#   )
# }

## ----load-genome, eval=FALSE--------------------------------------------------
# ecoli_pkg  <- "BSgenome.Ecoli.NCBI.ASM584v2"
# genome_obj <- "Ecoli"
# 
# suppressPackageStartupMessages(library(ecoli_pkg, character.only = TRUE))
# genome      <- base::get(genome_obj, envir = asNamespace(ecoli_pkg))
# genome_name <- ecoli_pkg
# chr_name    <- "U00096.3"
# chr_length  <- length(genome[[chr_name]])
# 
# cat("Chromosome:", chr_name, "\n")
# cat("Length:    ", format(chr_length, big.mark = ","), "bp\n")

## ----make-coord, eval=FALSE---------------------------------------------------
# bins_gc <- make_genomiccoord(
#   bsgenome     = genome,
#   chromosomes  = chr_name,
#   window       = 200L,
#   slide        = 200L,          # non-overlapping: slide == window
#   start        = 1,
#   end          = chr_length,
#   strand       = "+"
# )
# 
# cat("Total windows:", length(bins_gc), "\n")
# cat("Last window:  ", bins_gc[length(bins_gc)], "\n")

## ----coor-to-gr, eval=FALSE---------------------------------------------------
# input_new <- list(
#   pkg_name = genome_name,
#   seq      = bins_gc
# )
# 
# runtime1 <- system.time({
#   gr_batch <- to_genomic_ranges_fast(input_new)
# })
# 
# cat(sprintf(
#   "Coordinate resolution: %.2f s (elapsed)\n",
#   runtime1["elapsed"]
# ))

## ----tm-calculate, eval=FALSE-------------------------------------------------
# runtime2 <- system.time({
#   tm_ASM584v2 <- tm_calculate(
#     gr_batch,
#     method   = "tm_nn",
#     nn_table = "DNA_NN_SantaLucia_2004",
#     Na       = 50            # mM; standard PCR-like conditions
#   )
# })
# 
# cat(sprintf(
#   "Tm calculation: %.2f s (elapsed) for %s windows\n",
#   runtime2["elapsed"],
#   format(length(bins_gc), big.mark = ",")
# ))
# 
# 

## ----extract-tm, eval=FALSE---------------------------------------------------
# Tm <- as.data.frame(tm_ASM584v2$gr[, c(4, 6)])
# 
# # Quick summary
# summary(Tm[, c("Tm", "GC")])

## ----load-hotspots, eval=FALSE------------------------------------------------
# data(ecoli_rep_hotspots)
# 
# # Reference labels: replication origin (ori) and terminus (dif)
# label <- data.frame(
#   seqnames = genome_name,
#   start    = c(3925804, 1590777),
#   end      = c(3925804, 1590777),
#   label    = c("ori", "dif")
# )

## ----track-list, eval=FALSE---------------------------------------------------
# tracks <- list(
#   # --- Protein interaction ---
#   list(
#     type   = "rect",
#     data   = ecoli_rep_hotspots$all_peaks_IP_mutH,
#     col    = "black",
#     bg.col = "grey",
#     name   = "MutL-AR"
#   ),
# 
#   # --- Sequence thermodynamics ---
#   list(
#     type      = "line",
#     data      = Tm,
#     value_col = "GC",
#     name      = "GC content"
#   ),
#   list(
#     type      = "line",
#     data      = Tm,
#     value_col = "Tm",
#     name      = "Melting temp"
#   ),
# 
#   # --- Repeat / structural features ---
#   list(
#     type      = "line",
#     data      = ecoli_rep_hotspots$bins_rep,
#     value_col = "count",
#     name      = "Microsatellites"
#   ),
#   list(
#     type      = "line",
#     data      = ecoli_rep_hotspots$bins_cru,
#     value_col = "count",
#     name      = "Cruciform"
#   ),
# 
#   # --- ssDNA regions ---
#   list(
#     data = ecoli_rep_hotspots$ssdna,
#     name = "ssDNA"
#   ),
# 
#   # --- GATC methylation sites ---
#   list(
#     type      = "line",
#     data      = ecoli_rep_hotspots$bins_gatc,
#     value_col = "count",
#     name      = "GATC sites"
#   )
# )

## ----circos-full, eval=FALSE, fig.width=8, fig.height=8, fig.cap="Genome-wide Tm profile of *E. coli* K-12 MG1655 (ASM584v2). Concentric rings from outside in: MutL-AR ChIP-seq peaks (grey), GC content (%), melting temperature (°C), microsatellite density, cruciform-forming sequences, ssDNA regions, and GATC site density. Replication origin (ori, position 3,925,804) and terminus (dif, position 1,590,777) are labelled."----
# plot_circos_genome(
#   genome_name = genome_name,
#   genome_size = chr_length,
#   track_list  = tracks,
#   canvas.xlim = c(-1, 1),
#   canvas.ylim = c(-1, 1),
#   label       = label
# )

## ----circos-zoom, eval=FALSE, fig.width=7, fig.height=7, fig.cap="Zoomed circular view of the upper-right quadrant of the *E. coli* chromosome. Track order and colour scheme are identical to the full-genome plot above."----
# plot_circos_genome(
#   genome_name = genome_name,
#   genome_size = chr_length,
#   track_list  = tracks,
#   canvas.xlim = c(0.5, 1),
#   canvas.ylim = c(0,   1)
# )

## ----build-mutH-peaks, eval=FALSE---------------------------------------------
# mutH_peaks <- GRanges(
#   seqnames = ecoli_rep_hotspots$all_peaks_IP_mutH$chr,         # NC_000913.3
#   ranges   = IRanges(start = ecoli_rep_hotspots$all_peaks_IP_mutH$start,
#                      end   = ecoli_rep_hotspots$all_peaks_IP_mutH$end)
# )

## ----align-seqlevels, eval=FALSE----------------------------------------------
# # tm_ASM584v2 uses the BSgenome name "U00096.3"; remap the peaks to match.
# seqlevels(mutH_peaks) <- "U00096.3"

## ----annotate-mutH-peaks, eval=FALSE------------------------------------------
# mutH_peaks$peak_id <- paste0("mutH_", seq_along(mutH_peaks))
# 
# tm_annot <- integrate_granges(
#   gr_tm        = tm_ASM584v2$gr,
#   gr_features  = mutH_peaks,
#   strategy     = "overlap",
#   feature_cols = "peak_id",
#   keep_unmatched = TRUE   # keep non-overlapping tiles so they form the "non_peak" group
# )

## ----add-two-level-grouping-column, eval=FALSE--------------------------------
# tm_annot$in_mutH <- ifelse(is.na(tm_annot$peak_id), "non_peak", "peak")

## ----sanity-check, eval=FALSE-------------------------------------------------
# table(tm_annot$in_mutH)

## ----wilcoxon-on-tm-and-gc, eval=FALSE----------------------------------------
# res <- compare_groups(
#   gr           = tm_annot,
#   target       = c("Tm", "GC"),
#   method       = "wilcoxon",
#   group        = "in_mutH",
#   alternative  = "greater",   # tests "peak > non_peak"
#   posthoc      = FALSE     # only 2 groups, no post-hoc needed
# )
# res$results   # one row per target: test, statistic, p.value
# res$summary   # per-group n, mean, sd, median

## ----session-info-------------------------------------------------------------
sessionInfo()

