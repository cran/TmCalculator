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
library(TmCalculator)
library(BSgenome)
library(GenomicRanges)

## ----forge, eval=FALSE--------------------------------------------------------
# library(BSgenomeForge)
# 
# forgeBSgenomeDataPkgFromNCBI(
#   assembly_accession = "GCF_000005845.2",
#   pkg_maintainer     = "Junhui Li <ljh.biostat@gmail.com>",
#   destdir            = "."
# )
# 
# install.packages(
#   "./BSgenome.Ecoli.NCBI.ASM584v2",
#   repos = NULL,
#   type  = "source"
# )

## ----setup-ecoli-bsgenome, echo=FALSE, message=FALSE, warning=FALSE-----------
ecoli_pkg   <- "BSgenome.Ecoli.NCBI.ASM584v2"
genome_obj  <- "Ecoli"   # BSgenomeObjname in DESCRIPTION; not the package name

.ecoli_genome_ready <- function() {
  if (!requireNamespace(ecoli_pkg, quietly = TRUE)) return(FALSE)
  exists(genome_obj, envir = asNamespace(ecoli_pkg), inherits = FALSE)
}

if (!requireNamespace(ecoli_pkg, quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    utils::install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github(
    "JunhuiLi1017/BSgenome.Ecoli.NCBI.ASM584v2",
    upgrade = "never",
    quiet = TRUE
  )
}

if (!.ecoli_genome_ready()) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("BSgenomeForge", quietly = TRUE)) {
    BiocManager::install("BSgenomeForge", ask = FALSE, update = FALSE)
  }
  pkgdir <- BSgenomeForge::forgeBSgenomeDataPkgFromNCBI(
    assembly_accession = "GCF_000005845.2",
    pkg_maintainer     = "Junhui Li <ljh.biostat@gmail.com>",
    destdir            = tempdir()
  )
  utils::install.packages(pkgdir, repos = NULL, type = "source", quiet = TRUE)
  if (ecoli_pkg %in% loadedNamespaces()) {
    unloadNamespace(ecoli_pkg)
  }
}

if (!.ecoli_genome_ready()) {
  stop(
    "Could not load genome object '", genome_obj, "' from package '", ecoli_pkg, "'.\n",
    "Run the manual 'forge' chunk above, or forgeBSgenomeDataPkgFromNCBI() locally.",
    call. = FALSE
  )
}

## ----load-genome--------------------------------------------------------------
ecoli_pkg  <- "BSgenome.Ecoli.NCBI.ASM584v2"
genome_obj <- "Ecoli"

suppressPackageStartupMessages(library(ecoli_pkg, character.only = TRUE))
genome      <- base::get(genome_obj, envir = asNamespace(ecoli_pkg))
genome_name <- ecoli_pkg
chr_name    <- "U00096.3"
chr_length  <- length(genome[[chr_name]])

cat("Chromosome:", chr_name, "\n")
cat("Length:    ", format(chr_length, big.mark = ","), "bp\n")

## ----make-coord---------------------------------------------------------------
bins_gc <- make_genomiccoord(
  bsgenome    = genome_name,
  chromosomes = chr_name,
  window      = 200L,
  slide       = 200L,
  start       = 1,
  end         = chr_length,
  strand      = "+"
)

cat("Total windows:", length(bins_gc), "\n")

## ----coor-to-gr---------------------------------------------------------------
input_new <- list(pkg_name = genome_name, seq = bins_gc)
runtime1 <- system.time({
  gr_batch <- to_genomic_ranges_fast(input_new)
})

cat(sprintf(
  "Coordinate resolution: %.2f s (elapsed)\n",
  runtime1["elapsed"]
))

## ----tm-calculate-------------------------------------------------------------
runtime2 <- system.time({
  tm_ASM584v2 <- tm_calculate(
    gr_batch,
    method   = "tm_nn",
    nn_table = "DNA_NN_SantaLucia_2004",
    Na       = 50            # mM; standard PCR-like conditions
  )
})

cat(sprintf(
  "Tm calculation: %.2f s (elapsed) for %s windows\n",
  runtime2["elapsed"],
  format(length(bins_gc), big.mark = ",")
))

Tm <- as.data.frame(tm_ASM584v2$gr[, c("Tm", "GC")])
summary(Tm[, c("Tm", "GC")])

## ----track-list---------------------------------------------------------------
# Reference labels: replication origin (ori) and terminus (dif)
label <- data.frame(
  seqnames = genome_name,
  start    = c(3925804, 1590777),
  end      = c(3925804, 1590777),
  label    = c("ori", "dif")
)
data(ecoli_rep_hotspots)
tracks <- list(
  # Ideogram: MutL-AR peaks shown inside the chromosome bar
  list(type = "rect", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
       col = "#2C3E50", bg.col = "grey", name = "MutL-AR",
       legend_font_col = "#2C3E50", ideogram = TRUE, height = 0.5),

  # Sequence thermodynamics
  list(type = "line", data = Tm, value_col = "GC",
       name = "GC content", col = "#4A90E2",
       legend_font_col = "#4A90E2"),
  list(type = "line", data = Tm, value_col = "Tm",
       name = "Melting temp", col = "#E06666",
       legend_font_col = "#E06666", height = 2),

  # Repeat / structural features
  list(type = "line", data = ecoli_rep_hotspots$bins_rep,
       value_col = "count", name = "Microsatellites", col = "#2ECC71",
       legend_font_col = "#2ECC71"),
  list(type = "line", data = ecoli_rep_hotspots$bins_cru,
       value_col = "count", name = "Cruciform", col = "#3B3E6B",
       legend_font_col = "#3B3E6B"),

  # ssDNA regions
  list(data = ecoli_rep_hotspots$ssdna, name = "ssDNA",
       col = "#8E44AD", legend_font_col = "#8E44AD"),

  # GATC methylation sites
  list(type = "line", data = ecoli_rep_hotspots$bins_gatc,
       value_col = "count", name = "GATC sites", col = "#D35400",
       legend_font_col = "#D35400"),

  # Global highlight: translucent bands at MutL-AR peaks across all tracks
  list(type = "highlight", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
       col = "#F1C40F", alpha = 0.18)
)

## ----circular-full, fig.width=8, fig.height=8, fig.cap="Circular genome map of E. coli K-12 MG1655. Concentric rings from outside in: MutL-AR peaks (grey ideogram), GC content, melting temperature, microsatellite density, cruciform sequences, ssDNA regions, and GATC site density. Yellow highlight bands mark MutL-AR peak regions."----
plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks,
  circular    = TRUE,
  label       = label
)

## ----linear-full, fig.width=10, fig.height=5, fig.cap="Linear genome view of E. coli K-12 MG1655. MutL-AR peaks are drawn inside the chromosome ideogram bar."----
plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks
)

## ----zoom-single, fig.width=10, fig.height=5, fig.cap="Zoomed linear view of the 1.0-2.0 Mb region."----
plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks,
  zoom        = "U00096.3:1000000-2000000"
)

## ----zoom-multi-linear, fig.width=10, fig.height=8, fig.cap="Two zoomed regions displayed as stacked linear panels."----
plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks,
  zoom        = c("U00096.3:1000000-1200000",
                   "U00096.3:3000000-3200000")
)

## ----zoom-multi-circular, fig.width=8, fig.height=8, fig.cap="Two zoomed regions concatenated on a circular plot. Dashed lines separate the regions."----
plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks,
  circular    = TRUE,
  zoom        = c("U00096.3:1000000-1200000",
                   "U00096.3:3000000-3200000")
)

## ----circular-pan, fig.width=7, fig.height=7, fig.cap="Panned circular view showing the upper-right quadrant of the E. coli chromosome."----
plot_genome_track(
  genome_name   = genome_name,
  genome_size   = chr_length,
  track_list    = tracks,
  circular      = TRUE,
  canvas.xlim   = c(0.5, 1),
  canvas.ylim   = c(0,   1),
  circle.margin = c(0.05, 0.05)
)

## ----per-track-highlight, fig.width=8, fig.height=8, fig.cap="Circular plot with per-track highlight bands on the Microsatellites track."----
tracks_hl <- tracks
tracks_hl[[4]]$highlight <- list(
  data  = ecoli_rep_hotspots$bins_rep[1100:1200, ],
  col   = "black",
  alpha = 0.12
)

plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks_hl,
  circular    = TRUE
)

## ----build-mutH-peaks---------------------------------------------------------
mutH_peaks <- GRanges(
  seqnames = ecoli_rep_hotspots$all_peaks_IP_mutH$chr,
  ranges   = IRanges(start = ecoli_rep_hotspots$all_peaks_IP_mutH$start,
                     end   = ecoli_rep_hotspots$all_peaks_IP_mutH$end)
)
seqlevels(mutH_peaks) <- "U00096.3"

## ----annotate-mutH-peaks------------------------------------------------------
mutH_peaks$peak_id <- paste0("mutH_", seq_along(mutH_peaks))

tm_annot <- integrate_granges(
  gr_tm          = tm_ASM584v2$gr,
  gr_features    = mutH_peaks,
  strategy       = "overlap",
  feature_cols   = "peak_id",
  keep_unmatched = TRUE
)

tm_annot$in_mutH <- ifelse(is.na(tm_annot$peak_id), "non_peak", "peak")
table(tm_annot$in_mutH)

## ----wilcoxon-----------------------------------------------------------------
res <- compare_groups(
  gr          = tm_annot,
  target      = c("Tm", "GC"),
  method      = "wilcoxon",
  group       = "in_mutH",
  alternative = "greater",
  posthoc     = FALSE
)
res$results
res$summary

## ----session-info-------------------------------------------------------------
sessionInfo()

