## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo      = TRUE,
  message   = FALSE,
  warning   = FALSE,
  fig.align = "center",
  # html_vignette inherits fig.retina = 2 from rmarkdown, which renders every
  # figure at twice the nominal resolution and then scales it down in the
  # browser. The pixel count, and therefore the size of the base64 blob
  # embedded in the self-contained HTML, is four times larger for no visible
  # gain in a vignette. Setting it to 1 is what keeps the installed size of
  # doc/ within the limit R CMD check reports on.
  fig.retina = 1,
  dpi        = 72,
  fig.width = 7,
  fig.height = 6
)

## ----packages-----------------------------------------------------------------
library(TmCalculator)
library(BSgenome)
library(GenomicRanges)

## ----available-genomes, eval=FALSE--------------------------------------------
# BSgenome::available.genomes()          # is your assembly already packaged?
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

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

## ----setup-ecoli-bsgenome, message=FALSE, warning=FALSE-----------------------
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

## ----tm-profile---------------------------------------------------------------
runtime <- system.time({
  tm_ecoli <- tm_calculate(
    input_seq = genome_name,   # a BSgenome package name, a FASTA path, or sequences
    regions  = chr_name,     # the whole chromosome; "U00096.3:1-100000" would
    unit     = "region",     #   tile a sub-region instead
    window   = 200L,
    slide    = 200L,         # slide == window gives a non-overlapping tiling
    method   = "tm_nn",
    nn_table = "DNA_NN_Breslauer_1986",
    Na       = 50,           # mM; standard PCR-like conditions
    verbose  = FALSE
  )$gr                       # $gr is the profile; the object also carries $options
})

cat(sprintf(
  "Tm profile: %.2f s (elapsed) for %s windows\n",
  runtime[["elapsed"]], format(length(tm_ecoli), big.mark = ",")
))

Tm <- as.data.frame(tm_ecoli[, c("Tm", "GC")])
summary(Tm[, c("Tm", "GC")])

## ----profile-head-------------------------------------------------------------
head(tm_ecoli, 3)

## ----tm-calculate-seqs--------------------------------------------------------
tm_calculate(c("ACGTGCTAGCTAGCTAGC", "GGCCATATATGCGC"), method = "tm_nn", Na = 50)

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

## ----zoom-single, fig.width=10, fig.height=8, fig.cap="Zoomed linear view of the 0.1-0.5 Mb region. Seven tracks share one panel, so `track.gap` separates them and `axis.cex` shrinks the tick labels; without both the y-axis labels of adjacent tracks overlap."----
plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks,
  zoom        = "U00096.3:100000-500000",
  ## Seven tracks in a linear panel leave each one little vertical room. The
  ## gap is relative to the panel, so it has to grow with the track count.
  track.gap   = 0.03,
  axis.cex    = 0.55
)

## ----zoom-multi-circular, fig.width=10, fig.height=9, fig.cap="Two zoomed regions, 0.1-0.5 Mb and 3.6-4.5 Mb, concatenated around the circle with a gap between them."----
plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks,
  circular    = TRUE,
  zoom        = c("U00096.3:100000-500000",
                  "U00096.3:3600000-4500000")
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

## ----gc-vs-tm-data------------------------------------------------------------
## The GC-content prediction for the same windows, for comparison.
tm_gc_ecoli <- tm_calculate(genome_name, regions = chr_name, unit = "region",
                            window = 200L, slide = 200L,
                            method = "tm_gc", variant = "Schildkraut1965",
                            Na = 50, verbose = FALSE)$gr
GCmod <- as.data.frame(tm_gc_ecoli[, c("Tm", "GC")])

## Keep only windows whose GC percentage is an exact integer. A 200 bp
## window can only take GC values in steps of 0.5%, so binning 49.5% with
## 50.0% would put a real GC difference inside a group and some of the
## spread below could be attributed to GC rather than to composition.
## The same test also removes windows that contained an ambiguous base,
## whose GC denominator is smaller than the window.
prep <- function(d) {
  d <- d[!is.na(d$Tm) & !is.na(d$GC) & d$width == 200L, ]
  d <- d[d$GC == round(d$GC), ]
  d$GCi <- as.integer(d$GC)
  d
}
Cnn <- prep(as.data.frame(tm_ecoli[, c("Tm", "GC")]))
Cgc <- prep(GCmod)

cnt <- table(Cnn$GCi)
lv  <- sort(as.integer(names(cnt[cnt >= 30])))   # enough windows to be stable
Cnn <- Cnn[Cnn$GCi %in% lv, ]
Cgc <- Cgc[Cgc$GCi %in% lv, ]

modal   <- lv[which.max(cnt[as.character(lv)])]
at_mode <- Cnn$Tm[Cnn$GCi == modal]

cat(sprintf("At %d%% GC: %s windows, Tm from %.2f to %.2f (spread %.2f C)\n",
            modal, format(length(at_mode), big.mark = ","),
            min(at_mode), max(at_mode), diff(range(at_mode))))
cat(sprintf("GC-content formula predicts a single value: %.2f C\n",
            median(Cgc$Tm[Cgc$GCi == modal])))
cat(sprintf("Largest spread at any GC value: %.2f C\n",
            max(tapply(Cnn$Tm, Cnn$GCi, function(v) diff(range(v))))))

## ----gc-vs-tm-plot, fig.width=7, fig.height=5, fig.cap="Melting temperature against GC content for 200 bp windows of the E. coli chromosome. Each box holds windows of identical length, salt and GC content. The red line is the GC-content formula, which by construction has no width within a group; the boxes do, and that width is the composition effect the formula cannot represent."----
gc_line <- tapply(Cgc$Tm, Cgc$GCi, median)[as.character(lv)]

op <- par(mar = c(4.6, 4.6, 1.2, 1.2), las = 1, mgp = c(2.9, 0.7, 0))
boxplot(Tm ~ factor(GCi, levels = lv), data = Cnn, at = seq_along(lv),
        outline = FALSE, border = "#34495E", col = "#D6E4F0", lwd = 0.7,
        xaxt = "n", bty = "n", xlab = "GC content (%)",
        ylab = expression(paste("Melting temperature (", degree, "C)")))
sel <- seq(1, length(lv), by = 5)          # one label per box is unreadable
axis(1, at = seq_along(lv)[sel], labels = lv[sel])
lines(seq_along(lv), gc_line, col = "#C0392B", lwd = 2.2)
legend("topleft", bty = "n", cex = 0.85,
       legend = c("Nearest-neighbour (Breslauer 1986)",
                  "GC-content formula (Schildkraut 1965)"),
       col = c("#34495E", "#C0392B"), lwd = c(0.7, 2.2), seg.len = 1.4)
par(op)

## ----make-figure3, eval=FALSE-------------------------------------------------
# script <- system.file("scripts", "make_figure3.R", package = "TmCalculator")
# file.edit(script)      # --locus-kb, --n-groups, --pairs-in-peaks, --dpi

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
  gr_tm          = tm_ecoli,
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

