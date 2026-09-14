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

## ----make-coord---------------------------------------------------------------
runtime0 <- system.time({
  bins_gc <- make_genomiccoord(
    bsgenome    = genome_name,
    chromosomes = chr_name,
    window      = 200L,
    slide       = 200L,
    start       = 1,
    end         = chr_length,
    strand      = "+"
  )
})

cat("Total windows:", length(bins_gc), "\n")
cat(sprintf("Window generation: %.2f s (elapsed)\n", runtime0["elapsed"]))

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
    nn_table = "DNA_NN_Breslauer_1986",
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
tm_gc_ecoli <- tm_calculate(gr_batch, method = "tm_gc",
                            variant = "Schildkraut1965", Na = 50)
GCmod <- as.data.frame(tm_gc_ecoli$gr[, c("Tm", "GC")])

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
Cnn <- prep(as.data.frame(tm_ASM584v2$gr[, c("Tm", "GC")]))
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

## ----find-locus---------------------------------------------------------------
W <- as.data.frame(tm_ASM584v2$gr[, c("Tm", "GC")])
W <- W[!is.na(W$Tm) & !is.na(W$GC) & W$width == 200L, ]
W <- W[order(W$start), ]

## Which windows lie inside a MutL-AR peak. Matched on coordinates alone:
## the peak table and the window table both describe one contig, and
## matching on sequence names would fail silently if one calls it by the
## accession and the other by the genome package name.
peaks <- as.data.frame(ecoli_rep_hotspots$all_peaks_IP_mutH)
W$in_peak <- IRanges::overlapsAny(
  IRanges::IRanges(W$start, W$end),
  IRanges::IRanges(as.numeric(peaks$start), as.numeric(peaks$end)))

## The k GC values inside a locus whose windows are furthest apart in Tm.
## Each group contributes its two extreme windows; shading every window at
## that GC would fill the panel with stripes and hide the tracks.
top_groups <- function(d, k = 3L) {
  g  <- split(seq_len(nrow(d)), d$GC)
  g  <- g[lengths(g) >= 2L]
  if (!length(g)) return(list())
  sp <- vapply(g, function(ix) diff(range(d$Tm[ix])), numeric(1))
  o  <- order(sp, decreasing = TRUE)[seq_len(min(k, length(g)))]
  lapply(o, function(j) {
    ix <- g[[j]]
    list(gc = as.numeric(names(g)[j]), spread = unname(sp[j]),
         rows = ix[c(which.min(d$Tm[ix]), which.max(d$Tm[ix]))])
  })
}

n_loc  <- 100L                       # 100 x 200 bp = 20 kb
starts <- seq_len(nrow(W) - n_loc + 1L)

## Two restrictions on the candidate loci.
##
## Contiguity: a locus is n_loc consecutive ROWS of W, and W has had short
## and ambiguous windows dropped, so consecutive rows need not be adjacent
## on the chromosome. A locus spanning such a gap would be drawn with an
## axis covering far more than 20 kb and would quietly stop being a zoom.
##
## Peaks: the panel is meant to say something about the regions this study
## is about, not only about the model, so only loci centred on a MutL-AR
## peak are considered.
contig <- (W$start[seq(n_loc, nrow(W))] - W$start[starts]) == (n_loc - 1L) * 200L
starts <- starts[contig]
starts <- starts[W$in_peak[pmin(starts + n_loc %/% 2L, nrow(W))]]

## Scored on the SUM of the top three spreads rather than the single best
## one: a locus with one spectacular pair and nothing around it reads as a
## peculiar sequence, whereas three groups at three GC values read as a
## property of the model.
sc <- vapply(starts, function(i)
  sum(vapply(top_groups(W[i:(i + n_loc - 1L), ]), function(g) g$spread,
             numeric(1))), numeric(1))

i0 <- starts[which.max(sc)]
d  <- W[i0:(i0 + n_loc - 1L), ]
gs <- top_groups(d)

do.call(rbind, lapply(gs, function(g)
  data.frame(GC = g$gc, Tm_low = min(d$Tm[g$rows]), Tm_high = max(d$Tm[g$rows]),
             difference = g$spread)))

## ----zoom-locus, fig.width=9, fig.height=4.5, fig.cap="High-magnification view of a 20 kb locus centred on a MutL-AR peak. Three pairs of windows are shaded, one pair per colour; the two windows of a pair have identical length and identical GC content and differ only in the arrangement of their bases. Yellow marks the MutL-AR peaks. Only the GC and Tm tracks are drawn: seven tracks in one linear panel leave each of them too little height to read at this scale."----
grp_cols <- c("#B7791F", "#2E8B57", "#7D3C98")

tracks_zoom <- c(
  list(
    ## The MutL-AR peaks are drawn as the ideogram, in the same black used
    ## in the whole-genome map, so the three panels agree on what a peak
    ## looks like. As the ideogram it also names the chromosome bar, which
    ## is why `genome_name` below is the track name rather than the contig.
    list(type = "rect", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
         col = "black", bg.col = "grey", name = "MutL-AR",
         legend_font_col = "black", ideogram = TRUE, height = 0.5),
    list(type = "line", data = Tm, value_col = "GC", name = "GC",
         col = "#4A90E2", legend_font_col = "#4A90E2", height = 1),
    list(type = "line", data = Tm, value_col = "Tm", name = "Tm",
         col = "#E06666", legend_font_col = "#E06666", height = 1.4),
    list(type = "highlight", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
         col = "#F1C40F", alpha = 0.18)),
  lapply(seq_along(gs), function(i)
    list(type = "highlight", data = d[gs[[i]]$rows, c("seqnames", "start", "end")],
         col = grp_cols[i], alpha = 0.40)))

plot_genome_track(
  genome_name = "MutL-AR",
  genome_size = chr_length,
  track_list  = tracks_zoom,
  zoom        = sprintf("%s:%d-%d", chr_name, min(d$start), max(d$end)),
  track.gap   = 0.06,
  axis.cex    = 0.8,
  legend.show = FALSE,
  ## Without this the panel carries no base positions at all: the default
  ## tick spacing is 500 kb for any view under 10 Mb, which places no tick
  ## inside a 20 kb window.
  base.tick.dist  = 5000,
  base.tick.units = TRUE
)

## plot_genome_track() drops highlight entries before building its legend,
## so the shaded groups would otherwise be unlabelled.
legend("topright", bty = "n", cex = 0.8, border = NA,
       legend = c("GC", "Tm", "MutL-AR peak",
                  vapply(gs, function(g)
                    sprintf("GC %.1f%%:  Tm %.1f-%.1f", g$gc,
                            min(d$Tm[g$rows]), max(d$Tm[g$rows])),
                    character(1))),
       fill = c("#4A90E2", "#E06666",
                adjustcolor("#F1C40F", alpha.f = 0.18),
                adjustcolor(grp_cols, alpha.f = 0.40)),
       text.col = c("#4A90E2", "#E06666", "#B7950B", grp_cols))

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

