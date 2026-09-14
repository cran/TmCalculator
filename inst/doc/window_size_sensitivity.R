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
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)

## ----setup-ecoli-bsgenome, message=FALSE, warning=FALSE-----------------------
ecoli_pkg  <- "BSgenome.Ecoli.NCBI.ASM584v2"
genome_obj <- "Ecoli"   # BSgenomeObjname in DESCRIPTION; not the package name

.ecoli_genome_ready <- function() {
  if (!requireNamespace(ecoli_pkg, quietly = TRUE)) return(FALSE)
  exists(genome_obj, envir = asNamespace(ecoli_pkg), inherits = FALSE)
}

if (!.ecoli_genome_ready()) {
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
  stop(
    "Could not load genome object '", genome_obj, "' from package '",
    ecoli_pkg, "'.\n",
    "See vignette(\"genome_wide_tm_ecoli\") for how to forge it locally.",
    call. = FALSE
  )
}

## ----load-genome--------------------------------------------------------------
suppressPackageStartupMessages(library(ecoli_pkg, character.only = TRUE))
genome      <- base::get(genome_obj, envir = asNamespace(ecoli_pkg))
genome_name <- ecoli_pkg
chr_name    <- "U00096.3"
chr_length  <- GenomeInfoDb::seqlengths(genome)[[chr_name]]

data(ecoli_rep_hotspots)

## MutL-AR peaks, as in the main case study
mutH_peaks <- GRanges(
  seqnames = ecoli_rep_hotspots$all_peaks_IP_mutH$chr,
  ranges   = IRanges(start = ecoli_rep_hotspots$all_peaks_IP_mutH$start,
                     end   = ecoli_rep_hotspots$all_peaks_IP_mutH$end)
)
seqlevels(mutH_peaks) <- chr_name
mutH_peaks$peak_id <- paste0("mutH_", seq_along(mutH_peaks))

## Reference labels: replication origin (ori) and terminus (dif)
label <- data.frame(
  seqnames = genome_name,
  start    = c(3925804, 1590777),
  end      = c(3925804, 1590777),
  label    = c("ori", "dif")
)

## ----window-sensitivity-------------------------------------------------------
window_sizes <- c(50L, 100L, 200L, 500L)

sens <- lapply(window_sizes, function(w) {
  bins_w <- make_genomiccoord(
    bsgenome = genome_name, chromosomes = chr_name,
    window = w, slide = w, start = 1, end = chr_length,
    strand = "+", verbose = FALSE
  )
  gr_w <- to_genomic_ranges_fast(list(pkg_name = genome_name, seq = bins_w))
  tm_w <- tm_calculate(gr_w, method = "tm_nn",
                       nn_table = "DNA_NN_Breslauer_1986", Na = 50)$gr

  ann <- integrate_granges(gr_tm = tm_w, gr_features = mutH_peaks,
                           strategy = "overlap", feature_cols = "peak_id",
                           keep_unmatched = TRUE)
  ann$in_mutH <- ifelse(is.na(ann$peak_id), "non_peak", "peak")
  cg <- compare_groups(gr = ann, target = c("Tm", "GC"),
                       method = "wilcoxon", group = "in_mutH",
                       alternative = "greater", posthoc = FALSE)
  list(gr = tm_w, ann = ann, test = cg)
})
names(sens) <- paste0("w", window_sizes)

## ----sensitivity-distribution-------------------------------------------------
dist_tbl <- do.call(rbind, lapply(seq_along(window_sizes), function(i) {
  g <- sens[[i]]$gr
  data.frame(
    window_bp = window_sizes[i],
    n_windows = length(g),
    Tm_mean   = mean(g$Tm, na.rm = TRUE),
    Tm_sd     = stats::sd(g$Tm, na.rm = TRUE),
    Tm_IQR    = stats::IQR(g$Tm, na.rm = TRUE),
    GC_mean   = mean(g$GC, na.rm = TRUE)
  )
}))
knitr::kable(dist_tbl, digits = 3,
             caption = "Tm/GC distribution by window size.")

## ----sensitivity-density, fig.width=7, fig.height=4.5, fig.cap="Tm distribution at each window size. The distribution shifts to higher Tm as the window lengthens and narrows at the same time; dotted vertical lines mark the group means. Colours match the multi-scale tracks below (darkest = 50 bp, lightest = 500 bp)."----
## Palettes are defined here because they are used both by this panel and by
## the multi-scale track figures further down.
## darkest = finest window (50 bp), lightest = coarsest (500 bp)
scale_cols_gc <- c("#1A5276", "#2E86C1", "#5DADE2", "#AED6F1")
scale_cols_tm <- c("#7B241C", "#CB4335", "#EC7063", "#F5B7B1")

## Densities are normalised, so the ~93,000 windows at 50 bp and the ~9,300
## at 500 bp can be compared directly on one pair of axes.
tm_by_size <- lapply(sens, function(x) x$gr$Tm[is.finite(x$gr$Tm)])
dens <- lapply(tm_by_size, stats::density)

xlim <- range(vapply(dens, function(z) range(z$x), numeric(2)))
ylim <- c(0, max(vapply(dens, function(z) max(z$y), numeric(1))))

op <- par(mar = c(4.2, 4.4, 0.8, 0.8), las = 1)
plot(NA, xlim = xlim, ylim = ylim, bty = "n",
     xlab = expression(italic(T)[m] ~ "(" * degree * "C)"),
     ylab = "Density")
for (i in seq_along(dens)) {
  lines(dens[[i]], col = scale_cols_tm[i], lwd = 2)
  abline(v = dist_tbl$Tm_mean[i], col = scale_cols_tm[i], lwd = 1, lty = 3)
}
legend("topleft", bty = "n", lwd = 2, seg.len = 1.6,
       col = scale_cols_tm[seq_along(window_sizes)],
       legend = sprintf("%d bp: mean %.1f, SD %.2f",
                        dist_tbl$window_bp, dist_tbl$Tm_mean, dist_tbl$Tm_sd))
par(op)

## Spread of the finest window relative to the coarsest. This is the ratio
## quoted in the manuscript, so both come from the same computation.
sd_ratio <- dist_tbl$Tm_sd[1] / dist_tbl$Tm_sd[nrow(dist_tbl)]

## ----sensitivity-correlation--------------------------------------------------
grid_bp <- 1000L
stopifnot(all(grid_bp %% window_sizes == 0L))   # exact aggregation for each

gr_ref  <- sens[["w500"]]$gr
ref_agg <- tapply(gr_ref$Tm, (start(gr_ref) - 1L) %/% grid_bp,
                  mean, na.rm = TRUE)

cor_tbl <- vapply(c("w50", "w100", "w200"), function(k) {
  g      <- sens[[k]]$gr
  agg    <- tapply(g$Tm, (start(g) - 1L) %/% grid_bp, mean, na.rm = TRUE)
  common <- intersect(names(agg), names(ref_agg))
  stats::cor(agg[common], ref_agg[common], use = "complete.obs")
}, numeric(1))
round(cor_tbl, 3)

## ----sensitivity-tests--------------------------------------------------------
## Per-window-size test statistics
test_results <- do.call(rbind, lapply(seq_along(window_sizes), function(i) {
  data.frame(window_bp = window_sizes[i], sens[[i]]$test$results)
}))
test_results

## Per-window-size group summaries, plus the median Tm effect size
test_summary <- do.call(rbind, lapply(seq_along(window_sizes), function(i) {
  ann <- sens[[i]]$ann
  med <- tapply(ann$Tm, ann$in_mutH, stats::median, na.rm = TRUE)
  data.frame(window_bp   = window_sizes[i],
             n_peak      = sum(ann$in_mutH == "peak"),
             Tm_med_diff = unname(med["peak"] - med["non_peak"]),
             sens[[i]]$test$summary)
}))
test_summary

## ----sensitivity-table--------------------------------------------------------
sens_tbl <- do.call(rbind, lapply(seq_along(window_sizes), function(i) {
  g   <- sens[[i]]$gr
  ann <- sens[[i]]$ann
  med <- tapply(ann$Tm, ann$in_mutH, stats::median, na.rm = TRUE)
  res <- sens[[i]]$test$results
  w   <- window_sizes[i]
  data.frame(
    `Window (bp)`        = w,
    Windows              = length(g),
    `Tm mean (C)`        = mean(g$Tm, na.rm = TRUE),
    `Tm SD`              = stats::sd(g$Tm, na.rm = TRUE),
    `GC mean (%)`        = mean(g$GC, na.rm = TRUE),
    ## The reference profile correlates with itself by construction.
    `r vs 1 kb grid`     = if (w == max(window_sizes)) NA_real_
                           else unname(cor_tbl[[paste0("w", w)]]),
    `MutL-AR windows`    = sum(ann$in_mutH == "peak"),
    `Tm peak - bg (C)`   = unname(med["peak"] - med["non_peak"]),
    ## Formatted as character: these p values span sixty-five orders of
    ## magnitude, and rounding them to a fixed number of decimals prints
    ## every one of them as zero.
    `p (Tm)`             = format(res$p.value[res$target == "Tm"],
                                  digits = 2, scientific = TRUE),
    check.names = FALSE, stringsAsFactors = FALSE)
}))
knitr::kable(sens_tbl, digits = c(0, 0, 2, 2, 2, 3, 0, 3, 0),
             caption = paste("Window-size sensitivity. GC mean is invariant;",
                             "the Tm distribution shifts and narrows; the",
                             "MutL-AR effect size varies within a few tenths",
                             "of a degree while the p value tracks the number",
                             "of windows."))

## ----sensitivity-tracks, fig.width=8, fig.height=8, fig.cap="Multi-scale integration. Concentric rings from outside in: MutL-AR peaks (ideogram), GC content at 50/100/200/500 bp (blues), Tm at 50/100/200/500 bp (reds), and microsatellite density (green). Yellow bands mark MutL-AR peaks."----
## scale_cols_gc / scale_cols_tm are defined in the density-panel chunk above
sens_dfs <- lapply(sens, function(x) as.data.frame(x$gr[, c("Tm", "GC")]))

tracks_scale <- c(
  list(list(type = "rect", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
            col = "#2C3E50", bg.col = "grey", name = "MutL-AR",
            legend_font_col = "#2C3E50", ideogram = TRUE, height = 0.5)),
  lapply(seq_along(window_sizes), function(i)
    list(type = "line", data = sens_dfs[[i]], value_col = "GC",
         name = paste0("GC ", window_sizes[i], " bp"),
         col = scale_cols_gc[i], legend_font_col = scale_cols_gc[i])),
  lapply(seq_along(window_sizes), function(i)
    list(type = "line", data = sens_dfs[[i]], value_col = "Tm",
         name = paste0("Tm ", window_sizes[i], " bp"),
         col = scale_cols_tm[i], legend_font_col = scale_cols_tm[i])),
  list(list(type = "line", data = ecoli_rep_hotspots$bins_rep,
            value_col = "count", name = "Microsatellites", col = "#2ECC71",
            legend_font_col = "#2ECC71"),
       list(type = "highlight",
            data = ecoli_rep_hotspots$all_peaks_IP_mutH,
            col = "#F1C40F", alpha = 0.18))
)

plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks_scale,
  circular    = TRUE,
  label       = label
)

## ----sensitivity-zoom, fig.width=10, fig.height=7.5, fig.cap="The 0.1-0.3 Mb region, within the interval shown in the genome-wide figure of the main case study, with Tm computed at 50, 100, 200 and 500 bp (dark to light). The 50 bp profile resolves single-window fluctuations that the 500 bp profile averages away, but all four trace the same underlying landscape. Yellow bands mark MutL-AR peaks."----
## Only the four Tm profiles are drawn. Ten line tracks in one linear panel
## leave each too little vertical space for its own axis labels, which then
## collide, and everything except Tm is answering a different question: the
## GC layers duplicate what the summary table already gives numerically, and
## the microsatellite density is the subject of the main case study's figure,
## not of this one. Removing both leaves five tracks with room to be read.
## The full ten-track version is the circular figure above.
##
## The interval is the one used in the main case study's zoomed panel, so the
## two figures can be read against each other rather than against different
## parts of the chromosome.
tracks_zoom <- c(
  list(list(type = "rect", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
            col = "#2C3E50", bg.col = "grey", name = "MutL-AR",
            legend_font_col = "#2C3E50", ideogram = TRUE, height = 0.6)),
  lapply(seq_along(window_sizes), function(i)
    list(type = "line", data = sens_dfs[[i]], value_col = "Tm",
         name = paste0("Tm ", window_sizes[i], " bp"),
         col = scale_cols_tm[i], legend_font_col = scale_cols_tm[i],
         height = 1.2)),
  list(list(type = "highlight",
            data = ecoli_rep_hotspots$all_peaks_IP_mutH,
            col = "#F1C40F", alpha = 0.18))
)

plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks_zoom,
  zoom        = "U00096.3:100000-300000",
  track.gap   = 0.03,
  axis.cex    = 0.55
)

## ----sensitivity-effect-sizes, echo=FALSE-------------------------------------
eff <- unique(test_summary[, c("window_bp", "Tm_med_diff")])
gc_rows <- test_summary[test_summary$target == "GC", ]

## ----session-info-------------------------------------------------------------
sessionInfo()

