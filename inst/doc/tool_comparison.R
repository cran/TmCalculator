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
  fig.width = 9,
  fig.height = 4.2
)

## ----reproduce, eval=FALSE----------------------------------------------------
# ## 1. Confirm the three tools agree before timing anything.
# ##    If they do not, the benchmark measures default-value differences.
# system("Rscript inst/scripts/bench_crosstool.R --calibrate --outdir bench200")
# 
# ## 2. rmelting accepts one sequence per call, so it is given its own,
# ##    smaller set of sizes; results merge into the same CSV.
# system(paste("Rscript inst/scripts/bench_crosstool.R",
#              "--tools rmelting --sizes 100,1000 --reps 3 --outdir bench200"))
# 
# system(paste("Rscript inst/scripts/bench_crosstool.R",
#              "--tools TmCalculator,Biopython",
#              "--sizes 1000,10000,100000 --reps 3 --outdir bench200"))
# 
# ## 3. Summary table and figure.
# system("Rscript inst/scripts/plot_crosstool.R --indir bench200")

## ----load---------------------------------------------------------------------
## The measurements are shipped with the package rather than recomputed here,
## since reproducing them needs a Java runtime, rmelting and a Python
## interpreter with Biopython, none of which is a dependency. If the file is
## absent -- a source tree in which the benchmark has not been run -- the
## remaining sections are skipped rather than failing the build.
f <- system.file("extdata", "crosstool_bench.csv", package = "TmCalculator")
has_data <- nzchar(f) && file.exists(f) &&
            length(readLines(f, n = 1L, warn = FALSE)) > 0L

## ----nodata, eval=!has_data, echo=FALSE, results='asis'-----------------------
# ## message() is swallowed by the global message = FALSE, which would leave a
# ## build with no data looking simply empty. Write the explanation into the
# ## document itself instead.
# cat("> **The benchmark results are not present in this build.**\n>\n",
#     "> `system.file(\"extdata\", \"crosstool_bench.csv\")` returned nothing,\n",
#     "> so the tables and figures below are omitted. To populate them, run\n",
#     "> `inst/scripts/bench_crosstool.R`, copy `bench_crosstool.csv` and\n",
#     "> `consistency.csv` into `inst/extdata/` as `crosstool_bench.csv` and\n",
#     "> `crosstool_consistency.csv`, and **reinstall the package** --\n",
#     "> `system.file()` resolves against the installed copy, not the source\n",
#     "> tree.\n", sep = "")

## ----load2, eval=has_data-----------------------------------------------------
S <- utils::read.csv(f, stringsAsFactors = FALSE)
S <- S[S$ok & !is.na(S$compute_s), , drop = FALSE]

## The benchmark file accumulates across invocations, so it can hold rows
## taken before and after a change to the package. Averaging those together
## gives a table that describes no build at all, so refuse instead.
if ("pkg_version" %in% names(S)) {
  vs <- sort(unique(stats::na.omit(S$pkg_version)))
  if (length(vs) > 1L)
    stop("crosstool_bench.csv mixes TmCalculator versions (",
         paste(vs, collapse = ", "), "); re-run the benchmark with --fresh.")
} else {
  stop("crosstool_bench.csv predates the split of start-up from file I/O. ",
       "Re-run inst/scripts/bench_crosstool.R and copy the result into ",
       "inst/extdata/, then reinstall.")
}

## Median and range over repetitions. With three repetitions a standard
## deviation is not a meaningful estimate, and repeated runs of the same
## configuration on a laptop have differed by nearly 30%, so the spread is
## reported as the range actually observed.
agg <- do.call(rbind, lapply(split(S, list(S$tool, S$n), drop = TRUE), function(d) {
  data.frame(tool = d$tool[1], n = d$n[1], reps = nrow(d),
             compute_med = stats::median(d$compute_s),
             compute_lo  = min(d$compute_s),
             compute_hi  = max(d$compute_s),
             seqs_per_s  = d$n[1] / stats::median(d$compute_s),
             rss_med     = stats::median(d$rss_gb),
             rss_lo      = min(d$rss_gb), rss_hi = max(d$rss_gb),
             startup_med = stats::median(d$startup_s),
             ## Reading the input and writing the results exist only because
             ## each tool runs as a separate process, so that peak memory can
             ## be attributed to it. Reported, never added to start-up, never
             ## drawn: a user passes data in memory.
             io_med      = stats::median(d$io_s),
             stringsAsFactors = FALSE)
}))
agg <- agg[order(agg$tool, agg$n), ]

## Start-up cannot depend on the input. If it does, something that scales
## with n is being counted as start-up.
for (t in unique(agg$tool)) {
  d <- agg[agg$tool == t, ]
  if (nrow(d) > 1L && max(d$startup_med) / min(d$startup_med) > 1.5)
    warning(t, " start-up varies ",
            sprintf("%.1fx", max(d$startup_med) / min(d$startup_med)),
            " across input sizes; it should be constant.", call. = FALSE)
}
knitr::kable(agg, digits = 4, row.names = FALSE,
             caption = "Median and range over repetitions at each input size.")

## ----consistency, eval=has_data-----------------------------------------------
fc <- system.file("extdata", "crosstool_consistency.csv", package = "TmCalculator")
if (nzchar(fc)) {
  knitr::kable(utils::read.csv(fc), digits = 10, row.names = FALSE,
               caption = "Deviation in Tm relative to TmCalculator on identical input.")
}

## ----figure, eval=has_data, fig.cap="Compute time (A) and peak resident set size (B) against input size. Panel A is linear on both axes: above a few thousand sequences each tool is proportional to its input, so each is a straight line and the ratio of the slopes is the ratio of the throughputs, read directly off the picture. On logarithmic axes that ratio would carry the same visual weight as every other ratio in the plot, including the reversed one at the smallest input. Panel A omits any tool not measured across the whole range; panel B keeps all of them, the spread there being small enough to show."----
tools <- sort(unique(agg$tool))
pal <- c("#1B5E9C", "#C0392B", "#5C6B73")[seq_along(tools)]; names(pal) <- tools
pch <- c(16, 17, 15)[seq_along(tools)];                      names(pch) <- tools

draw_range <- function(x, lo, hi, col) {
  v <- is.finite(lo) & is.finite(hi) & hi / pmax(lo, 1e-12) > 1.02
  if (any(v)) arrows(x[v], lo[v], x[v], hi[v], code = 3, angle = 90,
                     length = 0.03, col = col)
}

op <- par(mfrow = c(1, 2), mar = c(4.3, 4.4, 2.2, 0.8), las = 1, cex = 0.85)
xr <- range(agg$n)

## A tool measured at only the smallest sizes is left off the time axis: its
## cost is orders of magnitude larger, which would flatten the others onto the
## axis, and drawing two points as a line invites interpolation through sizes
## at which it was never run.
tt   <- sort(unique(agg$tool[agg$n == max(agg$n)]))
aggt <- agg[agg$tool %in% tt, ]

plot(NA, xlim = c(0, max(aggt$n)), ylim = c(0, max(aggt$compute_hi) * 1.05),
     bty = "n", xlab = "Sequences", ylab = "Compute time (s)", xaxt = "n")
axis(1, at = pretty(c(0, max(aggt$n))),
     labels = format(pretty(c(0, max(aggt$n))), big.mark = ",",
                     scientific = FALSE, trim = TRUE))
mtext("A", side = 3, adj = 0, font = 2, line = 0.8, cex = 1.1)
for (t in tt) {
  d <- aggt[aggt$tool == t, ]; d <- d[order(d$n), ]
  lines(d$n, d$compute_med, col = pal[t], lwd = 2)
  draw_range(d$n, d$compute_lo, d$compute_hi, pal[t])
  points(d$n, d$compute_med, col = pal[t], pch = pch[t], cex = 1.05)
}
legend("topleft", bty = "n", legend = tt, col = pal[tt],
       pch = pch[tt], lwd = 2, cex = 0.9)

## Expand to whole decades: range() alone puts the extreme points on the frame
yr2 <- range(c(agg$rss_lo, agg$rss_hi))
yr2 <- c(10^floor(log10(yr2[1])), 10^ceiling(log10(yr2[2])))
plot(NA, xlim = xr, ylim = yr2, log = "xy", bty = "n", yaxt = "n",
     xlab = "Sequences", ylab = "Peak resident set size (GB)")
ticks <- 10^seq(log10(yr2[1]), log10(yr2[2]))
axis(2, at = ticks, labels = format(ticks, scientific = FALSE, drop0trailing = TRUE))
axis(2, at = as.numeric(outer(2:9, ticks)), labels = FALSE, tcl = -0.2)
mtext("B", side = 3, adj = 0, font = 2, line = 0.8, cex = 1.1)
for (t in tools) {
  d <- agg[agg$tool == t, ]; d <- d[order(d$n), ]
  if (nrow(d) > 1) lines(d$n, d$rss_med, col = pal[t], lwd = 2)
  draw_range(d$n, d$rss_lo, d$rss_hi, pal[t])
  points(d$n, d$rss_med, col = pal[t], pch = pch[t], cex = 1.1)
}
par(op)

## ----decompose, eval=has_data-------------------------------------------------
## Two points at the top of the range separate the constant part of a call
## from the part that scales with the input.
fit <- do.call(rbind, lapply(split(agg, agg$tool), function(d) {
  d <- d[order(d$n), ]
  if (nrow(d) < 2L) return(NULL)
  k <- nrow(d)
  slope <- (d$compute_med[k] - d$compute_med[k - 1L]) / (d$n[k] - d$n[k - 1L])
  data.frame(tool = d$tool[1],
             fixed_cost_s = d$compute_med[k] - slope * d$n[k],
             marginal_seqs_per_s = 1 / slope, stringsAsFactors = FALSE)
}))
knitr::kable(fit, digits = 3, row.names = FALSE,
             caption = "Fixed cost per call and marginal throughput.")

## The size at which two tools take equal time follows from those two
## numbers, and is the one figure that keeps either end of panel A from
## being over-read.
if (nrow(fit) >= 2L) {
  o <- order(fit$marginal_seqs_per_s, decreasing = TRUE)
  a <- fit[o[1], ]; b <- fit[o[2], ]
  crossover <- (a$fixed_cost_s - b$fixed_cost_s) /
               (1 / b$marginal_seqs_per_s - 1 / a$marginal_seqs_per_s)
  cat(sprintf("%s overtakes %s at about %s sequences.\n",
              a$tool, b$tool, signif(crossover, 2)))
}

## ----walltime, eval=has_data, fig.width=7.5, fig.height=4.6, fig.cap="Elapsed time for one invocation, separated into start-up and compute. Start-up is the constant part: an R session with the Bioconductor packages attached costs about two seconds whatever the input, which is more than an entire Biopython run at the smaller sizes. Only tools measured across the whole range are shown."----
wall <- aggt[order(aggt$n, aggt$tool), ]
m <- rbind(compute = wall$compute_med, `start-up` = wall$startup_med)
colnames(m) <- paste(wall$tool, wall$n)
grp   <- cumsum(c(1, diff(wall$n) != 0))
space <- ifelse(c(TRUE, diff(grp) != 0), 1.1, 0.18)
fill  <- c("#1B5E9C", "#BFD3E6")

op <- par(mar = c(5.6, 4.6, 2.4, 0.8), las = 1, cex = 0.85)
bp <- barplot(m, space = space, col = fill, border = NA, ylab = "",
              names.arg = rep("", ncol(m)), ylim = c(0, max(colSums(m)) * 1.18))
u <- par("usr")
text(bp, u[3] - 0.02 * (u[4] - u[3]), labels = wall$tool, srt = 40, adj = 1,
     xpd = TRUE, cex = 0.72)
mtext(format(unique(wall$n), big.mark = ","), side = 1, line = 3.6,
      at = tapply(bp, grp, mean), cex = 0.85)
mtext("Sequences", side = 1, line = 4.7, cex = 0.85)
mtext("Elapsed time (s)", side = 2, line = 3.2, las = 0, cex = 0.85)
legend("topleft", bty = "n", fill = fill, border = NA,
       legend = c("compute", "start-up"), cex = 0.9)

## Where the TOTAL, not the compute time, crosses over.
if (nrow(fit) >= 2L) {
  o <- order(fit$marginal_seqs_per_s, decreasing = TRUE)
  a <- fit[o[1], ]; b <- fit[o[2], ]
  su <- vapply(list(a, b), function(z)
    stats::median(agg$startup_med[agg$tool == z$tool]), numeric(1))
  n_eq <- ((a$fixed_cost_s + su[1]) - (b$fixed_cost_s + su[2])) /
          (1 / b$marginal_seqs_per_s - 1 / a$marginal_seqs_per_s)
  cat(sprintf("Including start-up, %s overtakes %s at about %s sequences.\n",
              a$tool, b$tool, format(signif(n_eq, 2), big.mark = ",")))
}
par(op)

## ----batch, eval=FALSE--------------------------------------------------------
# rmelting::melting(sequence = c("ACGTACGTACGTACGTACGTACGT",
#                                "GGCCGGCCGGCCGGCCGGCCGGCC"),
#                   nucleic.acid.conc = 1.25e-8, hybridisation.type = "dnadna",
#                   Na.conc = 0.05, method.nn = "san04", correction.ion = "san96")
# #> Error: 'sequence' should be a character vector of length 1.

## ----session-info-------------------------------------------------------------
sessionInfo()

