## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo    = TRUE,
  eval    = FALSE,   # the hg38 chunks need BSgenome.Hsapiens.UCSC.hg38 and minutes
  message = FALSE,
  warning = FALSE,
  fig.retina = 1,
  dpi        = 72
)

## ----install, eval=FALSE------------------------------------------------------
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

## ----quick-start--------------------------------------------------------------
# library(TmCalculator)
# hg38 <- "BSgenome.Hsapiens.UCSC.hg38"
# 
# tm <- tm_calculate(hg38, window = 200, slide = 200)$gr
# tm
# ## GRanges object with 14687412 ranges and 2 metadata columns:
# ##       seqnames      ranges strand |        Tm        GC

## ----quick-start-parallel-----------------------------------------------------
# tm <- tm_calculate(hg38, window = 200, slide = 200,
#                    BPPARAM = BiocParallel::SnowParam(workers = 5))$gr

## ----sources------------------------------------------------------------------
# tm_calculate(hg38, window = 200)         # an installed BSgenome package
# tm_calculate("contigs.fa.gz", window = 200)   # a FASTA file, gzip is fine
# tm_calculate(oligos)                     # a character vector of sequences

## ----regions------------------------------------------------------------------
# ## Chromosomes or records, by name or by number. On a BSgenome the chr
# ## prefix is added or removed to match the genome, so the same code works on
# ## UCSC and on Ensembl; FASTA record names are matched exactly.
# tm_calculate(hg38, regions = 1:22, window = 200)               # autosomes
# tm_calculate(hg38, regions = paste0("chr", c(1:22, "X", "Y")), # no chrM
#              window = 200)
# tm_calculate("contigs.fa.gz", regions = c("contig_7", "contig_9"),
#              window = 200)
# 
# ## Coordinate intervals. Commas and scientific notation are accepted, so a
# ## number pasted out of a genome browser works as it stands.
# tm_calculate(hg38, regions = c("chr1:1-10e6", "chrX:5,000,000-6,000,000"),
#              window = 200)
# 
# ## A mixture, when some chromosomes are wanted whole and others in part.
# tm_calculate(hg38, regions = c(1:20, "chr21:1-10e6", "X", "Y"), window = 200)
# 
# ## A GRanges, when the regions come from an annotation.
# prom <- promoters(genes(TxDb.Hsapiens.UCSC.hg38.knownGene),
#                   upstream = 1000, downstream = 500)
# tm_calculate(hg38, regions = prom, window = 50, slide = 25)
# 
# ## A GRanges carrying its own sequences is a source rather than a query, and
# ## then regions selects by overlap: a seqname takes every range on it, an
# ## interval takes the ranges it meets. Whole ranges come back, not clipped
# ## pieces of them, since their sequences are already fixed.
# tm_calculate(probes_gr, regions = "chr7")
# tm_calculate(probes_gr, regions = "chr7:1-1e6")
# 
# ## window = NULL, the default, gives one window per region, which is what
# ## short records call for: a FASTA of array probes, primers or synthetic
# ## oligos returns one Tm per record. It is refused for a region over 1 Mb,
# ## where a single Tm would mean nothing.
# tm_calculate("probes.fa", window = NULL)
# tm_calculate(oligos, BPPARAM = BiocParallel::SnowParam(5))

## ----unit---------------------------------------------------------------------
# tm_calculate(hg38, window = 200, unit = "segment",   # the default: 73 tasks
#              segment_size = 50e6, BPPARAM = BiocParallel::SnowParam(5))
# tm_calculate(hg38, window = 200, unit = "region",    # 24 tasks
#              BPPARAM = BiocParallel::SnowParam(5))

## ----sweep-data, eval=TRUE----------------------------------------------------
# Read the shipped summaries rather than transcribing numbers, so the table
# cannot drift from the measurements it describes.
read_sweep <- function(file, env) {
  d <- utils::read.csv(system.file("extdata", file, package = "TmCalculator"),
                       stringsAsFactors = FALSE)
  d$Environment <- env
  d
}
sweep <- rbind(read_sweep("bench_hg38_laptop.csv",  "Laptop"),
               read_sweep("bench_hg38_cluster.csv", "Compute node"))
sweep <- sweep[order(sweep$Environment != "Laptop", sweep$n_workers), ]

## ----sweep-table, eval=TRUE---------------------------------------------------
knitr::kable(
  data.frame(
    Environment = sweep$Environment,
    Workers     = sweep$n_workers,
    `Wall time (s)` = sprintf("%.1f [%.1f-%.1f]", sweep$wall_s, sweep$lo, sweep$hi),
    Speedup     = sprintf("%.2f", sweep$speedup),
    `Peak RSS per worker (GB)` = sprintf("%.2f", sweep$peak_worker_gb),
    check.names = FALSE),
  row.names = FALSE,
  caption = paste("Median of three repetitions, observed range in brackets.",
                  "Wall time includes worker start-up."))

## ----sweep-figure, eval=TRUE, fig.width=8, fig.height=3.6, fig.cap="Wall-clock time and peak memory per worker against worker count, on a six-core 16 GB laptop (solid, filled) and a compute node given six slots and 16 GB (dashed, open). Points are medians of three repetitions; bars give the observed range."----
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.2, 0.8), mgp = c(2.6, 0.7, 0))
for (what in c("wall", "rss")) {
  ys <- if (what == "wall") sweep$wall_s else sweep$peak_worker_gb
  plot(range(sweep$n_workers), range(0, ys * 1.05), type = "n",
       xlab = "Workers", xaxt = "n", adj = 0,
       ylab = if (what == "wall") "Wall clock (s)" else "Peak resident size per worker (GB)",
       main = if (what == "wall") "A" else "B")
  axis(1, at = sort(unique(sweep$n_workers)))
  for (e in unique(sweep$Environment)) {
    d <- sweep[sweep$Environment == e, ]
    y <- if (what == "wall") d$wall_s else d$peak_worker_gb
    solid <- e == "Laptop"
    if (what == "wall")
      arrows(d$n_workers, d$lo, d$n_workers, d$hi, angle = 90, code = 3,
             length = 0.03, col = "grey40")
    lines(d$n_workers, y, lty = if (solid) 1 else 2, col = "grey20")
    points(d$n_workers, y, pch = if (solid) 19 else 1, col = "grey20")
  }
  if (what == "wall")
    legend("topright", c("Laptop, 6 cores", "Compute node, 6 slots"),
           lty = c(1, 2), pch = c(19, 1), bty = "n", cex = 0.85, col = "grey20")
}
par(op)

## ----sweep-your-machine-------------------------------------------------------
# for (n in 1:6) {
#   t <- system.time(
#     tm_calculate(hg38, regions = "chr21", window = 200, slide = 200,
#                  method = "tm_nn", segment_size = 10e6,
#                  BPPARAM = BiocParallel::SnowParam(n), verbose = FALSE))
#   cat(sprintf("%d workers: %.1f s\n", n, t[["elapsed"]]))
# }

## ----downstream---------------------------------------------------------------
# tm_annot <- integrate_granges(gr_tm = tm, gr_features = atac_peaks,
#                               strategy = "overlap", weight = "overlap")
# compare_groups(tm_annot, value_cols = "Tm", group_col = "class")

## ----tm-calculate-------------------------------------------------------------
# tm_calculate(c("ACGTGCTAGCTAGCTAGC", "GGCCATATATGCGC"), method = "tm_nn", Na = 50)

## ----session-info, eval=TRUE--------------------------------------------------
sessionInfo()

