## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo    = TRUE,
  eval    = FALSE,   # requires BSgenome.Hsapiens.UCSC.hg38 and long runtimes
  message = FALSE,
  warning = FALSE,
  fig.retina = 1,
  dpi        = 72
)

## ----libraries----------------------------------------------------------------
# library(TmCalculator)
# library(BiocParallel)
# 
# pkg <- "BSgenome.Hsapiens.UCSC.hg38"
# suppressPackageStartupMessages(library(pkg, character.only = TRUE))
# genome <- get(pkg, envir = asNamespace(pkg))

## ----measuring----------------------------------------------------------------
# gc(reset = TRUE)                      # reset the "max used" high-water mark
# ## ... run the code to be measured ...
# gc()                                  # read peak from the "max used" column
# 
# ps::ps_memory_info()[["rss"]] / 1e9   # current process resident set, GB

## ----chr21-warmup-------------------------------------------------------------
# chr_len21 <- GenomeInfoDb::seqlengths(genome)[["chr21"]]
# 
# t21 <- system.time({
#   bins21 <- make_genomiccoord(bsgenome = pkg, chromosomes = "chr21",
#                               window = 200L, slide = 200L,
#                               start = 1, end = chr_len21, strand = "+")
#   gr21   <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins21),
#                                    method = "preload_chr")
#   tm21   <- tm_calculate(gr21, method = "tm_nn",
#                          nn_table = "DNA_NN_SantaLucia_2004", Na = 50)
# })
# t21["elapsed"]        # ~10-15 s cold, ~6 s warm, on the test machine
# head(tm21$gr)

## ----chr1-serial--------------------------------------------------------------
# chr_len <- GenomeInfoDb::seqlengths(genome)[["chr1"]]
# 
# t_coord <- system.time({
#   bins <- make_genomiccoord(bsgenome = pkg, chromosomes = "chr1",
#                             window = 200L, slide = 200L,
#                             start = 1, end = chr_len, strand = "+")
# })
# 
# t_extract <- system.time({
#   gr_batch <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins),
#                                      method = "preload_chr")
# })
# 
# base::gc(reset = TRUE)
# t_tm <- system.time({
#   tm_chr1 <- tm_calculate(gr_batch, method = "tm_nn",
#                           nn_table = "DNA_NN_SantaLucia_2004", Na = 50)
# })
# base::gc()   # "max used" = peak R memory during the Tm step
# 
# rbind(t_coord, t_extract, t_tm)[, "elapsed"]
# ## t_coord t_extract      t_tm
# ##     5.5      15.4      30.8     (freshly booted, otherwise idle machine)

## ----chr1-parallel------------------------------------------------------------
# system.time({
#   tm_serial <- tm_calculate(gr_batch, method = "tm_nn", Na = 50)
# })
# ## elapsed ~ 30 s
# 
# ## 1.0.x only: BPPARAM no longer exists
# system.time({
#   tm_snow <- tm_calculate(gr_batch, method = "tm_nn", Na = 50,
#                           BPPARAM = SnowParam(workers = 5))
# })
# ## elapsed ~ 55 s on the idle machine (33-80 s across repeated sessions)
# ## -- NEVER faster than serial

## ----genome-parallel----------------------------------------------------------
# chrs <- paste0("chr", c(1:22, "X", "Y"))
# ## Better: sort largest-first for load balance -- see "Choosing the
# ## worker count" below for the one-liner.
# n_workers <- 5   # see "Choosing the worker count" below
# 
# runtime <- system.time({
#   res_list <- bplapply(chrs, function(chr, pkg) {
#     ## SnowParam workers are fresh R processes: load packages HERE,
#     ## not in the manager session.
#     suppressPackageStartupMessages(library(TmCalculator))
#     suppressPackageStartupMessages(library(pkg, character.only = TRUE))
# 
#     genome  <- get(pkg, envir = asNamespace(pkg))
#     chr_len <- GenomeInfoDb::seqlengths(genome)[[chr]]
# 
#     bins <- make_genomiccoord(bsgenome = pkg, chromosomes = chr,
#                               window = 200L, slide = 200L,
#                               start = 1, end = chr_len, strand = "+",
#                               verbose = FALSE)
#     gr  <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins),
#                                   method = "preload_chr")
#     out <- tm_calculate(gr, method = "tm_nn",
#                         nn_table = "DNA_NN_SantaLucia_2004",
#                         Na = 50)$gr        # serial inside the worker
# 
#     ## Drop sequence columns before returning: Tm/GC are what we keep,
#     ## and this cuts per-chromosome serialization from ~500 MB to a few MB.
#     out$sequence   <- NULL
#     out$complement <- NULL
# 
#     ## Record this worker's peak resident memory (GB) for budgeting.
#     attr(out, "worker_rss_gb") <- ps::ps_memory_info()[["rss"]] / 1e9
#     out
#   }, pkg = pkg, BPPARAM = SnowParam(workers = n_workers))
# 
#   tm_genome <- unlist(GenomicRanges::GRangesList(res_list))
# })
# 
# runtime
# sapply(res_list, attr, "worker_rss_gb")   # per-worker memory check
# length(tm_genome)                         # 14,687,330 windows
# summary(tm_genome$Tm)
# ## Min. 47.2  1st Qu. 68.9  Median 71.5  Mean 72.0  3rd Qu. 74.8  Max. 101.1

## ----strategies---------------------------------------------------------------
# sl   <- GenomeInfoDb::seqlengths(genome)[paste0("chr", c(1:22, "X", "Y"))]
# chrs <- names(sort(sl, decreasing = TRUE))          # largest first
# 
# tasks_chrom <- lapply(chrs, function(ch)
#   list(chr = ch, start = 1L, end = as.integer(sl[[ch]])))
# 
# ## Segment boundaries must be multiples of `slide`, or the window grid
# ## shifts between segments and the result stops matching a whole-
# ## chromosome run. 50 Mb = 250,000 x 200 bp.
# seg_size  <- 50e6
# tasks_seg <- unlist(lapply(names(sl), function(ch) {
#   st <- seq(1, sl[[ch]], by = seg_size)
#   lapply(st, function(s)
#     list(chr = ch, start = as.integer(s),
#          end = as.integer(min(s + seg_size - 1, sl[[ch]]))))
# }), recursive = FALSE)
# 
# length(tasks_chrom)   # 24
# length(tasks_seg)     # 73
# 
# n <- 5
# BPPARAM_static  <- SnowParam(workers = n)
# BPPARAM_dynamic <- SnowParam(workers = n, tasks = length(tasks_chrom))
# BPPARAM_segment <- SnowParam(workers = n, tasks = length(tasks_seg))

## ----sweep-data, eval=TRUE----------------------------------------------------
# The two sweeps ship with the package as the CSV files the benchmark
# scripts wrote, and everything below is derived from them. Nothing is
# transcribed, so the tables, the figure and the numbers quoted in the text
# cannot drift apart, and re-running a sweep means replacing a file.
cols <- c("env", "strategy", "n_workers", "n_tasks", "rep", "wall_s",
          "startup_s", "wall_compute_s", "work_s", "hard_floor", "max_rss_gb")
read_sweep <- function(file, env) {
  d <- utils::read.csv(system.file("extdata", file, package = "TmCalculator"),
                       stringsAsFactors = FALSE)
  stopifnot(length(unique(d$n_windows)) == 1L)   # identical windows throughout
  d$env <- env
  d[, cols]
}
raw <- rbind(read_sweep("bench_parallel_strategy.csv", "laptop"),
             read_sweep("bench_parallel_cluster.csv",  "node"))

# Median over repetitions; the observed range is kept for the wall clock,
# which is the quantity the range bars in the figure show.
sweep <- do.call(rbind, lapply(
  split(raw, list(raw$env, raw$strategy, raw$n_workers), drop = TRUE),
  function(d) data.frame(
    env = d$env[1], strategy = d$strategy[1], workers = d$n_workers[1],
    tasks = d$n_tasks[1], n_rep = nrow(d),
    wall_s = median(d$wall_s), wall_lo = min(d$wall_s), wall_hi = max(d$wall_s),
    startup_s = median(d$startup_s), compute_s = median(d$wall_compute_s),
    work_s = median(d$work_s), longest_task_s = median(d$hard_floor),
    rss_gb = median(d$max_rss_gb), stringsAsFactors = FALSE)))

# Derived rather than transcribed. Speedup is the one-worker wall clock of
# the SAME strategy in the SAME environment over the configuration's wall
# clock: segmenting changes how much total work there is, so a ratio formed
# against a configuration's own summed task times would credit a strategy
# for doing more work and the three could not be compared.
key    <- paste(sweep$env, sweep$strategy)
ser    <- sweep[sweep$workers == 1L, ]
serial      <- stats::setNames(ser$wall_s, paste(ser$env, ser$strategy))
serial_work <- stats::setNames(ser$work_s, paste(ser$env, ser$strategy))
sweep$speedup    <- serial[key] / sweep$wall_s
sweep$efficiency <- sweep$speedup / sweep$workers
sweep$work_ratio <- sweep$work_s / serial_work[key]
sweep$strategy <- factor(sweep$strategy, levels = c("static", "dynamic", "segment"))
sweep$env      <- factor(sweep$env, levels = c("laptop", "node"))
sweep <- sweep[order(sweep$env, sweep$strategy, sweep$workers), ]
rownames(sweep) <- NULL

## ----sweep-table, eval=TRUE---------------------------------------------------
sweep_table <- function(s, caption) {
  tab <- data.frame(
    Strategy   = ifelse(duplicated(s$strategy), "", as.character(s$strategy)),
    Workers    = s$workers,
    Tasks      = s$tasks,
    `Wall (s)` = sprintf("%.1f [%.1f-%.1f]", s$wall_s, s$wall_lo, s$wall_hi),
    `Start-up (s)` = round(s$startup_s, 1),
    `Total task time (s)` = round(s$work_s, 1),
    `Work ratio` = round(s$work_ratio, 2),
    `Longest task (s)` = round(s$longest_task_s, 1),
    Speedup    = round(s$speedup, 2),
    Efficiency = round(s$efficiency, 2),
    `Peak RSS per worker (GB)` = round(s$rss_gb, 2),
    check.names = FALSE, stringsAsFactors = FALSE)
  knitr::kable(tab, row.names = FALSE, caption = caption)
}
sweep_table(sweep[sweep$env == "laptop", ],
  paste("Laptop (6-core Intel Core i7, 16 GB): task-partitioning strategy",
        "and worker count. Wall time is the median of three repetitions with",
        "the observed range in brackets."))

## ----sweep-table-node, eval=TRUE----------------------------------------------
sweep_table(sweep[sweep$env == "node", ],
  paste("Compute node (two 20-core Intel Xeon Gold 6230, six slots, 48 GB",
        "limit, shared with other jobs): the same sweep."))

## ----sweep-figure, eval=TRUE, fig.width=10.5, fig.height=3.9, fig.cap="Parallel performance under three task-partitioning strategies in two environments. Colour is the strategy; solid lines with filled points are the laptop, dashed lines with open points the compute node. (A) Speedup against worker count, relative to the one-worker run of the same strategy and environment; the dashed grey line marks linear speedup. (B) Total task time, summed over all tasks and measured inside the workers; values above the serial reference indicate that a configuration performed more work than the serial run rather than dividing the same work among more processes. (C) Peak resident set size of the heaviest worker."----
# Three panels, one argument, now made twice. Panel A shows that on the
# laptop every strategy turns over before the cores run out, which invites
# the usual explanation of a ragged schedule. Panel B rules that out: the
# total task time grows with the worker count, so the later configurations
# are doing more work, not dividing a fixed amount badly. Panel C gives the
# reason and the remedy together. The node, overlaid, shows what happens to
# each of the three when memory stops being scarce.
pal  <- c(static = "#1B5E9C", dynamic = "#C0392B", segment = "#5C6B73")
ltys <- c(laptop = "solid", node = "22")   # "22": a short dash, legible in a key
pchs <- c(laptop = 16, node = 1)
lv   <- levels(sweep$strategy)
ev   <- levels(sweep$env)
wk   <- sort(unique(sweep$workers))

op <- par(mfrow = c(1, 3), mar = c(4.4, 4.5, 2.2, 0.8), las = 1,
          cex = 0.8, mgp = c(2.8, 0.7, 0))

series <- function(col, ylab, ref = NULL, ref_lab = NULL, ymax = NULL,
                   diagonal = FALSE) {
  plot(NA, xlim = range(wk),
       ylim = c(0, if (is.null(ymax)) max(sweep[[col]]) * 1.05 else ymax),
       bty = "n", xaxt = "n", xlab = "Workers", ylab = ylab)
  axis(1, at = wk)
  if (diagonal) {
    abline(a = 0, b = 1, lty = 2, col = "grey55")
    text(max(wk), max(wk), "linear", adj = c(1.1, -0.4), cex = 0.75,
         col = "grey40")
  }
  if (!is.null(ref)) {
    abline(h = ref, lty = 2, col = "grey55")
    text(max(wk), ref, ref_lab, adj = c(1.1, -0.5), cex = 0.75, col = "grey40")
  }
  for (e in ev) for (st in lv) {
    d <- sweep[sweep$env == e & sweep$strategy == st, ]
    lines(d$workers, d[[col]], col = pal[st], lwd = 2, lty = ltys[e])
    points(d$workers, d[[col]], col = pal[st], pch = pchs[e], lwd = 1.6)
  }
}

series("speedup", "Speedup", ymax = max(wk) * 1.02, diagonal = TRUE)
mtext("A", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)
# One legend, two columns: environment keys on the left, strategy colours
# on the right. legend() fills column-major, so the shorter column is padded
# with a blank entry whose lty, pch and col are NA.
legend("topleft", bty = "n", cex = 0.8, ncol = 2, seg.len = 2.6,
       legend = c(ev, "", lv),
       lty = c(unname(ltys[ev]), NA, rep("solid", length(lv))),
       pch = c(unname(pchs[ev]), NA, rep(NA, length(lv))),
       col = c("grey25", "grey25", NA, unname(pal[lv])), lwd = 2)

series("work_s", "Total task time (s)",
       ref = median(sweep$work_s[sweep$workers == 1]), ref_lab = "serial")
mtext("B", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)

series("rss_gb", "Peak resident set size per worker (GB)")
mtext("C", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)

par(op)

## ----per-chrom, eval=TRUE-----------------------------------------------------
tk <- utils::read.csv(system.file("extdata", "bench_parallel_strategy_tasks.csv.gz",
                                  package = "TmCalculator"))
tk <- tk[tk$strategy == "dynamic", ]        # one task per chromosome

## Median over repetitions, one column per worker count.
w <- stats::reshape(
  stats::aggregate(secs ~ chr + n_workers, data = tk, FUN = stats::median),
  idvar = "chr", timevar = "n_workers", direction = "wide")
names(w) <- sub("^secs\\.", "w", names(w))

w$inflation <- w$w6 / w$w1                  # 6 workers vs serial
w <- w[order(-w$w1), ]
knitr::kable(w, row.names = FALSE, digits = 1,
             caption = "Laptop, dynamic dispatch: median seconds per chromosome at each worker count.")

## ----segment-parallel---------------------------------------------------------
# chr_len <- GenomeInfoDb::seqlengths(genome)[["chr1"]]
# 
# ## Segment length must be a multiple of `slide` so the window grid stays
# ## aligned across segment boundaries (50 Mb = 250,000 x 200).
# seg_starts <- seq(1, chr_len, by = 50e6)
# seg <- data.frame(start = seg_starts,
#                   end   = pmin(seg_starts + 50e6 - 1, chr_len))
# 
# t_seg <- system.time({
#   res_seg <- bplapply(seq_len(nrow(seg)), function(i, seg, pkg) {
#     suppressPackageStartupMessages(library(TmCalculator))
#     suppressPackageStartupMessages(library(pkg, character.only = TRUE))
#     bins <- make_genomiccoord(bsgenome = pkg, chromosomes = "chr1",
#                               window = 200L, slide = 200L,
#                               start = seg$start[i], end = seg$end[i],
#                               strand = "+", trim_N = "none",
#                               verbose = FALSE)
#     gr <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins),
#                                  method = "preload_chr")
#     out <- tm_calculate(gr, method = "tm_nn",
#                         nn_table = "DNA_NN_SantaLucia_2004", Na = 50)$gr
#     out$sequence <- NULL; out$complement <- NULL
#     out
#   }, seg = seg, pkg = pkg, BPPARAM = SnowParam(workers = 5))
# 
#   tm_chr1_seg <- sort(unlist(GenomicRanges::GRangesList(res_seg)))
# })
# 
# t_seg["elapsed"]      # 38 s measured, vs ~52 s for the serial FULL pipeline
# length(tm_chr1_seg)   # 1,152,300

## ----workers------------------------------------------------------------------
# mem_gb    <- 16                                     # your machine
# cores     <- parallel::detectCores(logical = FALSE) # physical cores
# per_worker_gb <- 2                                  # 50 Mb segment task
# n_workers <- min(cores - 1L, floor((mem_gb - 4) / per_worker_gb))
# n_workers
# ## 16 GB / 6 cores  -> 5 workers, which is what the sweep found fastest
# ## 32 GB / 8 cores  -> 7 workers
# ## 64 GB / 10 cores -> 9 workers

## ----sort-chrs----------------------------------------------------------------
# sl   <- GenomeInfoDb::seqlengths(genome)
# chrs <- names(sort(sl[paste0("chr", c(1:22, "X", "Y"))], decreasing = TRUE))

## ----sessioninfo, eval=TRUE---------------------------------------------------
sessionInfo()

