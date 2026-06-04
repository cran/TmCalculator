## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6
)

library(TmCalculator)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(IRanges)

## ----install, eval = FALSE----------------------------------------------------
# suppressMessages({
#   library(BiocManager)
#   BiocManager::install(c("Biostrings", "BSgenome", "BSgenome.Hsapiens.UCSC.hg38"))
# })
# install.packages("TmCalculator")
# pak::pkg_install("Gviz")
# pak::pkg_install("JunhuiLi1017/TmCalculator@dev")

## ----basic_usage, message=FALSE-----------------------------------------------
seqs <- c("ATGCGCGAAAT", "GCTAGCTAGCT")
names(seqs) <- c("chr1:1-11:+:seq1", "chr2:1-11:+:seq2")
result <- tm_calculate(seqs, method = "tm_wallace")
result$gr

## ----genomic_input, message=FALSE---------------------------------------------
coords <- c(
  "chr1:1898000-1898050:+:BSgenome.Hsapiens.UCSC.hg38:seq1",
  "chr2:2563000-2563050:-:BSgenome.Hsapiens.UCSC.hg38:seq2"
)
result <- tm_calculate(coords, method = "tm_nn")
result$gr

## ----genomic_input_options, message=FALSE-------------------------------------
result$options

## ----fasta_processing, message=FALSE------------------------------------------
fasta_file <- system.file("extdata", "example1.fasta", package = "TmCalculator")
result_fa <- tm_calculate(fasta_file, method = "tm_nn")
result_fa$gr

## ----complement_mismatch, message=FALSE---------------------------------------
seqs <- c("GCTAGCCGACAATGGGCAGATAGTAGAAA")
comp_seqs <- c("CGATCGGCTATTACCCGTCTATCATCTTT")
result <- tm_calculate(
  input_seq = seqs,
  complement_seq = comp_seqs,
  method = "tm_nn"
)
result$gr

## ----complement_dangling_end, message=FALSE-----------------------------------
seqs <- c("GCTAGCCGACAATGGGCAGATAGTAGAAA")
comp_seqs <- c("GATCGGCTATTACCCGTCTATCATCTTT")
result <- tm_calculate(
  input_seq = seqs,
  complement_seq = comp_seqs,
  method = "tm_nn",
  shift = -1
)
result$gr

## ----genome_tiling, eval = FALSE----------------------------------------------
# coords <- make_genomiccoord(
#   bsgenome    = BSgenome.Hsapiens.UCSC.hg38,
#   chromosomes = "chr1",
#   window      = 200L,
#   slide       = 50L,
#   trim_N      = "ends",
#   verbose     = FALSE
# )
# gr <- to_genomic_ranges_fast(
#   list(pkg_name = "BSgenome.Hsapiens.UCSC.hg38", seq = head(coords, 20))
# )
# tm_tiled <- tm_calculate(gr, method = "tm_nn")

## ----tm_viz_data, message=FALSE-----------------------------------------------
set.seed(123)

chr1_starts <- sort(sample(1:249250621, 100))
chr1_ends   <- chr1_starts + sample(50:200, 100, replace = TRUE)
chr2_starts <- sort(sample(1:243199373, 50))
chr2_ends   <- chr2_starts + sample(50:200, 50, replace = TRUE)

tm_results <- GRanges(
  seqnames = Rle(c(rep("chr1", 100), rep("chr2", 50))),
  ranges = IRanges(
    start = c(chr1_starts, chr2_starts),
    end   = c(chr1_ends, chr2_ends)
  ),
  strand = Rle(sample(c("+", "-"), 150, replace = TRUE)),
  Tm = c(runif(100, 60, 80), runif(50, 60, 80)),
  GC = c(runif(100, 30, 70), runif(50, 30, 70))
)
genome(seqinfo(tm_results)) <- "hg38"

gr_tm <- GRanges(
  seqnames = c("chr1", "chr2", "chr1", "chr2", "chr1"),
  ranges = IRanges(
    start = c(100, 200, 300, 400, 150),
    end   = c(150, 250, 350, 450, 200)
  ),
  Tm = c(65.5, 68.2, 70.1, 63.8, 72.0)
)

## ----plot_karyotype, eval=FALSE-----------------------------------------------
# track_list <- list(
#   list(type = "points", data = tm_results, value_col = "Tm",
#        col = "#e41a1c", name = "Tm (deg C)"),
#   list(type = "bars", data = tm_results, value_col = "GC",
#        col = "#377eb8", name = "GC (%)")
# )
# 
# plot_karyotype_genome(
#   genome_assembly = "hg38",
#   track_list      = track_list,
#   chromosomes     = c("chr1", "chr2"),
#   title_name      = "Tm and GC - hg38"
# )
# 
# plot_karyotype_genome(
#   genome_assembly = "hg38",
#   track_list = list(
#     list(type = "line", data = tm_results[seqnames(tm_results) == "chr1"],
#          value_col = "Tm", col = "#984ea3", name = "Tm (deg C)")
#   ),
#   chromosomes = "chr1",
#   zoom        = "chr1:50000000-150000000",
#   title_name  = "Tm - chr1 zoom"
# )

## ----plot_genome_tracks, eval=FALSE-------------------------------------------
# plot_tm_genome_tracks(
#   tm_results,
#   chromosome_to_plot = "chr1",
#   genome_assembly    = "hg38",
#   track_type         = "gradient",
#   zoom               = "chr1:10000000-80000000"
# )
# 
# p_regions <- plot_tm_genome_tracks(
#   tm_results,
#   chromosome_to_plot = "chr1",
#   genome_assembly    = "hg38",
#   x_axis             = "regions",
#   track_type         = "points",
#   color_by           = "Tm"
# )
# print(p_regions)

## ----plot_heatmap, eval=FALSE-------------------------------------------------
# plot_tm_heatmap(gr_tm, genome_assembly = "hg19", plot_type = "karyogram")
# plot_tm_heatmap(gr_tm, genome_assembly = "hg19", plot_type = "faceted")
# plot_tm_heatmap(gr_tm, genome_assembly = "hg19", x_axis = "regions")

## ----plot_linear, fig.height=5------------------------------------------------
plot_tm_linear(
  gr_tm,
  color_by     = "Tm",
  facet_by_chr = FALSE,
  add_line     = TRUE
)

## ----plot_interactive, eval=FALSE---------------------------------------------
# plot_tm_heatmap_interactive(gr_tm, genome_assembly = "hg19", plot_type = "karyogram")
# plot_tm_genome_tracks_interactive(
#   tm_results,
#   chromosome_to_plot = "chr1",
#   genome_assembly    = "hg38",
#   x_axis             = "regions"
# )
# plot_tm_karyotype_interactive(
#   tm_results,
#   genome_assembly = "hg38",
#   chromosomes     = c("chr1", "chr2")
# )

## ----launch_shiny, eval = FALSE-----------------------------------------------
# TmCalculatorShiny::TmCalculator_shiny()

