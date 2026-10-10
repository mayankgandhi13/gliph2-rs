# Compares gliph2rs::gliph2() with turboGliph::gliph2() in one R session, on
# the same input and reference, under both parameter sets.
# Usage: Rscript hpc/compare_r_binding.R INPUT.tsv REFERENCE.tsv N_CORES
# Prints one JSON line per parameter set.

suppressPackageStartupMessages({
  library(gliph2rs)
  library(turboGliph)
})
# turboGliph parallelises with doParallel workers; stringdist would otherwise
# start its own OpenMP pool in every worker and oversubscribe the node.
options(sd_num_thread = 1)
Sys.setenv(OMP_NUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
input <- read.delim(args[1], colClasses = "character")
reference <- read.delim(args[2], colClasses = "character")[, c("CDR3b", "TRBV")]
n_cores <- as.integer(args[3])

# A group is identified by its type, tag and sorted member set.
key <- function(type, tag, members) {
  paste(type, tag, vapply(members, function(m) paste(sort(unique(m)), collapse = ","), ""))
}

settings <- list(
  paper = list(lcminp = 0.001, lcminove = 10, min_seq_length = 8, all_aa_interchangeable = TRUE),
  default = list(lcminp = 0.01, lcminove = c(1000, 100, 10), min_seq_length = 0, all_aa_interchangeable = FALSE)
)

for (name in names(settings)) {
  p <- settings[[name]]

  t0 <- Sys.time()
  rs <- do.call(gliph2, c(list(
    cdr3_sequences = input[, c("CDR3b", "TRBV")], refdb_beta = reference, n_cores = n_cores
  ), p))
  t_rs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t0 <- Sys.time()
  tg <- do.call(turboGliph::gliph2, c(list(
    cdr3_sequences = input[, c("CDR3b", "TRBV")], refdb_beta = reference, result_folder = "",
    motif_distance_cutoff = 3, kmer_mindepth = 3, structboundaries = TRUE, boundary_size = 3,
    motif_length = c(2, 3, 4), discontinuous_motifs = FALSE, global_vgene = FALSE,
    cluster_min_size = 2, n_cores = n_cores
  ), p))
  t_tg <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  a <- key(tg$cluster_properties$type, tg$cluster_properties$tag,
           strsplit(tg$cluster_properties$members, " "))
  b <- key(rs$cluster_properties$type, rs$cluster_properties$tag, rs$cluster_list)

  cat(sprintf(
    '{"params":"%s","n":%d,"turbogliph_groups":%d,"gliph2rs_groups":%d,"identical":%s,"only_turbogliph":%d,"only_gliph2rs":%d,"t_gliph2rs":%.2f,"t_turbogliph":%.2f,"speedup":%.1f}\n',
    name, nrow(input), length(a), length(b),
    tolower(identical(sort(a), sort(b))),
    length(setdiff(a, b)), length(setdiff(b, a)), t_rs, t_tg, t_tg / t_rs
  ))
}
