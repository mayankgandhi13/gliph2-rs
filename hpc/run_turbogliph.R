# Runs turboGliph::gliph2 and writes its convergence groups in the same layout
# as `gliph2-rs gliph2` (type, tag, cluster_size, members).
# Usage: Rscript hpc/run_turbogliph.R INPUT.tsv REFERENCE.tsv OUT.tsv N_CORES PARAMS
#   PARAMS = "default" (turboGliph code defaults) or "paper" (GLIPH2 parameter files)

suppressPackageStartupMessages(library(turboGliph))
# turboGliph parallelises with doParallel workers; stringdist would otherwise
# start its own OpenMP pool in every worker and oversubscribe the node.
options(sd_num_thread = 1)
Sys.setenv(OMP_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
input <- read.delim(args[1], colClasses = "character")
reference <- read.delim(args[2], colClasses = "character")[, c("CDR3b", "TRBV")]
n_cores <- as.integer(args[4])

common <- list(
  cdr3_sequences = input[, c("CDR3b", "TRBV")], refdb_beta = reference, result_folder = "",
  motif_distance_cutoff = 3, kmer_mindepth = 3, structboundaries = TRUE, boundary_size = 3,
  motif_length = c(2, 3, 4), discontinuous_motifs = FALSE, global_vgene = FALSE,
  cluster_min_size = 2, n_cores = n_cores
)
params <- switch(args[5],
  default = list(lcminp = 0.01, lcminove = c(1000, 100, 10), min_seq_length = 0, all_aa_interchangeable = FALSE),
  paper = list(lcminp = 0.001, lcminove = 10, min_seq_length = 8, all_aa_interchangeable = TRUE),
  stop("PARAMS must be default or paper")
)

t0 <- Sys.time()
res <- do.call(gliph2, c(common, params))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cp <- res$cluster_properties
sorted_members <- vapply(strsplit(cp$members, " "), function(m) paste(sort(unique(m)), collapse = " "), "")
out <- data.frame(type = cp$type, tag = cp$tag, cluster_size = cp$cluster_size, members = sorted_members)
write.table(out, args[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf('{"tool":"turboGliph","params":"%s","input":%d,"groups":%d,"t_total":%.2f}\n',
            args[5], nrow(input), nrow(out), elapsed))
