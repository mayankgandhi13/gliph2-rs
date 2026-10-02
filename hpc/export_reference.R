# Exports each GLIPH reference repertoire in reference_list.RData to a TSV
# that gliph2-rs reads (CDR3b, TRBV columns).
# Usage: Rscript hpc/export_reference.R reference_list.RData OUT_DIR

args <- commandArgs(trailingOnly = TRUE)
env <- new.env()
load(args[1], envir = env)
str(mget(ls(env), envir = env), max.level = 2)

refs <- get(ls(env)[1], envir = env)
for (name in names(refs)) {
  ref <- refs[[name]]
  # Entries are either a data frame or a list whose first element is the
  # CDR3b/TRBV table (alongside V-gene and length frequency tables).
  tab <- if (is.data.frame(ref)) ref else ref[[1]]
  if (!is.data.frame(tab)) next
  cdr3 <- grep("^cdr3", names(tab), ignore.case = TRUE, value = TRUE)[1]
  v <- grep("^(trbv|v)", names(tab), ignore.case = TRUE, value = TRUE)[1]
  out <- data.frame(CDR3b = tab[[cdr3]], TRBV = if (is.na(v)) NA else tab[[v]])
  path <- file.path(args[2], paste0(name, ".tsv"))
  write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(name, nrow(out), "->", path, "\n")
}
