#' Convergence groups with the GLIPH2 algorithm
#'
#' Groups CDR3 beta sequences that are likely to share antigen specificity,
#' following the GLIPH2 method as implemented by `turboGliph::gliph2()`, but
#' computed in Rust. Local groups come from motifs enriched against a
#' reference repertoire; global groups from CDR3s whose interiors differ at
#' one position. Groups are reported separately, not merged.
#'
#' Results are deterministic: repeated calls, and any value of `n_cores`,
#' give identical groups.
#'
#' @param cdr3_sequences A character vector of CDR3b amino-acid sequences, or
#'   a data frame with a `CDR3b` column and optionally `TRBV`, `patient`
#'   (or `subject`) and `counts`.
#' @param refdb_beta The reference repertoire: a character vector of CDR3b
#'   sequences, or a data frame with a `CDR3b` column.
#' @param lcminp Maximum motif p-value (rounded to two significant figures
#'   before comparison, as GLIPH2 does).
#' @param lcminove Minimum fold enrichment, either one number or one per
#'   entry of `motif_length`.
#' @param kmer_mindepth Minimum motif occurrences in the sample.
#' @param motif_length Motif lengths to test.
#' @param discontinuous_motifs Also test motifs with one wildcard position.
#' @param motif_distance_cutoff Motif occurrences whose start positions differ
#'   by less than this are chained into one position range.
#' @param boundary_size Residues ignored at each end of the CDR3.
#' @param min_seq_length Minimum CDR3 length; raised to `2 * boundary_size + 1`.
#' @param accept_sequences_with_C_F_start_end Keep only CDR3s matching `^C.*F$`.
#' @param all_aa_interchangeable If `FALSE`, global edges need a BLOSUM62
#'   score of at least 0 between the two residues.
#' @param global_vgene Require the same TRBV gene for global edges.
#' @param cluster_min_size Drop groups smaller than this.
#' @param local_similarities,global_similarities Compute each group type.
#' @param min_subjects Keep only groups spanning at least this many donors.
#' @param global_max_p,global_max_q Keep only global groups at or below this
#'   p-value or Benjamini-Hochberg q-value. GLIPH2 applies no such filter, so
#'   both default to `NULL`.
#' @param case,control Donor conditions to compare (e.g. `"CMV+"` and
#'   `"CMV-"`), read from the text after the last `:` in the donor column.
#'   When both are given, every group gets a donor-level Fisher p-value
#'   (`assoc.p`) and q-value (`assoc.q`).
#' @param assoc_min_donors Only test groups seen in at least this many
#'   labelled donors.
#' @param assoc_max_p,assoc_max_q Keep only groups at or below this
#'   association p- or q-value.
#' @param assoc_permute_seed Shuffle the case/control labels among donors with
#'   this seed, as a calibration control.
#' @param n_cores Threads to use; `0` lets Rust choose.
#'
#' @return A list with
#'   \describe{
#'     \item{cluster_properties}{one row per convergence group: `type`, `tag`,
#'       `cluster_size`, `unique_cdr3_sample`, `unique_cdr3_ref`, `OvE`,
#'       `fisher.score`, `fdr.q`, `n_subjects`, the association columns when
#'       `case`/`control` are given, and `members` (space-separated CDR3s).}
#'     \item{selected_motifs}{the enriched motifs behind the local groups.}
#'     \item{cluster_list}{`members` split into a character vector per group,
#'       named by tag.}
#'     \item{association}{donor counts and the number of tested and
#'       significant groups, or `NULL`.}
#'     \item{sample_size,reference_size}{sequences that passed filtering.}
#'   }
#'
#' @section Differences from turboGliph:
#' Group membership matches `turboGliph::gliph2()` exactly (verified on
#' 2,000-100,000 CDR3 subsets of the Emerson 2017 cohort). The
#' simulation-based cluster scores (network size, CDR3 length, V gene, clonal
#' expansion, HLA) are not computed, so `total.score` and friends are absent.
#'
#' @examples
#' # Two CDR3s whose interiors differ only at one position, where isoleucine
#' # and valine are interchangeable under BLOSUM62.
#' cdr3 <- c("CASGIETQYQF", "CASGVETQYQF", "CASSLNTEAFF")
#' reference <- c("CASSLAPGATNEKLFF", "CASSFGREQYF", "CASSQEGTEAFF")
#' res <- gliph2(cdr3, reference, local_similarities = FALSE)
#' res$cluster_properties[, c("type", "tag", "cluster_size", "members")]
#'
#' \dontrun{
#' # Motif groups need a real reference repertoire (GLIPH2 ships several).
#' res <- gliph2(my_cdr3s, refdb_beta = ref_CD48)
#'
#' # CMV association across donors labelled "HIP00110:CMV-" and so on
#' res <- gliph2(cohort, ref_CD48, case = "CMV+", control = "CMV-",
#'               assoc_max_q = 0.05)
#' }
#' @export
gliph2 <- function(cdr3_sequences,
                   refdb_beta,
                   lcminp = 0.01,
                   lcminove = c(1000, 100, 10),
                   kmer_mindepth = 3,
                   motif_length = c(2, 3, 4),
                   discontinuous_motifs = FALSE,
                   motif_distance_cutoff = 3,
                   boundary_size = 3,
                   min_seq_length = 0,
                   accept_sequences_with_C_F_start_end = TRUE,
                   all_aa_interchangeable = FALSE,
                   global_vgene = FALSE,
                   cluster_min_size = 2,
                   local_similarities = TRUE,
                   global_similarities = TRUE,
                   min_subjects = 0,
                   global_max_p = NULL,
                   global_max_q = NULL,
                   case = NULL,
                   control = NULL,
                   assoc_min_donors = 10,
                   assoc_max_p = NULL,
                   assoc_max_q = NULL,
                   assoc_permute_seed = NULL,
                   n_cores = 0) {
  input <- .as_tcr_frame(cdr3_sequences, "cdr3_sequences")
  reference <- .as_tcr_frame(refdb_beta, "refdb_beta")$CDR3b
  reference <- .keep_amino_acid(reference, "refdb_beta")
  input <- input[.is_amino_acid(input$CDR3b), , drop = FALSE]
  if (!nrow(input)) {
    stop("`cdr3_sequences` has no valid CDR3 amino-acid sequences.", call. = FALSE)
  }
  if (!length(reference) && isTRUE(local_similarities)) {
    stop("`refdb_beta` is empty; motif enrichment needs a reference repertoire.", call. = FALSE)
  }
  if (is.null(case) != is.null(control)) {
    stop("`case` and `control` must be given together.", call. = FALSE)
  }

  params <- list(
    lcminp = as.numeric(lcminp),
    lcminove = as.numeric(lcminove),
    kmer_mindepth = as.numeric(kmer_mindepth),
    motif_length = as.numeric(motif_length),
    discontinuous_motifs = isTRUE(discontinuous_motifs),
    motif_distance_cutoff = as.numeric(motif_distance_cutoff),
    boundary_size = as.numeric(boundary_size),
    min_seq_length = as.numeric(min_seq_length),
    accept_sequences_with_C_F_start_end = isTRUE(accept_sequences_with_C_F_start_end),
    all_aa_interchangeable = isTRUE(all_aa_interchangeable),
    global_vgene = isTRUE(global_vgene),
    cluster_min_size = as.numeric(cluster_min_size),
    local_similarities = isTRUE(local_similarities),
    global_similarities = isTRUE(global_similarities),
    min_subjects = as.numeric(min_subjects),
    global_max_p = .num_or_null(global_max_p),
    global_max_q = .num_or_null(global_max_q),
    case = case,
    control = control,
    assoc_min_donors = as.numeric(assoc_min_donors),
    assoc_max_p = .num_or_null(assoc_max_p),
    assoc_max_q = .num_or_null(assoc_max_q),
    assoc_permute_seed = .num_or_null(assoc_permute_seed),
    n_cores = as.numeric(n_cores)
  )

  res <- .Call(
    "wrap__gliph2_impl",
    input$CDR3b, input$TRBV, input$patient, input$counts,
    reference, params,
    PACKAGE = "gliph2rs"
  )

  groups <- .as_frame(res$groups)
  tags <- if (nrow(groups)) groups$tag else character()
  list(
    cluster_properties = groups,
    selected_motifs = .as_frame(res$motifs),
    cluster_list = stats::setNames(strsplit(groups$members, " ", fixed = TRUE), tags),
    association = if (is.null(res$association)) NULL else as.data.frame(res$association),
    sample_size = res$n_sample,
    reference_size = res$n_reference,
    input_size = res$n_input_cdr3
  )
}

#' GLIPH2 parameters from the published parameter files
#'
#' The settings distributed with GLIPH2 (`local_min_pvalue = 0.001`,
#' `local_min_OVE = 10`, `cdr3_length_cutoff = 8`,
#' `all_aa_interchangeable = 1`), as a list to splice into [gliph2()] with
#' `do.call`. The defaults of [gliph2()] follow turboGliph's code instead.
#'
#' @return A named list of arguments.
#' @examples
#' str(gliph2_paper_params())
#' @export
gliph2_paper_params <- function() {
  list(
    lcminp = 0.001,
    lcminove = 10,
    min_seq_length = 8,
    all_aa_interchangeable = TRUE
  )
}

# Accepts a character vector or a data frame and returns a frame with the
# four columns the Rust side expects, using "" for missing values.
.as_tcr_frame <- function(x, arg) {
  if (is.character(x) || is.factor(x)) {
    x <- data.frame(CDR3b = as.character(x), stringsAsFactors = FALSE)
  }
  if (!is.data.frame(x)) {
    stop("`", arg, "` must be a character vector or a data frame.", call. = FALSE)
  }
  col <- function(names, default) {
    hit <- intersect(tolower(names), tolower(colnames(x)))
    if (!length(hit)) {
      return(rep(default, nrow(x)))
    }
    v <- x[[colnames(x)[match(hit[1], tolower(colnames(x)))]]]
    v <- as.character(v)
    v[is.na(v)] <- default
    v
  }
  cdr3 <- col(c("CDR3b", "cdr3", "cdr3aa", "CDR3.beta", "junction_aa"), "")
  if (!length(cdr3)) {
    stop("`", arg, "` has no rows.", call. = FALSE)
  }
  if (all(cdr3 == "")) {
    stop("`", arg, "` has no CDR3b column.", call. = FALSE)
  }
  counts <- suppressWarnings(as.numeric(col(c("counts", "count", "templates"), "1")))
  counts[is.na(counts)] <- 1
  data.frame(
    CDR3b = toupper(cdr3),
    TRBV = col(c("TRBV", "v", "v_gene", "v_call"), ""),
    patient = col(c("patient", "subject", "subject:condition", "donor", "sample_id"), ""),
    counts = counts,
    stringsAsFactors = FALSE
  )
}

# The 20 standard amino acids, matching the Rust parser.
.is_amino_acid <- function(x) grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", x)

.keep_amino_acid <- function(x, arg) {
  ok <- .is_amino_acid(x)
  if (any(!ok)) {
    warning(sum(!ok), " of ", length(x), " `", arg,
            "` sequences dropped: not amino-acid CDR3s.", call. = FALSE)
  }
  x[ok]
}

.as_frame <- function(cols) {
  as.data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)
}

.num_or_null <- function(x) if (is.null(x)) NULL else as.numeric(x)
