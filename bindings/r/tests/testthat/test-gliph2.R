reference <- local({
  aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
  # 400 distinct reference CDR3s; the first 300 carry the G%ETQ interior.
  vapply(seq_len(400), function(i) {
    x <- aa[(i - 1) %% 20 + 1]
    y <- aa[((i - 1) %/% 20) %% 20 + 1]
    if (i <= 300) paste0("C", x, y, "G", x, "ETQYQF") else paste0("CASSL", x, y, "GYEQYF")
  }, character(1))
})

test_that("global groups use the CDR3 interior and respect BLOSUM", {
  seqs <- c("CASGIETQYQF", "CASGVETQYQF", "CASGWETQYQF")
  strict <- gliph2(seqs, reference, local_similarities = FALSE)
  tags <- strict$cluster_properties$tag
  # I and V are BLOSUM62-compatible; W is not, so it stays out.
  expect_true("G%ETQ_IV" %in% tags)
  expect_false(any(grepl("W", tags, fixed = TRUE)))

  loose <- gliph2(seqs, reference,
    local_similarities = FALSE,
    all_aa_interchangeable = TRUE
  )
  expect_true("G%ETQ_IVW" %in% loose$cluster_properties$tag)
})

test_that("members are returned and match cluster_size", {
  res <- gliph2(c("CASGIETQYQF", "CASGVETQYQF"), reference, local_similarities = FALSE)
  cp <- res$cluster_properties
  expect_gt(nrow(cp), 0)
  expect_equal(lengths(res$cluster_list), cp$cluster_size, ignore_attr = TRUE)
  expect_setequal(res$cluster_list[["G%ETQ_IV"]], c("CASGIETQYQF", "CASGVETQYQF"))
})

test_that("enriched motifs give local groups", {
  # WWW appears in the interior of 4 sample CDR3s and never in the reference.
  seqs <- c(
    "CASWWWDEKGQYF", "CASWWWDEKDEKGQYF", "CASWWWDEKDEKDEKGQYF",
    "CASWWWQQQGQYF", "CASSLAAGYEQYF"
  )
  res <- gliph2(seqs, reference, lcminove = 1, global_similarities = FALSE)
  expect_true("WWW" %in% res$selected_motifs$motif)
  local <- res$cluster_properties
  expect_true(all(local$type == "local"))
  expect_true(any(grepl("^WWW_", local$tag)))
})

test_that("donor association separates case-only groups", {
  donors <- c(paste0("D", 1:6, ":CMV+"), paste0("D", 7:12, ":CMV-"))
  input <- data.frame(
    CDR3b = c(
      rep("CASGQETQYQF", 3), rep("CASGPETQYQF", 3), # 6 case donors
      rep("CASWAWWWYQF", 3), rep("CASWCWWWYQF", 3), # 3 case, 3 control
      rep("CASSLKKKYQF", 3)                         # the other 3 control donors
    ),
    patient = c(donors[1:6], donors[c(1:3, 7:9)], donors[10:12]),
    stringsAsFactors = FALSE
  )
  res <- gliph2(input, reference,
    local_similarities = FALSE, all_aa_interchangeable = TRUE,
    case = "CMV+", control = "CMV-", assoc_min_donors = 2
  )
  cp <- res$cluster_properties
  expect_equal(res$association$n_case, 6)
  expect_equal(res$association$n_control, 6)

  case_only <- cp[cp$tag == "G%ETQ_PQ", ]
  expect_equal(case_only$case_donors, 6)
  expect_equal(case_only$control_donors, 0)
  expect_equal(case_only$assoc.p, 1 / choose(12, 6)) # only table at least this extreme

  shared <- cp[cp$tag == "W%WWW_AC", ]
  expect_equal(c(shared$case_donors, shared$control_donors), c(3, 3))
  expect_gt(shared$assoc.p, 0.5)

  filtered <- gliph2(input, reference,
    local_similarities = FALSE, all_aa_interchangeable = TRUE,
    case = "CMV+", control = "CMV-", assoc_min_donors = 2, assoc_max_p = 0.01
  )
  expect_equal(filtered$cluster_properties$tag, "G%ETQ_PQ")
})

test_that("results do not depend on thread count", {
  seqs <- unique(c(
    vapply(1:200, function(i) sprintf("CASSL%sGQ%sETQYF", LETTERS[i %% 20 + 1], LETTERS[i %% 7 + 1]), character(1)),
    "CASGIETQYQF", "CASGVETQYQF"
  ))
  one <- gliph2(seqs, reference, n_cores = 1)$cluster_properties
  many <- gliph2(seqs, reference, n_cores = 4)$cluster_properties
  expect_identical(one, many)
})

test_that("bad input is rejected", {
  expect_error(gliph2(character(), reference), "no rows")
  expect_error(gliph2(c("XXXX", "1234"), reference), "no valid CDR3")
  expect_error(gliph2(c("CASSLGQETQYF"), reference, case = "CMV+"), "together")
  expect_error(gliph2(data.frame(x = 1), reference), "CDR3b")
})
