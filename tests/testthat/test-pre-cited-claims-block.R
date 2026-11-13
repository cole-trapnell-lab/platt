test_that("claims_for_cell_type returns the cell type's entry or empty", {
  claims <- c(notochord = "- n claim [@PMID1]", podocyte = NA)
  expect_equal(platt:::claims_for_cell_type(claims, "notochord"), "- n claim [@PMID1]")
  expect_equal(platt:::claims_for_cell_type(claims, "podocyte"), "")
  expect_equal(platt:::claims_for_cell_type(claims, "somite"), "")
  expect_equal(platt:::claims_for_cell_type(NULL, "notochord"), "")
  expect_equal(platt:::claims_for_cell_type(list(notochord = "x"), "notochord"), "x")
})

test_that("with_pre_cited_claims puts gene then cell type claims in one block", {
  out <- platt:::with_pre_cited_claims("CONTEXT", "GENE [@PMID1]", "CT [@PMID2]")
  expect_equal(out, "<pre_cited_claims>\nGENE [@PMID1]\n\nCT [@PMID2]\n</pre_cited_claims>\n\nCONTEXT")
  expect_equal(
    platt:::with_pre_cited_claims("CONTEXT", NULL, "CT [@PMID2]"),
    "<pre_cited_claims>\nCT [@PMID2]\n</pre_cited_claims>\n\nCONTEXT"
  )
  expect_equal(
    platt:::with_pre_cited_claims("CONTEXT", "GENE [@PMID1]", ""),
    "<pre_cited_claims>\nGENE [@PMID1]\n</pre_cited_claims>\n\nCONTEXT"
  )
  expect_equal(platt:::with_pre_cited_claims("CONTEXT", "", ""), "CONTEXT")
  expect_equal(platt:::with_pre_cited_claims("CONTEXT", NULL, ""), "CONTEXT")
})
