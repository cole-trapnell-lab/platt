test_that("KA' and KA'' cell type notes are read from distinct directories", {
  expect_equal(
    ai_note_dirname("Kolmer-Agduhr neuron, spinal cord (KA')"),
    "kolmer-agduhr-neuron-spinal-cord-ka-prime"
  )
  expect_equal(
    ai_note_dirname("Kolmer-Agduhr neuron, spinal cord (KA'')"),
    "kolmer-agduhr-neuron-spinal-cord-ka-double-prime"
  )
  for (ct in c("epithelial cell, Kupffer's vesicle", "podocyte, pronephric glomerulus")) {
    expect_equal(ai_note_dirname(ct), blogdown:::dash_filename(ct))
  }
})
