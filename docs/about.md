
The function `fit_genotype_ccm()`
* `genotype`
* `ccs`
* `ctrl_ids`
* `perturbation_col`
```

lmx1b_cds = load_monocle_objects("~/OneDrive/UW/Trapnell/hooke_manuscript/supplement/data_processing/prdm_projected_comb_cds/")

lmx1b_ccs = new_cell_count_set(lmx1b_cds, sample_group = "embryo", cell_group = "cell_type")

lmx1b_ccm = fit_genotype_ccm("lmx1b", 
                            lmx1b_ccs, 
                            ctrl_ids = c("control"), 
                            perturbation_col = "condition")

```


```

make_contrast(lmx1b_ccm)

```