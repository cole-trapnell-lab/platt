# Running DEGs over a graph

![](assets/degs_over_graph.png)


See explanation of gene patterns [here](https://cole-trapnell-lab.github.io/platt/patterns/): 

The function `compare_genes_within_state_graph()`:

* `ccs`- a Hooke `cell_count_set` object
* `graph`
* `gene_ids`
* `cores`

```
pf_graph_degs = compare_genes_over_graph(pf_ccs,
                                         pf_graph)
```


# Running DEGs within each perturbation

```
chem10_cds = load_monocle_objects("/net/trapnell/vol1/home/elizab9/projects/projects/CHEMFISH/manuscript/data/chem10_projected_comb_cds_v2.0.2_remove_outliers")

chem10_ccs = new_cell_count_set(chem10_cds,
                         sample_group = "embryo",
                         cell_group = "cell_type")
```

The function `compare_genes_within_state_graph()`: 

* `ccs`- a Hooke `cell_count_set` object
* `perturbation_col` - column name of the perturbations
* `control_ids` - list of control ids 
* `cell_groups` - subset of cell groups to run DEGs on 
* `perturbations` - defaults to perturbation
* `cores`

```
genes_within_cell_state = compare_genes_within_state_graph(chem10_ccs, 
                                                           perturbation_col = "drug_target", 
                                                           control_ids = c("Control"), 
                                                           perturbations = c("TGFB", "Shh", "BMP", "Notch", "FGF", "Wnt", "RA"),
                                                           cores = 6)
```