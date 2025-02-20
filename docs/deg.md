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
                                         pf_cell_state_graph)
```

The output of this table will look like this: 

| cell_state                     | gene_class_scores        |
|--------------------------------|-------------------------|
| pectoral fin condensate        | `<tibble [14,899 × 5]>`   |
| pectoral fin distal mesenchyme | `<tibble [14,899 × 5]>`   |
| pectoral fin central cells     | `<tibble [14,899 × 5]>`   |
| pectoral fin bud mesoderm      | `<tibble [14,899 × 5]>`   |
| pectoral fin cleithrum         | `<tibble [14,899 × 5]>`   |
| pectoral fin bud progenitor    | `<tibble [14,899 × 5]>`   |

To unnest the dataframe: 

```
pf_graph_degs %>% 
  tidyr::unnest(gene_class_scores) %>% 
  filter(pattern_activity_score > 1) %>%
  filter(interpretation == "Selectively activated")
  
```

| cell_state              | gene_id             | data    | interpretation       | pattern_activity_score | gene_short_name |
|-------------------------|---------------------|---------|----------------------|------------------------|-----------------|
| pectoral fin condensate | ENSDARG000000062…  | `<tibble>` | Selectively activated      | 1.02                   | ell2            |
| pectoral fin condensate | ENSDARG000000099…  | `<tibble>` | Selectively activated      | 2.16                   | slc38a5a        |
| pectoral fin condensate | ENSDARG000000106…  | `<tibble>` | Selectively activated      | 1.86                   | clic2           |
| pectoral fin condensate | ENSDARG000000116…  | `<tibble>` | Selectively activated      | 1.19                   | slc26a2         |
| pectoral fin condensate | ENSDARG000000124…  | `<tibble>` | Selectively activated      | 3.77                   | col11a2         |
| pectoral fin condensate | ENSDARG000000309…  | `<tibble>` | Selectively activated      | 2.00                   | mybl1           |

We can check some of these markers by plotting them either in the UMAP space or on our platt graph:


```
plot_cells(pf_ccs@cds, genes = c())
```

![](assets/degs_over_graph.png)

```
plot_gene_expr(pf_cell_state_graph, genes = c())
```

![](assets/degs_over_graph.png)

_For more information about plotting on a Platt graph, see our [plotting page](https://cole-trapnell-lab.github.io/platt/plotting)._

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