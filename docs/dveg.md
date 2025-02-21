# Deviantly expressed genes (DvEGs)

<p align="center">
![](assets/DvEG.png){width=75%}
</p>

Genes are “deviantly expressed” (DvEG) in each perturbation, when they are upregulated during a transition in wild-type, but are under- or overexpressed in perturbed cells undergoing that same transition relative to controls. 

```

dvegs = calculate_dvegs(perturb_degs, 
                        perturbation_table, 
                        ref_abundances, 
                        ref_degs)

```



