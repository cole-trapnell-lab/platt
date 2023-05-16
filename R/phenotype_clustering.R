
library(tidyr)
library(dplyr)
library(monocle3)
#library(data.table)
library(future)
library(purrr)
library(ggplot2)
library(msigdbr)
library(ggtext)


drive_get("https://drive.google.com/file/d/1-05lHha1Z2dE7ufwC2AbZ_1A8GZHSriX") %>%
  drive_download(overwrite = T)

sing_mut_df = data.table::fread("clean_zfin_single-mut_with-ids_phenotype_df.csv", sep = ",", data.table = F, stringsAsFactors = F)

struct_df = sing_mut_df %>%
  group_by(aff_struct_super_1) %>%
  tally(name = "aff_count")

struct_df %>%
  arrange(aff_count) %>% dim()

filt_structs = struct_df %>%
  filter(aff_count > 4) %>%
  pull(aff_struct_super_1)

sel_df = sing_mut_df %>%
  filter(aff_struct_super_1 %in% filt_structs & phen_tag == "abnormal") %>%
  select(gene, aff_struct_super_1, val)

mat = reshape2::acast(sel_df, gene ~ aff_struct_super_1, fill = 0, value.var = "val")

meta_df = as.data.frame(rownames(mat))
colnames(meta_df) = c("gene")
meta_df = meta_df %>%
  left_join(sing_mut_df %>%
              group_by(gene) %>%
              tally(name = "aff_structures"),
            by = "gene")

rownames(meta_df) = meta_df$gene

rowdat = as.data.frame(rownames(t(mat)))
colnames(rowdat) = c("structure")
rownames(rowdat) = rowdat$structure

phen_cds = new_cell_data_set(expression_data = t(mat),
                             cell_metadata = meta_df,
                             gene_metadata = rowdat)

phen_cds = preprocess_cds(phen_cds, method = "LSI") %>%
  reduce_dimension(umap.n_neighbors = 5L, preprocess_method = "LSI") %>%
  cluster_cells(resolution = 1e-3)
colData(phen_cds)$umap1 = reducedDim(x = phen_cds,
                                     type = "UMAP")[,1]
colData(phen_cds)$umap2 = reducedDim(x = phen_cds,
                                     type = "UMAP")[,2]
colData(phen_cds)$group = clusters(phen_cds)

plot_cells(phen_cds, cell_size = 1,
           color_cells_by = "cluster", group_label_size = 4)
