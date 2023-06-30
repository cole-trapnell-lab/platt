subcds_from_coldata = function(cds, coldata, partiton_id) {
  
  sub_coldata = coldata %>% filter(partition == partiton_id)
  sub_cds = cds[, colData(cds)$cell %in% sub_coldata$cell]
  col_names = sub_coldata %>% colnames
  umap_cols = col_names[grepl("partition_umap", col_names)]
  reducedDims(sub_cds)[["UMAP"]] = sub_coldata %>% select(all_of(umap_cols)) %>% as.matrix
  colData(sub_cds)$cluster = sub_coldata$cluster
  return(sub_cds)
}

detect_outlier_cells = function(cds, prefix , k=10) {
  
  # build annoy index
  cds = make_cds_nn_index(cds, reduction_method = "UMAP")
  
  # save it
  # save_transform_models(cds, paste0(prefix, "mt_umap_nn_models"))
  
  # throw them away or just NA everything?
  
  query_ann = cds@reduce_dim_aux$UMAP$nn_index$annoy$nn_index$annoy_index
  query_dims = reducedDims(cds)[["UMAP"]]
  query_res = uwot:::annoy_search(query_dims, k = k + 1, ann = query_ann)
  
  df = lapply(1:nrow(query_dims), function(i) {
    # remember to ignore first index
    neighbor_indices = query_res$idx[i,2:(k+1)]
    neighbor_dists = query_res$dist[i,2:(k+1)]
    data.frame("idx"=i,
               "mean_nn_dist" = mean(neighbor_dists),
               "max_nn_dist" = max(neighbor_dists),
               "median_nn_dist" = median(neighbor_dists),
               "min_nn_dist" = min(neighbor_dists)
    )
  }) %>% bind_rows()
  
  colData(cds)$mean_nn_dist = df$mean_nn_dist
  colData(cds)$median_nn_dist = df$median_nn_dist
  colData(cds)$min_nn_dist = df$min_nn_dist
  colData(cds)$max_nn_dist = df$max_nn_dist
  
  # fig = plot_cells_3d(cds, color_cells_by = "mean")
  # saveWidget(fig, paste0("plots/", prefix,"_mean_nn_dist.html"))
  
  # p = colData(cds) %>% as.data.frame %>% ggplot(aes(mean)) + geom_density() + geom_vline(xintercept = 0.1)
  # ggsave(p, paste0("plots", prefix, "_mean_nn_dist.png"))
  
  return(cds)
}

drop_outlier_cells = function(cds, pct = 99) {
  
  cds = detect_outlier_cells(cds)
  max_dist_quantiles = quantile(colData(cds)$max_nn_dist, probs = seq(0, 1, 0.01), na.rm=TRUE)
  percentile = paste0(pct, "%")
  max_dist_thresh = max_dist_quantiles[percentile]
  colData(cds)$outlier = colData(cds)$max_nn_dist > max_dist_thresh
  # colData(cds)$outlier %>% sum()
  cds = cds[,colData(cds)$outlier == FALSE]
  return(cds)
}


# what is the distance of that cell to the centroid
# of all cells with that label
dist_from_other_cell_type_sub_cells = function(ref_coldata, centroids_per_cts, i, cell_type_label) {
  curr_loc = ref_coldata[i,c("umap3d_1","umap3d_2", "umap3d_3")]
  centroid_loc = centroids_per_cts[cell_type_label,]
  eucl.dist = dist(rbind(curr_loc,centroid_loc) %>% as.matrix)
  return(eucl.dist[[1]])
}

# this currently needs to be run on the full cds in the global space

clean_up_labels = function(cds,
                           old_colname,
                           new_colname = old_colname,
                           return_na = FALSE,
                           k = 10,
                           max_dist_fct = 10) {
  
  # build a knn
  cds = make_cds_nn_index(cds, reduction_method = "UMAP")
  query_ann = cds@reduce_dim_aux$UMAP$nn_index$annoy$nn_index$annoy_index
  query_dims = reducedDims(cds)[["UMAP"]]
  query_res = uwot:::annoy_search(query_dims, k = k + 1, ann = query_ann)
  ref_coldata = colData(cds) %>% as.data.frame()
  ref_coldata = ref_coldata %>% mutate(rn = row_number())
  
  ref_coldata$umap3d_1 = reducedDim(x = cds, type = "UMAP")[,1]
  ref_coldata$umap3d_2 = reducedDim(x = cds, type = "UMAP")[,2]
  ref_coldata$umap3d_3 = reducedDim(x = cds,type = "UMAP")[,3]
  
  
  # where are the rest of those labels
  centroids_per_cts = reducedDims(cds)[["UMAP"]] %>% as.data.frame
  centroids_per_cts$aggregate_col = colData(cds)[[old_colname]]
  centroids_per_cts = aggregate(.~aggregate_col, data = centroids_per_cts, FUN=mean)
  centroids_per_cts = centroids_per_cts %>% tibble::column_to_rownames("aggregate_col")
  names(centroids_per_cts) = c("umap3d_1","umap3d_2", "umap3d_3")
  
  #
  new_labels = sapply(ref_coldata$rn, function(i) {
    # print(i)
    self_label = ref_coldata[query_res$idx[i,1], old_colname]
    neighbor_indices = query_res$idx[i,2:(k+1)]
    neighbor_labels = ref_coldata[neighbor_indices, old_colname]
    neighbor_dists = query_res$dist[i,2:(k+1)]
    # mean_dist = mean(neighbor_dists)
    max_dist = max(neighbor_dists)
    min_dist = min(neighbor_dists)
    majority_label = monocle3:::which_mode(neighbor_labels)
    eucl.dist = dist_from_other_cell_type_sub_cells(ref_coldata, centroids_per_cts, i, self_label)
    
    if (is.na(self_label)){
      self_label = ""
    }
    
    if (self_label != "") {
      eucl.dist = dist_from_other_cell_type_sub_cells(ref_coldata, centroids_per_cts, i, self_label)
    } else {
      # we can keep these as ""
      eucl.dist = 0
    }
    
    
    # only change it if far from like-things and clear majority
    if (eucl.dist > max_dist_fct*max_dist & !is.na(majority_label)) {
      majority_label
    } else if (return_na) {
      NA_character_
    } else {
      self_label
    }
    
  })
  
  colData(cds)[[new_colname]] = new_labels
  return(cds)
}
