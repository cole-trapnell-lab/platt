
# Calculate the probability vector
#' @noRd
makeprobsvec <- function(p) {
  phat <- p/sum(p)
  phat[is.na(phat)] = 0
  phat
}

# Calculate the probability matrix for a relative abundance matrix
#' @noRd
makeprobs <- function(a) {
  colSums<-apply(a,2,sum)
  b <- Matrix::t(Matrix::t(a)/colSums)
  b[is.na(b)] = 0
  b
}

# Calculate the Shannon entropy based on the probability vector
# shannon.entropy <- function(p) {
#   if (min(p) < 0 || (p) <=0)
#     return(Inf)
#   p.norm <- p[p>0]/sum(p)
#   -sum(log2(p.norm)*p.norm)
# }

#' @noRd
shannon_entropy <- function(p) {
  #if (Matrix::rowMin(p) < 0 || (p) <=0)
  #  return(Inf)
  p.norm <- p / Matrix::rowSums(p)
  lg_pnorm = log2(p.norm) * p.norm
  lg_pnorm[p.norm == 0] = 0
  SE = -Matrix::rowSums(lg_pnorm)
  return (SE)
}

# Calculate the Jensen-Shannon distance for two probability distribution
#' @noRd
js_dist_to_pattern <- function (x, pattern)
{
  stopifnot(ncol(x) == length(pattern))
  avg_x_pattern = sweep(x, 2, pattern, "+") / 2
  JSdiv = shannon_entropy(avg_x_pattern) -
    (shannon_entropy(x) + shannon_entropy(matrix(pattern, nrow=1))) * 0.5
  JSdiv[is.infinite(JSdiv)] = 1
  JSdiv[JSdiv < 0] = 0
  JSdist <- sqrt(JSdiv)
  pattern_match_score = 1 - JSdist
  return(pattern_match_score)
}

# Measure the degree of upregulation
#' @noRd
measure_upregulation_effect <- function (self_estimate, other_estimates)
{
  effect_size = self_estimate - matrixStats::rowMaxs(other_estimates)
  return(effect_size)
}

measure_downregulation_effect <- function (self_estimate, other_estimates)
{
  effect_size = matrixStats::rowMaxs(other_estimates) - self_estimate
  return(effect_size)
}

measure_maintenance_effect <- function (self_estimate, other_estimates)
{
  effect_size = 1/matrixStats::rowMaxs(abs(other_estimates - as.vector(self_estimate)))
  return(effect_size)
}


#' get the parent(s) of a state in a state transition graph
#' @noRd
get_parents = function(state_graph, cell_state){
  parents = igraph::neighbors(state_graph, cell_state, mode="in")
  if (length(parents) > 0)
    return (parents$name)
  else
    return (c())
}

#' get the children of a state in a state transition graph
#' @noRd
get_children = function(state_graph, cell_state){
  children = igraph::neighbors(state_graph, cell_state, mode="out")
  if (length(children) > 0)
    return (children$name)
  else
    return (c())
}

#' get the siblings of a state in a state transition graph
#' @noRd
get_siblings = function(state_graph, cell_state){
  parents = get_parents(state_graph, cell_state)
  if (length(parents) > 0){
    siblings = igraph::neighbors(state_graph, parents, mode="out")
    siblings = setdiff(siblings$name, cell_state) #exclude self
    return(siblings)
  } else{
    return (c())
  }
  
}


#' @noRd
score_genes_for_expression_pattern <- function(cell_state, gene_patterns, state_graph, estimate_matrix, state_term="cell_group", cores=1){
  
  parents = get_parents(state_graph, cell_state) #igraph::neighbors(state_graph, cell_state, mode="in")
  parents = intersect(parents, colnames(estimate_matrix))
  
  children = get_children(state_graph, cell_state)#igraph::neighbors(state_graph, cell_state, mode="out")
  children = intersect(children, colnames(estimate_matrix))
  
  siblings = get_siblings(state_graph, cell_state)#igraph::neighbors(state_graph, parents, mode="out")
  siblings = intersect(siblings, colnames(estimate_matrix))
  
  #states_in_contrast = c(cell_state, parents, children, siblings) %>% unique()
  
  expr_df = tibble(gene_id=row.names(estimate_matrix))
  
  self_estimates = estimate_matrix[gene_patterns$gene_id, cell_state, drop=FALSE]
  parents_estimates = estimate_matrix[gene_patterns$gene_id, c(parents), drop=FALSE]
  parents_and_sib_estimates = estimate_matrix[gene_patterns$gene_id, c(parents, siblings), drop=FALSE]
  
  self_and_parent = exp(estimate_matrix[gene_patterns$gene_id, c(cell_state, parents), drop=FALSE])
  self_and_parent = self_and_parent / Matrix::rowSums(self_and_parent) #normalize so that rows sum to 1
  
  self_parent_sibs = exp(estimate_matrix[gene_patterns$gene_id, c(cell_state, parents, siblings), drop=FALSE])
  self_parent_sibs = self_parent_sibs / Matrix::rowSums(self_parent_sibs) #normalize so that rows sum to 1
  
  num_parents = length(parents)
  num_siblings = length(siblings)
  gene_patterns_and_scores = gene_patterns %>%
    tidyr::unnest(interpretation) %>%
    #group_by(interpretation) %>%
    mutate(pattern_match_score =
             case_when(interpretation %in% c("Absent") ~ 0,
                       interpretation %in% c("Maintained") ~ js_dist_to_pattern(self_and_parent,c(1, rep(1, num_parents))),
                       interpretation %in% c("Selectively maintained", "Specifically maintained") ~ js_dist_to_pattern(self_parent_sibs, c(1, rep(1, num_parents), rep(0, num_siblings))),
                       interpretation %in% c("Upregulated", "Activated") ~ js_dist_to_pattern(self_and_parent, c(1, rep(0, num_parents))),
                       interpretation %in% c("Selectively upregulated", "Specifically upregulated", "Selectively activated", "Specifically activated") ~ js_dist_to_pattern(self_parent_sibs, c(1, rep(0, num_parents), rep(0, num_siblings))),
                       interpretation %in% c("Downregulated", "Dectivated") ~ js_dist_to_pattern(self_and_parent, c(0, rep(1, num_parents))),
                       interpretation %in% c("Selectively downregulated", "Specifically downregulated", "Selectively deactivated", "Specifically deactivated") ~ js_dist_to_pattern(self_parent_sibs, c(0, rep(1, num_parents), rep(0, num_siblings))),
                       TRUE ~ 0)
    ) %>%
    #group_by(interpretation) %>%
    mutate(pattern_activity_score =
             case_when(interpretation %in% c("Absent") ~ 0,
                       interpretation %in% c("Maintained") ~ measure_maintenance_effect(self_estimates, parents_estimates),
                       interpretation %in% c("Selectively maintained", "Specifically maintained") ~ measure_maintenance_effect(self_estimates, parents_and_sib_estimates),
                       interpretation %in% c("Upregulated", "Activated") ~ measure_upregulation_effect(self_estimates, parents_estimates),
                       interpretation %in% c("Selectively upregulated", "Specifically upregulated", "Selectively activated", "Specifically activated") ~ measure_upregulation_effect(self_estimates, parents_and_sib_estimates),
                       interpretation %in% c("Downregulated", "Dectivated") ~ measure_downregulation_effect(self_estimates, parents_estimates),
                       interpretation %in% c("Selectively downregulated", "Specifically downregulated", "Selectively deactivated", "Specifically deactivated") ~ measure_downregulation_effect(self_estimates, parents_and_sib_estimates),
                       TRUE ~ 0)
    )
  return(gene_patterns_and_scores)
}


#' @noRDS
compare_genes_within_cell_state <- function(cell_state, state_graph){
  
  
  estimates_for_parents = Matrix::rowSums(pb_coeffs$coefficients[,control_ids, drop=F])
  stderrs_for_parents = sqrt(Matrix::rowSums(pb_coeffs$stdev.unscaled[,control_ids, drop=F]^2))
  
  
}

#' @noRd
classify_genes_in_cell_state <- function(cell_state, state_graph, estimate_matrix, stderr_matrix, state_term="cell_group", log_fc_thresh=1, abs_expr_thresh = 1e-3, sig_thresh=0.05, cores=1){
  #expr_self = expr_mat[,cell_state]
  
  parents = get_parents(state_graph, cell_state) #igraph::neighbors(state_graph, cell_state, mode="in")
  parents = intersect(parents, colnames(estimate_matrix))
  
  children = get_children(state_graph, cell_state)#igraph::neighbors(state_graph, cell_state, mode="out")
  children = intersect(children, colnames(estimate_matrix))
  
  siblings = get_siblings(state_graph, cell_state)#igraph::neighbors(state_graph, parents, mode="out")
  siblings = intersect(siblings, colnames(estimate_matrix))
  
  states_in_contrast = c(cell_state, parents, children, siblings) %>% unique()
  
  expr_df = tibble(gene_id=row.names(estimate_matrix))
  
  message("      examining coefficients ", cell_state)
  
  expr_df$expr_self = pnorm(estimate_matrix[,cell_state] - log(abs_expr_thresh), sd = stderr_matrix[,cell_state], lower.tail=FALSE)
  expr_df$expr_self = p.adjust(expr_df$expr_self, method="BH") < sig_thresh
  
  expr_df$expressed_in_parents = NA
  expr_df$expressed_in_siblings = NA
  expr_df$higher_than_parents = NA
  expr_df$lower_than_parents = NA
  expr_df$higher_than_all_siblings = NA
  expr_df$lower_than_all_siblings = NA
  expr_df$higher_than_siblings = NA
  expr_df$lower_than_siblings = NA
  expr_df$expressed_in_children = NA
  expr_df$higher_than_all_children = NA
  expr_df$lower_than_all_children = NA
  expr_df$higher_than_children = NA
  expr_df$lower_than_children = NA
  
  if (length(parents) > 0){
    expressed_in_parents_mat = pnorm(estimate_matrix[,parents, drop=F] - log(abs_expr_thresh), sd = stderr_matrix[,parents, drop=F], lower.tail=FALSE)
    expressed_in_parents_mat = apply(expressed_in_parents_mat, 2, p.adjust, method="BH")
    
    expressed_in_parents_mat = expressed_in_parents_mat < sig_thresh
    expr_df$expressed_in_parents = Matrix::rowSums(expressed_in_parents_mat) > 0
    
    higher_than_parents_stat = -t(sweep(t(estimate_matrix[,parents, drop=F]), 2, as.numeric(estimate_matrix[,cell_state]) , `-`))
    higher_than_parents_pval = pnorm(higher_than_parents_stat,
                                     sd = sqrt(sweep(t(stderr_matrix[,parents, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    higher_than_parents_pval = apply(higher_than_parents_pval, 2, p.adjust, method="BH")
    
    higher_than_parents_mat = abs(higher_than_parents_stat) > log_fc_thresh & higher_than_parents_pval < sig_thresh
    expr_df$higher_than_parents = Matrix::rowSums(higher_than_parents_mat) > 0
    
    lower_than_parents_pval = pnorm(-higher_than_parents_stat,
                                    sd = sqrt(sweep(t(stderr_matrix[,parents, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    lower_than_parents_pval = apply(lower_than_parents_pval, 2, p.adjust, method="BH")
    
    lower_than_parents_mat = abs(higher_than_parents_stat) > log_fc_thresh & lower_than_parents_pval < sig_thresh
    expr_df$lower_than_parents = Matrix::rowSums(lower_than_parents_mat) > 0
  }else{
    expr_df$expressed_in_parents = NA
    expr_df$expressed_in_siblings = NA
    expr_df$higher_than_parents = NA
    expr_df$lower_than_parents = NA
    expr_df$higher_than_all_siblings = NA
    expr_df$lower_than_all_siblings = NA
    expr_df$higher_than_siblings = NA
    expr_df$lower_than_siblings = NA
    expr_df$expressed_in_children = NA
    expr_df$higher_than_all_children = NA
    expr_df$lower_than_all_children = NA
    expr_df$higher_than_children = NA
    expr_df$lower_than_children = NA
  }
  
  if (length(siblings) > 0){
    expressed_in_siblings_mat = pnorm(estimate_matrix[,siblings, drop=F] - log(abs_expr_thresh), sd = stderr_matrix[,siblings, drop=F], lower.tail=FALSE)
    expressed_in_siblings_mat = apply(expressed_in_siblings_mat, 2, p.adjust, method="BH")
    
    expressed_in_siblings_mat = expressed_in_siblings_mat < sig_thresh
    expr_df$expressed_in_siblings = Matrix::rowSums(expressed_in_siblings_mat) > 0
    
    higher_than_siblings_stat = -t(sweep(t(estimate_matrix[,siblings, drop=F]), 2, as.numeric(estimate_matrix[,cell_state]) , `-`))
    higher_than_siblings_pval = pnorm(higher_than_siblings_stat,
                                      sd = sqrt(sweep(t(stderr_matrix[,siblings, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    higher_than_siblings_pval = apply(higher_than_siblings_pval, 2, p.adjust, method="BH")
    
    higher_than_siblings_mat = abs(higher_than_siblings_stat) > log_fc_thresh & higher_than_siblings_pval < sig_thresh
    expr_df$higher_than_all_siblings = Matrix::rowSums(higher_than_siblings_mat) == ncol(higher_than_siblings_pval)
    expr_df$higher_than_siblings = Matrix::rowSums(higher_than_siblings_mat) > 0
    
    lower_than_siblings_pval = pnorm(-higher_than_siblings_stat,
                                     sd = sqrt(sweep(t(stderr_matrix[,siblings, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    lower_than_siblings_pval = apply(lower_than_siblings_pval, 2, p.adjust, method="BH")
    
    lower_than_siblings_mat = abs(higher_than_siblings_stat) > log_fc_thresh & lower_than_siblings_pval < sig_thresh
    expr_df$lower_than_all_siblings = Matrix::rowSums(lower_than_siblings_mat) == ncol(lower_than_siblings_mat)
    expr_df$lower_than_siblings = Matrix::rowSums(lower_than_siblings_mat) > 0
    
    
  }else{
    expr_df$expressed_in_siblings = NA
    expr_df$higher_than_all_siblings = NA
    expr_df$lower_than_all_siblings = NA
    expr_df$higher_than_siblings = NA
    expr_df$lower_than_siblings = NA
  }
  
  if (length(children) > 0){
    expressed_in_children_mat = pnorm(estimate_matrix[,children, drop=F] - log(abs_expr_thresh), sd = stderr_matrix[,children, drop=F], lower.tail=FALSE)
    expressed_in_children_mat = apply(expressed_in_children_mat, 2, p.adjust, method="BH")
    
    expressed_in_children_mat = expressed_in_children_mat < sig_thresh
    expr_df$expressed_in_children = Matrix::rowSums(expressed_in_children_mat) > 0
    
    higher_than_children_stat = -t(sweep(t(estimate_matrix[,children, drop=F]), 2, as.numeric(estimate_matrix[,cell_state]) , `-`))
    higher_than_children_pval = pnorm(higher_than_children_stat,
                                      sd = sqrt(sweep(t(stderr_matrix[,children, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    higher_than_children_pval = apply(higher_than_children_pval, 2, p.adjust, method="BH")
    
    higher_than_children_mat = abs(higher_than_children_stat) > log_fc_thresh & higher_than_children_pval < sig_thresh
    expr_df$higher_than_all_children = Matrix::rowSums(higher_than_children_mat) == ncol(higher_than_children_pval)
    expr_df$higher_than_children = Matrix::rowSums(higher_than_children_mat) > 0
    
    lower_than_children_pval = pnorm(-higher_than_children_stat,
                                     sd = sqrt(sweep(t(stderr_matrix[,children, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    lower_than_children_pval = apply(lower_than_children_pval, 2, p.adjust, method="BH")
    
    lower_than_children_mat = abs(higher_than_children_stat) > log_fc_thresh & lower_than_children_pval < sig_thresh
    expr_df$lower_than_all_children = Matrix::rowSums(lower_than_children_mat) == ncol(lower_than_children_mat)
    expr_df$lower_than_children = Matrix::rowSums(lower_than_children_mat) > 0
    
    
  }else{
    expr_df$expressed_in_children = NA
    expr_df$higher_than_all_children = NA
    expr_df$lower_than_all_children = NA
    expr_df$higher_than_children = NA
    expr_df$lower_than_children = NA
  }
  
  expr_df = expr_df %>% tidyr::nest(data = !gene_id)
  
  message("      interpreting patterns")
  interpret_expression_pattern = function(pat_df){
    if (pat_df$expr_self){
      if (is.na(pat_df$expressed_in_parents)){
        # no parents, therefore no siblings
        #return ("Maintained")
        if (is.na(pat_df$expressed_in_children)){
          return("Maintained")
        } else {
          # no parent, but there are children
          if (pat_df$expressed_in_children == FALSE | pat_df$higher_than_all_children){
            # Higher than parent, and higher than children
            return("Precursor-specific")
          }
          else if (pat_df$higher_than_children){
            # no parent, higher than children
            return("Precursor-specific")
          }
          else if(pat_df$lower_than_children){
            # no parent, but lower than children
            return("Precursor-depleted")
          }
          else { # no parent same as children
            return("Maintained")
          }
        }
      }else if (pat_df$expressed_in_parents){
        # Expressed in self and parent
        if (is.na(pat_df$expressed_in_siblings)){
          # Expressed in self and parent and there are no siblings
          if (pat_df$higher_than_parents){
            if (is.na(pat_df$expressed_in_children)){
              return("Upregulated")
            } else {
              # there are children
              if (pat_df$expressed_in_children == FALSE | pat_df$higher_than_all_children){
                # Higher than parent, and higher than siblings
                return("Transiently upregulated")
              }
              else if(pat_df$lower_than_all_children){
                # lower than children
                return("Increasingly upregulated")
              }
              else { # same as children
                return("Upregulated")
              }
            }
          }
          else if(pat_df$lower_than_parents){
            if (is.na(pat_df$expressed_in_children)){
              return("Downregulated")
            } else {
              # there are children
              if (pat_df$lower_than_all_children){
                # Lower than parent, and lower than children
                return("Decreasingly downregulated")
              }
              else if(pat_df$lower_than_all_children){
                # lower than children
                return("Transiently downregulated")
              }
              else { # same as children
                return("Downregulated")
              }
            }
          }else{
            if (is.na(pat_df$expressed_in_children)){
              return("Maintained")
            } else {
              # same as parent, and there are children
              if (pat_df$expressed_in_children == FALSE | pat_df$higher_than_all_children){
                # Higher than parent, and higher than children
                return("Precursor-specific")
              }
              else if(pat_df$lower_than_all_children){
                # no parent, but lower than children
                return("Precursor-depleted")
              }
              else { # no parent same as children
                return("Maintained")
              }
            }
          }
        } else {
          # Expressed in self and parent and there are siblings
          if (pat_df$higher_than_parents){
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Higher than parent, and higher than siblings
              return("Specifically upregulated")
            }
            else if (pat_df$higher_than_siblings){
              # Higher than parent, and higher than siblings
              return("Selectively upregulated")
            }
            else if(pat_df$lower_than_siblings){
              # Higher than parent, but lower than siblings
              return("Upregulated")
            }
            else { # higher than parent, same as siblings
              return("Upregulated")
            }
          }
          else if(pat_df$lower_than_parents){
            if (pat_df$expressed_in_siblings == TRUE & pat_df$lower_than_all_siblings){
              # Lower than parent, and higher than siblings
              return("Specifically downregulated")
            }
            else if (pat_df$expressed_in_siblings == TRUE & pat_df$lower_than_siblings){
              # Lower than parent, and higher than some siblings
              return("Selectively downregulated")
            }
            else { # lower than parent, same as or higher than siblings
              return("Downregulated")
            }
          }
          else { # same as parent
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Same as parent, and higher than all siblings
              return("Specifically maintained")
            }
            else if (pat_df$higher_than_siblings){
              # Same as parent, and higher than some siblings
              return("Selectively maintained")
            }
            else if(pat_df$lower_than_all_siblings){
              # Same as parent, but lower than siblings
              return("Maintained")
            }
            else { # same as parent, same as siblings
              return("Maintained")
            }
          }
        }
        
      }else{
        # expressed in self but not in parent
        if (is.na(pat_df$expressed_in_siblings)){
          # Expressed in self, not in parent and there are no siblings
          if (pat_df$higher_than_parents)
            return("Activated")
          else if(pat_df$lower_than_parents)
            return("Downregulated") # shouldn't happen
          else
            return("Activated") # might happen if its above threshold but not significantly above parent (and parent is below thresh)
        } else {
          # expressed in self, not in parent, and there are siblings
          if (pat_df$higher_than_parents){ # Expressed in self and higher than parent and there are siblings
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Higher than parent, and higher than all siblings
              return("Specifically activated")
            } else if (pat_df$higher_than_siblings){
              # Higher than parent, and higher than some siblings
              return("Selectively activated")
            }
            else if(pat_df$lower_than_all_siblings){
              # Higher than parent (which is off), but lower than all siblings
              return("Activated")
            }
            if(pat_df$lower_than_siblings){
              # Higher than parent (which is off), but lower than some siblings
              return("Activated")
            }
            else { # Higher than parent (which is off), same as siblings
              return("Activated")
            }
          }
          else if(pat_df$lower_than_parents){
            # gene is expressed, lower in the parent (which is off)
            if (pat_df$higher_than_all_siblings){
              # Lower than parent, and higher than all siblings
              return("Absent") # shouldn't happen
            }
            else if (pat_df$higher_than_siblings){
              # Lower than parent, and higher than some siblings
              return("Absent") # shouldn't happen
            }
            else if(pat_df$lower_than_all_siblings){
              # Lower than parent and  lower than all siblings
              return("Absent")
            }
            else if(pat_df$lower_than_siblings){
              # Lower than parent and  lower than some siblings
              return("Absent")
            }
            else { # Lower than parent, same as siblings
              return("Absent")
            }
          }
          else { # same as parent (which is off)
            if (pat_df$higher_than_all_siblings){
              # Same as parent, and higher than all siblings
              return("Absent")
            }
            else if (pat_df$higher_than_siblings){
              # Same as parent, and higher than some siblings
              return("Absent")
            }
            else if(pat_df$lower_than_all_siblings){
              # Same as parent, but lower than all siblings
              return("Absent")
            }
            else if(pat_df$lower_than_siblings){
              # Same as parent, but lower than some siblings
              return("Absent")
            }
            else { # same as parent, same as siblings
              return("Absent")
            }
          }
        }
      }
      return ("Expressed")
    }else{
      # Not expressed in self
      if (is.na(pat_df$expressed_in_parents)){
        # no parents, therefore no siblings
        return ("Absent")
      }else if (pat_df$expressed_in_parents){
        # Not expressed in self, but expressed in parents
        if (is.na(pat_df$expressed_in_siblings)){
          # Not expressed in self, expressed parent and there are no siblings
          if(pat_df$lower_than_parents)
            return("Deactivated")
          else
            return("Absent") # shouldn't happen
        } else {
          # Not expressed in self, expressed in parent and there are siblings
          if(pat_df$lower_than_parents){
            # Lower than parent
            if(pat_df$lower_than_all_siblings){
              # Lower than parent and  lower than siblings
              return("Specifically deactivated")
            }
            else if(pat_df$lower_than_siblings){
              # Lower than parent and  lower than siblings
              return("Selectively deactivated")
            }
            return("Deactivated")
          }
          else {
            #Not expressed in self, not lower than parent
            return ("Absent")
          }
        }
      }else{
        # Not expressed in self or parents
        return ("Absent")
      }
      return ("Absent")
    }
    return ("Absent")
    #match_row = match(data.frame(t(pat_df)), data.frame(t(interp_table)))
    #interpetation[match_row]
  }
  #debug(interpret_expression_pattern)
  expr_df = expr_df %>% mutate(interpretation = purrr::map(.f = purrr::possibly(
    interpret_expression_pattern, NA_character_), .x = data))
  message("      completed ", cell_state)
  return(expr_df)
}
#debug(classify_genes_in_cell_state)




#' helper function 
#' @param ccs
#' @param cell_states
#' @param label_nodes_by
#' @export
unnest_degs = function(ccs, 
                       cell_states, 
                       label_nodes_by) {
  
  cell_states = cell_states %>%
    filter(is.na(gene_classes) == FALSE) %>%
    tidyr::unnest(gene_class_scores) %>% 
    dplyr::select(cell_state, gene_id, interpretation, pattern_match_score, pattern_activity_score)
  
  cell_states = left_join(cell_states,
                          rowData(ccs@cds) %>% as_tibble %>% select(id, gene_short_name), by=c("gene_id"="id"))
  
  cell_states = left_join(cell_states,
                          collect_psg_node_metadata(ccs, color_nodes_by=NULL, group_nodes_by = NULL, label_nodes_by=label_nodes_by), 
                          by=c("cell_state"="id")) %>%
    dplyr::rename(cell_type=label_nodes_by)
  
  return(cell_states)
}


#' #' similar to fit_genotype_ccm
#' #' @param ccs
#' #' @param genotype
#' #' @param ctrl_ids
#' #' @param perturbation_col
#' #' @param interval_col
#' #' @param cell_groups
#' fit_genotype_deg = function(ccs,
#'                             genotype,
#'                             ctrl_ids = c("Control"),
#'                             perturbation_col = "perturbation",
#'                             interval_col = "timepoint",
#'                             assembly_time_start = 18,
#'                             assembly_time_stop= 72,
#'                             cell_groups = NULL,
#'                             cores = 1,
#'                             min_samples_detected = 2,
#'                             min_cells_per_pseudobulk = 3,
#'                             ... ) {
#'   
#'   expts = unique(colData(ccs)$expt)
#'   
#'   if (is.null(assembly_time_start)){
#'     knockout_time_start = min(colData(ccs)[[interval_col]])
#'   }else{
#'     knockout_time_start = assembly_time_start
#'   }
#'   
#'   if (is.null(assembly_time_stop)){
#'     knockout_time_stop = max(colData(subset_ccs)[[interval_col]])
#'   }else{
#'     knockout_time_stop = assembly_time_stop
#'   }
#'   
#'   num_knockout_timepoints = length(unique(colData(subset_ccs)[[interval_col]]))
#'   
#'   message(paste("\ttime range:", knockout_time_start, "to", knockout_time_stop))
#'   subset_ccs = ccs[,( replace_na(colData(ccs)[[perturbation_col]] == genotype, F) | colData(ccs)[[perturbation_col]] %in% ctrl_ids) & colData(ccs)$expt %in% expts]
#'   
#'   colData(subset_ccs)$knockout = colData(subset_ccs)[[perturbation_col]] == genotype
#'   subset_ccs = subset_ccs[,(colData(subset_ccs)[[interval_col]] >= knockout_time_start & colData(subset_ccs)[[interval_col]] <= knockout_time_stop)]
#'   
#'   colData(subset_ccs@cds)$knockout = colData(subset_ccs@cds)[[perturbation_col]] == genotype
#'   
#'   pb_cds = hooke:::pseudobulk_ccs_for_states(subset_ccs)
#'   
#'   # subset to genes that are expressed over a certain min value
#'   expr_over_thresh = normalized_counts(pb_cds, "size_only", pseudocount = 0)
#'   genes_to_test = which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
#'   pb_cds = pb_cds[genes_to_test,]
#'   
#'   pseudobulks_to_test = which(colData(pb_cds)$num_cells_in_group >= min_cells_per_pseudobulk)
#'   pb_cds = pb_cds[,pseudobulks_to_test]
#'   
#'   
#'   pb_cds = hooke:::add_covariate(subset_ccs, pb_cds, "knockout")
#'   
#'   
#'   if (is.null(cell_groups)) {
#'     cell_groups = rownames(counts(ccs))
#'   }
#'   
#'   if (is.null(perturbations)) {
#'     perturbations = unique(colData(pb_cds)[[perturbation_col]])
#'   }
#'   
#'   # cell_group_models = lapply(perturbations, function(perturbation) {
#'   cell_group_models = lapply(cell_groups, function(cell_group) {
#'     
#'     # perturb_pb_cds = pb_cds[, colData(pb_cds)[[perturbation_col]] == perturbation]
#'     cg_pb_cds = pb_cds[, colData(pb_cds)[[state_term]] == cell_group]
#'   
#'     # message(paste0("fitting regression models for ", perturbation))
#'     message(paste0("fitting regression models for ", cell_group))
#'   
#'     pb_group_models = fit_models(cg_pb_cds,
#'                                  model_formula_str=paste("~ 0 + knockout", ),
#'                                  weights=colData(pb_cds)$num_cells_in_group,
#'                                  cores=cores) %>% dplyr::select(gene_short_name, id, model, model_summary)
#'   
#'     message(paste0("fitting regression models for ", cell_group))
#'     
#'     pb_group_models = coefficient_table(pb_group_models) %>%
#'       dplyr::select(gene_short_name, id, term, estimate, std_err) %>%
#'       mutate(term = stringr::str_replace_all(term, state_term, ""))
#'     estimate_matrix = pb_group_models %>% dplyr::select(id, term, estimate)
#'     estimate_matrix = estimate_matrix %>% mutate(term = factor(term, levels=unique(colData(cg_pb_cds)[,"knockout"])))
#'     estimate_matrix = estimate_matrix %>% tidyr::pivot_wider(names_from=term, values_from=estimate, values_fill=0)
#'   
#'     gene_ids = estimate_matrix$id
#'     estimate_matrix$id = NULL
#'     estimate_matrix = as.matrix(estimate_matrix)
#'     row.names(estimate_matrix) = gene_ids
#'     colnames(estimate_matrix) = as.character(colnames(estimate_matrix))
#'   
#'     stderr_matrix = pb_group_models %>% dplyr::select(id, term, std_err)
#'     estimate_matrix = estimate_matrix %>% mutate(term = factor(term, levels=unique(colData(cg_pb_cds)[,"knockout"])))
#'     stderr_matrix = stderr_matrix %>% tidyr::pivot_wider(names_from=term, values_from=std_err, values_fill=0)
#'   
#'     gene_ids = stderr_matrix$id
#'     stderr_matrix$id = NULL
#'     stderr_matrix = as.matrix(stderr_matrix)
#'     row.names(stderr_matrix) = gene_ids
#'     colnames(stderr_matrix) = as.character(colnames(stderr_matrix))
#'   
#'     # states_to_assess = intersect(as.character(unique(colData(pb_cds)[,state_term])), unlist(igraph::V(state_graph)$name))
#'     # cell_states = tibble(cell_state = states_to_assess)
#'   
#'     # cell_states = cell_states %>%
#'     #   dplyr::mutate(gene_classes = purrr::map(.f = purrr::possibly(
#'     #     classify_genes_in_cell_state, NA_real_), .x = cell_state,
#'     #     state_graph, estimate_matrix, stderr_matrix, state_term,
#'     #     log_fc_thresh=log_fc_thresh,
#'     #     abs_expr_thresh = abs_expr_thresh,
#'     #     sig_thresh=sig_thresh,
#'     #     cores=cores))
#'     # 
#'     # cell_states = cell_states %>%
#'     #   filter(is.na(gene_classes) == FALSE) %>%
#'     #   dplyr::mutate(gene_class_scores = purrr::map2(.f = purrr::possibly(
#'     #     score_genes_for_expression_pattern, NA_real_),
#'     #     .x = cell_state,
#'     #     .y = gene_classes,
#'     #     state_graph,
#'     #     estimate_matrix))
#'       
#'     # cell_states
#' 
#'   })
#'   
#'   return(perturb_group_models)
#'   
#' }
#' 
#' #' to do: do i want to subset by time or do i not care
#' #' @param ccs
#' #' @param perturbation_col
#' #' @param control_ids
#' #' @param cell_groups
#' fit_mt_deg_models = function(ccs,
#'                              perturbation_col,
#'                              control_ids,
#'                              cell_groups = NULL,
#'                              num_threads = 1) {
#'   
#'   ccs@cds_coldata[["perturb_name"]] = ccs@cds_coldata[[perturbation_col]]
#'   
#'   perturb_df = ccs@cds_coldata %>%
#'     as_tibble() %>%
#'     dplyr::select(perturb_name) %>%
#'     filter(!perturb_name %in% ctrl_ids) %>%
#'     distinct()
#'   
#'   perturb_models_tbl = perturb_df %>%
#'     purrr::map(.f = fit_genotype_deg,
#'                .x = perturb_name,
#'                control_ids = control_ids,
#'                cell_groups = cell_groups)
#'   
#'   return(perturb_models_tbl)
#' }



#' Classify each gene's pattern of expression in each state in a state transition graph
#' @export
classify_genes_over_graph <- function(ccs,
                                      state_graph,
                                      gene_ids = NULL,
                                      group_nodes_by=NULL,
                                      assembly_group = NULL, 
                                      label_nodes_by="cell_state", 
                                      log_fc_thresh=1,
                                      abs_expr_thresh = 1e-3,
                                      sig_thresh=0.05,
                                      min_samples_detected = 2,
                                      min_cells_per_pseudobulk = 3,
                                      cores=1,
                                      ...){
  if (is.null(group_nodes_by)){
    pb_cds = hooke:::pseudobulk_ccs_for_states(ccs, cell_agg_fun="sum")
    state_term = "cell_group"
  }else{
    pb_cds = hooke:::pseudobulk_ccs_for_states(ccs, state_col = group_nodes_by, cell_agg_fun="sum")
    state_term = group_nodes_by
  }
  
  # if we want to run it by assembly group
  if (is.null(assembly_group) == FALSE) {
    pb_cds = hooke:::add_covariate(ccs, pb_cds, "assembly_group")
    pb_cds = pb_cds[, colData(pb_cds)$assembly_group == assembly_group]
    
    if (is(state_graph, "igraph")) {
      state_graph = igraph::as_data_frame(state_graph)
    }
    state_graph = state_graph[state_graph$assembly_group == assembly_group,]
  }
  
  
  if (!is(state_graph, "igraph")){
    state_graph = state_graph %>% igraph::graph_from_data_frame()
    
  }
  
  
  #cds_to_test = pb_cds[,as.character(colData(pb_cds)[,state_term]) %in% states_in_model]
  
  #colData(cds_to_test)[,state_term] = factor(as.character(colData(cds_to_test)[,state_term]), levels=states_in_model) # set the "self" state as the reference level
  
  #norm_expr_mat = normalized_counts(pb_cds, "size_only", pseudocount = 0)
  
  if (is.null(gene_ids) == FALSE){
    pb_cds = pb_cds[gene_ids,]
  }
  
  # expr_over_thresh = threshold_expression_matrix(normalized_counts(pb_cds, "size_only", pseudocount = 0), ...)
  expr_over_thresh = normalized_counts(pb_cds, "size_only", pseudocount = 0)
  genes_to_test = which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
  pb_cds = pb_cds[genes_to_test,]
  
  pseudobulks_to_test = which(colData(pb_cds)$num_cells_in_group >= min_cells_per_pseudobulk)
  
  message("fitting regression models")
  pb_cds = pb_cds[,pseudobulks_to_test]
  
  colData(pb_cds)$Size_Factor = colData(pb_cds)$num_cells_in_group
  pb_group_models = fit_models(pb_cds,
                               model_formula_str=paste("~ 0 + ", state_term),
                               # weights=colData(pb_cds)$num_cells_in_group,
                               cores=cores) %>% dplyr::select(gene_short_name, id, model, model_summary, status)
  
  message("      collecting coefficients")
  
  # pb_group_models = coefficient_table(pb_group_models) %>%
  #   dplyr::select(gene_short_name, id, term, estimate, std_err) %>%
  #   mutate(term = stringr::str_replace_all(term, state_term, ""))
  pb_group_models = collect_coefficients_for_shrinkage(pb_cds, pb_group_models, abs_expr_thresh)
  
  
  estimate_matrix = pb_group_models %>% dplyr::select(id, term, estimate)
  estimate_matrix = estimate_matrix %>% mutate(term = factor(term, levels=unique(colData(pb_cds)[,state_term])))
  estimate_matrix = estimate_matrix %>% tidyr::pivot_wider(names_from=term, values_from=estimate, values_fill=0)
  
  gene_ids = estimate_matrix$id
  estimate_matrix$id = NULL
  estimate_matrix = as.matrix(estimate_matrix)
  row.names(estimate_matrix) = gene_ids
  colnames(estimate_matrix) = as.character(colnames(estimate_matrix))
  
  stderr_matrix = pb_group_models %>% dplyr::select(id, term, std_err)
  stderr_matrix = stderr_matrix %>% mutate(term = factor(term, levels=unique(colData(pb_cds)[,state_term])))
  stderr_matrix = stderr_matrix %>% tidyr::pivot_wider(names_from=term, values_from=std_err, values_fill=0)
  
  gene_ids = stderr_matrix$id
  stderr_matrix$id = NULL
  stderr_matrix = as.matrix(stderr_matrix)
  row.names(stderr_matrix) = gene_ids
  colnames(stderr_matrix) = as.character(colnames(stderr_matrix))
  
  #p_val_matrix = pnorm(estimate_matrix - log(abs_expr_thresh), sd = stderr_matrix, lower.tail=FALSE)
  
  #expr_thresh_mat = p_val_matrix < sig_thresh
  
  #cell_states = tibble(cell_state = unlist(igraph::V(state_graph)$name))
  states_to_assess = intersect(as.character(unique(colData(pb_cds)[,state_term])), unlist(igraph::V(state_graph)$name))
  cell_states = tibble(cell_state = states_to_assess)
  
  cell_states = cell_states %>%
    dplyr::mutate(gene_classes = purrr::map(.f = purrr::possibly(
      classify_genes_in_cell_state, NA_real_), .x = cell_state,
      state_graph, estimate_matrix, stderr_matrix, state_term,
      log_fc_thresh=log_fc_thresh,
      abs_expr_thresh = abs_expr_thresh,
      sig_thresh=sig_thresh,
      cores=cores))
  
  cell_states = cell_states %>%
    filter(is.na(gene_classes) == FALSE) %>%
    dplyr::mutate(gene_class_scores = purrr::map2(.f = purrr::possibly(
      score_genes_for_expression_pattern, NA_real_),
      .x = cell_state,
      .y = gene_classes,
      state_graph,
      estimate_matrix))
  
  return(cell_states)
}


#' Classify each gene's pattern of expression in each state in a state transition graph
#' @param ccs a cell_count_set object
#' @param state_graph a graph to run the DEG testing over
#' @param cores 
#' @export
compare_genes_over_graph <- function(ccs,
                                    state_graph,
                                    gene_ids = NULL,
                                    group_nodes_by=NULL,
                                    assembly_group = NULL, 
                                    label_nodes_by="cell_state", 
                                    nuisance_model_formula_str = "0",
                                    log_fc_thresh=1,
                                    abs_expr_thresh = 1e-3,
                                    sig_thresh=0.05,
                                    min_samples_detected = 2,
                                    min_cells_per_pseudobulk = 3,
                                    cores=1,
                                    ...){
  if (is.null(group_nodes_by)){
    pb_cds = hooke:::pseudobulk_ccs_for_states(ccs, cell_agg_fun="sum")
    state_term = "cell_group"
  }else{
    pb_cds = hooke:::pseudobulk_ccs_for_states(ccs, state_col = group_nodes_by, cell_agg_fun="sum")
    state_term = group_nodes_by
  }
  
  # if we want to run it by assembly group
  if (is.null(assembly_group) == FALSE) {
    pb_cds = hooke:::add_covariate(ccs, pb_cds, "assembly_group")
    pb_cds = pb_cds[, colData(pb_cds)$assembly_group == assembly_group]
    
    if (is(state_graph, "igraph")) {
      state_graph = igraph::as_data_frame(state_graph)
    }
    state_graph = state_graph[state_graph$assembly_group == assembly_group,]
  }
  
  
  if (!is(state_graph, "igraph")){
    state_graph = state_graph %>% igraph::graph_from_data_frame()
    
  }
  
  if (is.null(gene_ids) == FALSE){
    pb_cds = pb_cds[gene_ids,]
  }
  
  # expr_over_thresh = threshold_expression_matrix(normalized_counts(pb_cds, "size_only", pseudocount = 0), ...)
  expr_over_thresh = normalized_counts(pb_cds, "size_only", pseudocount = 0)
  genes_to_test = which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
  pb_cds = pb_cds[genes_to_test,]
  
  pseudobulks_to_test = which(colData(pb_cds)$num_cells_in_group >= min_cells_per_pseudobulk)
  
  message("fitting regression models")
  pb_cds = pb_cds[,pseudobulks_to_test]
  
  colData(pb_cds)$Size_Factor = colData(pb_cds)$num_cells_in_group
  
  full_model_str = paste0("~ 0 + cell_group")
  nuisance_model_formula_str = stringr::str_replace_all(nuisance_model_formula_str, "~", "")
  
  if (nuisance_model_formula_str != "0") {
    full_model_str = paste0(full_model_str, " + ", nuisance_model_formula_str)
  }
  
  pb_group_models = fit_models(pb_cds,
                               model_formula_str = full_model_str,
                               cores=cores) %>% dplyr::select(gene_short_name, id, model, model_summary, status)
  
  message("      collecting coefficients")
  
  
  pb_coeffs = collect_coefficients_for_shrinkage(pb_cds, pb_group_models, abs_expr_thresh, term_to_keep = "cell_group") 
  
  
  states_to_assess = intersect(as.character(unique(colData(pb_cds)[,state_term])), unlist(igraph::V(state_graph)$name))
  cell_states = tibble(cell_state = states_to_assess)
  
  
  # cell_states = tibble(term = states_in_contrast)
  # cell_states = cell_states %>% 
  #   mutate(perturb_effects = purrr:::map(.f = contrast_helper, 
  #                                        .x = term,
  #                                        control_effects = estimates_for_parents,
  #                                        control_effects_se = stderrs_for_parents,
  #                                        PEM = pb_coeffs$coefficients, 
  #                                        PSEM = pb_coeffs$stdev.unscaled))
  # 
  # new_estimate_matrix = cell_states %>% tidyr::unnest(perturb_effects) %>%
  #   select(id, term, shrunken_lfc) %>% 
  #   pivot_wider(names_from = term, values_from = shrunken_lfc) %>% 
  #   column_to_rownames("id")
  # 
  # new_stderr_matrix = cell_states %>% tidyr::unnest(perturb_effects) %>% 
  #   select(id, term, shrunken_lfc_se) %>% 
  #   pivot_wider(names_from = term, values_from = shrunken_lfc_se) %>% 
  #   column_to_rownames("id")
  # 
  # p_value_matrix = cell_states %>% tidyr::unnest(perturb_effects) %>% 
  #   select(id, term, p_value) %>% 
  #   pivot_wider(names_from = term, values_from = p_value) %>% 
  #   column_to_rownames("id")
  
  
  cell_states = cell_states %>%
    dplyr::mutate(gene_classes = purrr::map(.f = purrr::possibly(
      compare_genes_in_cell_state, NA_real_), 
      .x = cell_state,
      state_graph = state_graph, 
      estimate_matrix = pb_coeffs$coefficients, 
      stderr_matrix = pb_coeffs$stdev.unscaled,
      state_term = state_term,
      log_fc_thresh=log_fc_thresh,
      abs_expr_thresh = abs_expr_thresh,
      sig_thresh=sig_thresh,
      cores=cores))
  
  
  cell_states = cell_states %>%
    filter(is.na(gene_classes) == FALSE) %>%
    dplyr::mutate(gene_class_scores = purrr::map2(.f = purrr::possibly(
      score_genes_for_expression_pattern, NA_real_),
      .x = cell_state,
      .y = gene_classes,
      state_graph,
      pb_coeffs$coefficients)) %>% 
    select(-gene_classes) # otherwise it is duplicated 
  
  return(cell_states)
  
}



#' #' @export
#' classify_genes_within_state_graph = function(ccs,
#'                                              perturbation_col = "perturbation", 
#'                                              control_ids = c("Control"), 
#'                                              cell_groups = NULL, 
#'                                              assembly_group = NULL, 
#'                                              state_graph = NULL,
#'                                              perturbations = NULL, 
#'                                              gene_ids = NULL,
#'                                              group_nodes_by = NULL,
#'                                              log_fc_thresh = 1,
#'                                              abs_expr_thresh = 1e-3,
#'                                              sig_thresh = 0.05,
#'                                              min_samples_detected = 2,
#'                                              min_cells_per_pseudobulk = 3,
#'                                              cores = 1,
#'                                              ...) {
#'   
#'   # to do make sure that ccs and state graph match 
#'   # check if all nodes in the state graph exist in the cds
#'   
#'   # if (is.null(state_graph) == FALSE) {
#'   #   
#'   # }
#'   
#'   expts = unique(colData(ccs)$expt)
#'   
#'   pb_cds = hooke:::pseudobulk_ccs_for_states(ccs, cell_agg_fun="sum")
#'   colData(pb_cds)$Size_Factor = colData(pb_cds)$num_cells_in_group
#'   
#'   pb_cds = hooke:::add_covariate(ccs, pb_cds, perturbation_col)
#'   colData(pb_cds)[["perturbation"]] = colData(pb_cds)[[perturbation_col]] 
#'   
#'   # to do: if we want to run it by assembly group, must provide a state graph 
#'   if (is.null(assembly_group) == FALSE & is.null(state_graph) == FALSE) {
#'     pb_cds = hooke:::add_covariate(ccs, pb_cds, "assembly_group")
#'     pb_cds = pb_cds[, colData(pb_cds)$assembly_group == assembly_group]
#'     
#'     if (is(state_graph, "igraph")) {
#'       state_graph = igraph::as_data_frame(state_graph)
#'     }
#'     state_graph = state_graph %>% filter(assembly_group == assembly_group)
#'   }
#'   
#'   # ability to subset by perturbation 
#'   if (!is.null(perturbations)) {
#'     pb_cds = pb_cds[, colData(pb_cds)[[perturbation_col]] %in% c(control_ids, perturbations)]
#'   } 
#'   
#'   # subset to genes that are expressed over a certain min value
#'   expr_over_thresh = normalized_counts(pb_cds, "size_only", pseudocount = 0)
#'   genes_to_test = which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
#'   pb_cds = pb_cds[genes_to_test,]
#'   
#'   pseudobulks_to_test = which(colData(pb_cds)$num_cells_in_group >= min_cells_per_pseudobulk)
#'   pb_cds = pb_cds[,pseudobulks_to_test]
#'   
#'   
#'   if (is.null(cell_groups)) {
#'     cell_groups = rownames(counts(ccs))
#'   }
#'   
#'   df = data.frame(cell_group = cell_groups) %>% 
#'     mutate(genes_within_cell_group = purrr::map(.f = purrr:::possibly(classify_genes_within_node,NA_real_),
#'                                                 .x = cell_group, 
#'                                                 pb_cds = pb_cds, 
#'                                                 ccs = ccs, 
#'                                                 cores = cores))
#'   
#'   return(df)
#'   
#' }

collect_coefficients_for_shrinkage <- function(cds, model_tbl, abs_expr_thresh, term_to_keep){
  
  model_tbl = model_tbl %>% dplyr::mutate(dispersion = purrr::map2(.f = purrr::possibly(
    extract_dispersion_helper, NA_real_), .x = model,
    .y = model_summary)) %>% tidyr::unnest(dispersion)
  
  mean_expr_helper <- function(model, new_data, type="response"){
    res = tryCatch({
      mean(stats::predict(model, newdata=new_data, type=type))
    }, error = function(e){
      NA
    })
    return(res)
  }
  
  model_tbl <- model_tbl %>%
    dplyr::mutate(mean_expr = purrr::map(model, mean_expr_helper, new_data=colData(cds
    )%>% as.data.frame)) %>% tidyr::unnest(mean_expr)
  
  disp_fit = mgcv::gam(log(dispersion)~s(log(mean_expr), bs="cs"), data=model_tbl)
  model_tbl$disp_fit = exp(predict(disp_fit))
  #model_tbl$disp_fit = 1
  
  model_tbl = update_summary(model_tbl, dispersion_type="fitted")
  
  
  raw_coefficient_table = coefficient_table(model_tbl)
  
  extract_extra_model_stats = function(model, newdata){
    if (class(model)[1] == "speedglm") {
      extra_stats = tibble(RSS=model$RSS, df.residual=model$df, mean_expr=mean(predict(model, newdata=newdata)))
      return(extra_stats)
    }else{
      extra_stats = tibble(RSS=NA_real_, df.residual=NA_real_, mean_expr=NA_real_)
      return(extra_stats)
    }
  }
  
  extra_model_stats = model_tbl %>%
    dplyr::mutate(extra_stats = purrr::map(.f = purrr::possibly(
      extract_extra_model_stats, NA_real_), .x = model, newdata=colData(cds) %>% as.data.frame)) %>%
    select(id, extra_stats) %>%
    tidyr::unnest(extra_stats)
  
  raw_coefficient_table = left_join(raw_coefficient_table, extra_model_stats) %>%
    left_join(model_tbl %>% select(id, disp_fit, dispersion)) %>%
    #coefficient_table(pb_group_models) %>%
    #dplyr::select(gene_short_name, id, term, estimate, std_err, p_value, status) %>%
    filter(grepl(term_to_keep, term)) %>% 
    mutate(term = stringr::str_replace_all(term, term_to_keep, ""))
  
  estimate_matrix = raw_coefficient_table %>% dplyr::select(id, term, estimate)
  estimate_matrix = estimate_matrix %>% mutate(term = factor(term, levels=unique(colData(cds)[,term_to_keep])))
  estimate_matrix = estimate_matrix %>% tidyr::pivot_wider(names_from=term, values_from=estimate, values_fill=0)
  
  gene_ids = estimate_matrix$id
  estimate_matrix$id = NULL
  estimate_matrix = as.matrix(estimate_matrix)
  row.names(estimate_matrix) = gene_ids
  colnames(estimate_matrix) = as.character(colnames(estimate_matrix))
  
  stderr_matrix = raw_coefficient_table %>% dplyr::select(id, term, std_err)
  stderr_matrix = stderr_matrix %>% mutate(term = factor(term, levels=unique(colData(cds)[,term_to_keep])))
  stderr_matrix = stderr_matrix %>% tidyr::pivot_wider(names_from=term, values_from=std_err, values_fill=0)
  
  gene_ids = stderr_matrix$id
  stderr_matrix$id = NULL
  stderr_matrix = as.matrix(stderr_matrix)
  row.names(stderr_matrix) = gene_ids
  colnames(stderr_matrix) = as.character(colnames(stderr_matrix))
  
  # collect the ids of any genes that threw an exception in fit_models and
  # set their estimates and std_errors to NA
  fail_gene_ids = model_tbl %>% filter(status == "FAIL") %>% pull(id)
  if (length(fail_gene_ids) > 0){
    estimate_matrix[fail_gene_ids,] = log(abs_expr_thresh)
    stderr_matrix[fail_gene_ids,] = Inf
  }
  
  sigma_df = raw_coefficient_table %>% dplyr::select(id, RSS, df.residual, mean_expr, disp_fit, dispersion) %>% distinct() %>% as.data.frame
  row.names(sigma_df) = sigma_df$id
  sigma_df$id = NULL
  sigma_df = sigma_df[row.names(estimate_matrix),]
  sigma_df$sigma = sigma_df$RSS / sigma_df$df.residual
  
  std_dev.unscaled = stderr_matrix^2 /  sigma_df$sigma
  
  coefs_for_shrinkage = list(coefficients = estimate_matrix,
                             stdev.unscaled = stderr_matrix,
                             sigma = sigma_df$sigma,
                             df.residual = sigma_df$df.residual,
                             trend=sigma_df$mean_expr,
                             est_dispersion=sigma_df$dispersion,
                             disp_fit = sigma_df$disp_fit)
  
  return(coefs_for_shrinkage)
}

#' @param ccs
#' @param perturbation_col
#' @param control_ids
#' @param nuisance_model_formula_str
#' @param cell_groups
#' @export
compare_genes_within_state_graph = function(ccs,
                                            perturbation_col = "perturbation", 
                                            control_ids = c("Control"), 
                                            nuisance_model_formula_str = "0", 
                                            cell_groups = NULL, 
                                            assembly_group = NULL, 
                                            state_graph = NULL,
                                            perturbations = NULL, 
                                            gene_ids = NULL,
                                            group_nodes_by = NULL,
                                            log_fc_thresh = 1,
                                            abs_expr_thresh = 1e-3,
                                            sig_thresh = 0.05,
                                            min_samples_detected = 2,
                                            min_cells_per_pseudobulk = 3,
                                            cores = 1,
                                            ...) {
  
  # to do make sure that ccs and state graph match 
  # check if all nodes in the state graph exist in the cds
  
  # if (is.null(state_graph) == FALSE) {
  #   
  # }
  
  expts = unique(colData(ccs)$expt)
  
  pb_cds = hooke:::pseudobulk_ccs_for_states(ccs, cell_agg_fun="sum")
  colData(pb_cds)$Size_Factor = colData(pb_cds)$num_cells_in_group
  
  # pb_cds = hooke:::add_covariate(ccs, pb_cds, perturbation_col)
  colData(pb_cds)[["perturbation"]] = colData(pb_cds)[[perturbation_col]] 
  colData(pb_cds)[["perturbation"]] = ifelse(colData(pb_cds)$perturbation %in% control_ids, 
                                             "Control", colData(pb_cds)$perturbation)
  
  # to do: if we want to run it by assembly group, must provide a state graph 
  if (is.null(assembly_group) == FALSE & is.null(state_graph) == FALSE) {
    pb_cds = hooke:::add_covariate(ccs, pb_cds, "assembly_group")
    pb_cds = pb_cds[, colData(pb_cds)$assembly_group == assembly_group]
    
    if (is(state_graph, "igraph")) {
      state_graph = igraph::as_data_frame(state_graph)
    }
    state_graph = state_graph %>% filter(assembly_group == assembly_group)
  }
  
  # ability to subset by perturbation 
  if (!is.null(perturbations)) {
    pb_cds = pb_cds[, colData(pb_cds)[[perturbation_col]] %in% c(control_ids, perturbations)]
  } 
  
  # subset to genes that are expressed over a certain min value
  expr_over_thresh = normalized_counts(pb_cds, "size_only", pseudocount = 0)
  genes_to_test = which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
  pb_cds = pb_cds[genes_to_test,]
  
  pseudobulks_to_test = which(colData(pb_cds)$num_cells_in_group >= min_cells_per_pseudobulk)
  pb_cds = pb_cds[,pseudobulks_to_test]
  
  
  if (is.null(cell_groups)) {
    cell_groups = rownames(counts(ccs))
  }
  
  df = data.frame(cell_group = cell_groups) %>% 
    mutate(genes_within_cell_group = purrr::map(.f = purrr:::possibly(compare_gene_expression_within_node,NA_real_),
                                                .x = cell_group, 
                                                pb_cds = pb_cds, 
                                                ccs = ccs,
                                                control_ids = c("Control"),
                                                cores = cores, 
                                                nuisance_model_formula_str = nuisance_model_formula_str))
  
  return(df)
  
}



extract_dispersion_helper = function(model, model_summary) {
  disp_res = if (class(model)[1] == "speedglm") {
    model_summary$dispersion
  } else if (class(model)[1] == "negbin"){
    model_summary$dispersion
  } else if (class(model) == "zeroinfl"){
    model_summary$dispersion
  } else {
    model_summary$dispersion
  }
}


update_summary <- function(model_tbl, dispersion_type=c("max", "fitted", "estimated"), min_dispersion=1){
  dispersion_type = match.arg(dispersion_type)
  model_tbl = model_tbl %>% mutate(
    disp_for_test = case_when(
      dispersion_type == "estimated" ~ dispersion,
      dispersion_type == "fitted" ~ disp_fit,
      dispersion_type == "max" ~ pmax(disp_fit, dispersion, min_dispersion),
      .default = dispersion),
    model_summary = purrr::map2(model, .y=disp_for_test, .f=summary))
  return(model_tbl)
}


compare_gene_expression_within_node <- function(cell_group, 
                                                ccs, 
                                                pb_cds,
                                                control_ids,
                                                nuisance_model_formula_str = "0",
                                                perturbation_ids = NULL,
                                                state_term ="cell_group",
                                                log_fc_thresh=1,
                                                abs_expr_thresh = 1e-3,
                                                sig_thresh=0.05, 
                                                cores=1) {
  
  # now fit models per cell group
  
  cg_pb_cds = pb_cds[, colData(pb_cds)[[state_term]] == cell_group]
  message(paste0("fitting regression models for ", cell_group))
  
  colData(cg_pb_cds)$Size_Factor = colData(cg_pb_cds)$num_cells_in_group
  
  full_model_str = paste0("~ 0 + perturbation")
  nuisance_model_formula_str = stringr::str_replace_all(nuisance_model_formula_str, "~", "")
  
  if (nuisance_model_formula_str != "0") {
    full_model_str = paste0(full_model_str, " + ", nuisance_model_formula_str)
  }
  
  pb_group_models = fit_models(cg_pb_cds,
                               model_formula_str = full_model_str,
                               cores=cores) %>% 
    dplyr::select(gene_short_name, id, model, model_summary, status)
  
  #pb_coeffs = collect_coefficients_for_limma(cg_pb_cds, pb_group_models, abs_expr_thresh) #coefficient_table(pb_group_models) %>%
  
  pb_coeffs = collect_coefficients_for_shrinkage(cg_pb_cds, pb_group_models, abs_expr_thresh, term_to_keep = "perturbation") #coefficient_table(pb_group_models) %>%
  
  estimates_for_controls = Matrix::rowSums(pb_coeffs$coefficients[,control_ids, drop=F])
  stderrs_for_controls = sqrt(Matrix::rowSums(pb_coeffs$stdev.unscaled[,control_ids, drop=F]^2))
  
  if (is.null(perturbation_ids)){
    perturbation_ids = setdiff(colnames(pb_coeffs$coefficients), control_ids)
  }else{
    perturbation_ids = intersect(perturbation_ids, colnames(pb_coeffs$coefficients))
  }
  
  #perturbation_ids = setdiff(colnames(estimate_matrix))
  cell_perturbations = tibble(term = unique(perturbation_ids))
  
  contrast_helper = function(perturb_id, 
                             PEM, 
                             PSEM, 
                             control_effects, 
                             control_effects_se){
    
    # ### Debugging:
    # to_include = control_effects > -10 & PEM[,perturb_id] > -10 
    # #to_include = row.names(PEM)
    # effect_est = PEM[to_include,perturb_id] - control_effects[to_include]
    # se_est = sqrt(PSEM[to_include,perturb_id]^2 + control_effects_se[to_include]^2)
    # #se_est[se_est > 2] = 2
    # #se_est[se_est < 0.5] = 0.5
    # shrunkren_res = ashr::ash(effect_est, se_est)
    # contrast_res = tibble(id = row.names(PEM[to_include,]),
    #                       raw_lfc = effect_est,
    #                       raw_lfc_se = se_est,
    #                       shrunken_lfc = shrunkren_res$result$PosteriorMean,
    #                       shrunken_lfc_se = shrunkren_res$result$PosteriorSD,
    #                       p_value = shrunkren_res$result$lfsr)
    # qplot(PSEM[to_include,perturb_id], control_effects_se[to_include], color=contrast_res$shrunken_lfc) + geom_abline()
    # qplot(raw_lfc, shrunken_lfc, data=contrast_res,color=contrast_res$p_value < 0.05) + geom_abline()
    # #qplot(effect_est, se_est)
    # ### end debugging
    
    
    effect_est = PEM[,perturb_id] - control_effects
    se_est = sqrt(PSEM[,perturb_id]^2 + control_effects_se^2)
    shrunkren_res = ashr::ash(effect_est, se_est,  method="fdr")
    
    contrast_res = tibble(id = row.names(PEM),
                          raw_lfc = effect_est,
                          raw_lfc_se = se_est,
                          shrunken_lfc = shrunkren_res$result$PosteriorMean,
                          shrunken_lfc_se = shrunkren_res$result$PosteriorSD,
                          p_value = shrunkren_res$result$lfsr)
    return(contrast_res)
  }
  #debug(contrast_helper)
  cell_perturbations = cell_perturbations %>% 
    mutate(perturb_effects = purrr:::map(.f = contrast_helper, 
                                         .x = term,
                                         control_effects = estimates_for_controls,
                                         control_effects_se = stderrs_for_controls,
                                         PEM = pb_coeffs$coefficients, 
                                         PSEM = pb_coeffs$stdev.unscaled))
  
  return(cell_perturbations) 
  
}


#' #' classify each gene's pattern of expression in each node of the state transition graph
#' #' @export
#' classify_genes_within_node <- function(cell_group, 
#'                                        ccs, 
#'                                        pb_cds, 
#'                                        state_term ="cell_group",
#'                                        log_fc_thresh=1,
#'                                        abs_expr_thresh = 1e-3,
#'                                        sig_thresh=0.05, 
#'                                        cores=1) {
#'   
#'   # now fit models per cell group
#'   
#'   cg_pb_cds = pb_cds[, colData(pb_cds)[[state_term]] == cell_group]
#'   message(paste0("fitting regression models for ", cell_group))
#'   
#'   
#'   pb_group_models = fit_models(cg_pb_cds,
#'                                model_formula_str=paste(paste0("~ 0 + perturbation")),
#'                                weights=colData(cg_pb_cds)$num_cells_in_group,
#'                                cores=cores) %>% 
#'     dplyr::select(gene_short_name, id, model, model_summary, status)
#'   
#'   
#'   pb_coeffs = shrunken_coefficient_table(cg_pb_cds, pb_group_models, abs_expr_thresh) #coefficient_table(pb_group_models) %>%
#'   dplyr::select(gene_short_name, id, term, estimate, std_err, p_value, status) %>%
#'     mutate(term = stringr::str_replace_all(term, "perturbation", ""))
#'   estimate_matrix = pb_coeffs %>% dplyr::select(id, term, estimate)
#'   estimate_matrix = estimate_matrix %>% mutate(term = factor(term, levels=unique(colData(cg_pb_cds)[,"perturbation"])))
#'   estimate_matrix = estimate_matrix %>% tidyr::pivot_wider(names_from=term, values_from=estimate, values_fill=0)
#'   
#'   gene_ids = estimate_matrix$id
#'   estimate_matrix$id = NULL
#'   estimate_matrix = as.matrix(estimate_matrix)
#'   row.names(estimate_matrix) = gene_ids
#'   colnames(estimate_matrix) = as.character(colnames(estimate_matrix))
#'   
#'   stderr_matrix = pb_coeffs %>% dplyr::select(id, term, std_err)
#'   stderr_matrix = stderr_matrix %>% mutate(term = factor(term, levels=unique(colData(cg_pb_cds)[,"perturbation"])))
#'   stderr_matrix = stderr_matrix %>% tidyr::pivot_wider(names_from=term, values_from=std_err, values_fill=0)
#'   
#'   gene_ids = stderr_matrix$id
#'   stderr_matrix$id = NULL
#'   stderr_matrix = as.matrix(stderr_matrix)
#'   row.names(stderr_matrix) = gene_ids
#'   colnames(stderr_matrix) = as.character(colnames(stderr_matrix))
#'   
#'   # collect the ids of any genes that threw an exception in fit_models and
#'   # set their estimates and std_errors to NA
#'   fail_gene_ids = pb_group_models %>% filter(status == "FAIL") %>% pull(id)
#'   if (length(fail_gene_ids) > 0){
#'     estimate_matrix[fail_gene_ids,] = log(abs_expr_thresh)
#'     stderr_matrix[fail_gene_ids,] = Inf
#'   }
#'   
#'   # Now shrink the std errors with limma
#'   # FIXME: this needs read degrees of freedom
#'   stderr_matrix[estimate_matrix < log(abs_expr_thresh)] = Inf
#'   estimate_matrix[estimate_matrix < log(abs_expr_thresh)] = log(abs_expr_thresh)
#'   
#'   # make a graph of control --> all perturbations
#'   cell_perturbations = tibble(perturbation = unique(colData(pb_cds)[,"perturbation"]))
#'   state_graph = data.frame("from" = cell_perturbations[cell_perturbations != "Control"])
#'   state_graph$to = "Control"
#'   # igraph defaults to first col > second col, so need to reverse the direction 
#'   state_graph = state_graph %>% igraph::graph_from_data_frame() %>% igraph::reverse_edges()
#'   
#'   cell_perturbations = cell_perturbations %>%
#'     dplyr::mutate(gene_classes = purrr::map(.f = purrr::possibly(
#'       classify_genes_in_cell_state, NA_real_), .x = perturbation,
#'       state_graph, estimate_matrix, stderr_matrix, state_term,
#'       log_fc_thresh=log_fc_thresh,
#'       abs_expr_thresh = abs_expr_thresh,
#'       sig_thresh=sig_thresh,
#'       cores=cores))
#'   
#'   cell_perturbations = cell_perturbations %>%
#'     filter(is.na(gene_classes) == FALSE) %>%
#'     dplyr::mutate(gene_class_scores = purrr::map2(.f = purrr::possibly(
#'       score_genes_for_expression_pattern, NA_real_),
#'       .x = perturbation,
#'       .y = gene_classes,
#'       state_graph,
#'       estimate_matrix))
#'   
#'   # cell_perturbations$cell_group = cell_group
#'   
#'   cell_perturbations = cell_perturbations %>%
#'     filter(is.na(gene_classes) == FALSE) %>%
#'     tidyr::unnest(gene_class_scores) %>% 
#'     dplyr::select(perturbation, gene_id, interpretation, pattern_match_score, pattern_activity_score) #%>% 
#'   # dplyr::filter(!interpretation %in% c("Absent", "Maintained", "Specifically maintained", "Selectively maintained"))
#'   
#'   
#'   cell_perturbations = left_join(cell_perturbations,
#'                                  rowData(ccs@cds) %>%
#'                                    as_tibble %>%
#'                                    select(id, gene_short_name), 
#'                                  by=c("gene_id"="id")) %>% 
#'     mutate(group = case_when(
#'       grepl(pattern ="down", interpretation) | grepl(pattern = "de", interpretation)  ~ "Down",
#'       grepl(pattern ="aintain", interpretation) ~ "Maintained",
#'       T ~ "Up"
#'     )) %>% 
#'     mutate(broad_interpretation = case_when(
#'       (grepl(pattern = "ly", interpretation) & group == "Up")  ~ "Sp/Sel Up",
#'       (grepl(pattern ="ly", interpretation) & group == "Down") ~ "Sp/Sel Down",
#'       (!grepl(pattern ="ly", interpretation) & group == "Up") ~ "Up",
#'       (!grepl(pattern ="ly", interpretation) & group == "Down") ~ "Down",
#'       (grepl(pattern = "ly", interpretation) & group == "Maintained")  ~ "Sp/Sel Maintained",
#'       (!grepl(pattern ="ly", interpretation) & group == "Maintained") ~ "Maintained",
#'     )) %>% 
#'     select(-group)
#'   
#'   cell_perturbations = left_join(cell_perturbations, pb_coeffs, 
#'                                  by = c("perturbation" = "term", 
#'                                         "gene_short_name", 
#'                                         "gene_id" = "id"))
#'   
#'   
#'   return(cell_perturbations) 
#'   
#' }


calc_gsea_enrichment_on_state_specific_genes <- function(gene_patterns_within_state_graph, 
                                                         gene_df, 
                                                         sig_thresh = 0.1) {
  
  gene_set_list = split(x = gene_df$gene_short_name, f = gene_df$gs_name)
  gene_ranking = gene_patterns_within_state_graph$pattern_activity_score[,1]
  names(gene_ranking) = gene_patterns_within_state_graph %>% pull(gene_short_name)
  gsea_res = fgsea::fgsea(pathways=gene_set_list, stats=gene_ranking) %>% as_tibble()
  gsea_res = gsea_res %>% filter(padj < sig_thresh)
  return(gsea_res)
}


calc_pathway_enrichment_on_state_specific_genes <- function(gene_df, msigdbr_t2g, sig_thresh = 0.1, ...){
  gene_symbols_vector = gene_df$gene_short_name
  enrich_res = clusterProfiler::enricher(gene = gene_symbols_vector, TERM2GENE = msigdbr_t2g, ...) %>% as_tibble()
  return(enrich_res)
}


#' Function to find degs within a cell type 
#' @param ccm output from fit_genotype_ccm 
#' @param ... filtering criteria for
fit_genotype_deg = function(ccm, 
                            cores = 1, 
                            ...) {
  
  # filter by cell types 
  sub_ccs = subset_ccs(ccm@ccs, ...)
  sub_pb_cds = pseudobulk_ccs_for_states(sub_ccs)
  sub_pb_cds = add_covariate(sub_ccs, 
                             sub_pb_cds, 
                             ccm@info$perturbation_col)
  # filter by perturbation 
  sub_pb_cds = sub_pb_cds[, colData(sub_pb_cds)[[ccm@info$perturbation_col]] %in% c(ccm@info$genotype, ccm@info$ctrl_ids)]
  colData(sub_pb_cds)$knockout = ifelse(colData(sub_pb_cds)[[ccm@info$perturbation_col]] %in% ccm@info$ctrl_ids, F, T)
  
  gene_fits <- fit_models(sub_pb_cds, 
                          model_formula_str = "~ knockout", 
                          weights=colData(sub_pb_cds)$num_cells_in_group, 
                          cores = cores)
  
  return(gene_fits)
  
}


compare_genes_in_cell_state <- function(cell_state, 
                                        state_graph, 
                                        estimate_matrix, 
                                        stderr_matrix,
                                        state_term="cell_group", log_fc_thresh=1, 
                                        abs_expr_thresh = 1e-3, sig_thresh=0.05, cores=1) {
  
  
  #expr_self = expr_mat[,cell_state]
  
  parents = get_parents(state_graph, cell_state) #igraph::neighbors(state_graph, cell_state, mode="in")
  parents = intersect(parents, colnames(estimate_matrix))
  
  children = get_children(state_graph, cell_state)#igraph::neighbors(state_graph, cell_state, mode="out")
  children = intersect(children, colnames(estimate_matrix))
  
  siblings = get_siblings(state_graph, cell_state)#igraph::neighbors(state_graph, parents, mode="out")
  siblings = intersect(siblings, colnames(estimate_matrix))
  
  states_in_contrast = c(cell_state, parents, children, siblings) %>% unique()
  
  expr_df = tibble(gene_id=row.names(estimate_matrix))
  
  message("      examining coefficients ", cell_state)
  
  expr_df$expr_self = pnorm(estimate_matrix[,cell_state] - log(abs_expr_thresh), sd = stderr_matrix[,cell_state], lower.tail=FALSE)
  expr_df$expr_self = p.adjust(expr_df$expr_self, method="BH") < sig_thresh
  
  expr_df$expressed_in_parents = NA
  expr_df$expressed_in_siblings = NA
  expr_df$higher_than_parents = NA
  expr_df$lower_than_parents = NA
  expr_df$higher_than_all_siblings = NA
  expr_df$lower_than_all_siblings = NA
  expr_df$higher_than_siblings = NA
  expr_df$lower_than_siblings = NA
  expr_df$expressed_in_children = NA
  expr_df$higher_than_all_children = NA
  expr_df$lower_than_all_children = NA
  expr_df$higher_than_children = NA
  expr_df$lower_than_children = NA
  
  
  contrast_helper = function(parents, 
                             children, 
                             PEM, 
                             PSEM){
    
    
    parent_effects = Matrix::rowSums(PEM[,parents, drop=F])
    parent_effects_se = sqrt(Matrix::rowSums(PSEM[,parents, drop=F]^2))
    
    child_effects = Matrix::rowSums(PEM[,children, drop=F])
    child_effects_se = sqrt(Matrix::rowSums(PSEM[,children, drop=F]^2))
    
    # effect_est = PEM[,perturb_id] - parent_effects
    # se_est = sqrt(PSEM[,perturb_id]^2 + parent_effects_se^2)
    effect_est = child_effects - parent_effects
    se_est = sqrt(child_effects_se^2 + parent_effects_se^2)
    shrunkren_res = ashr::ash(effect_est, se_est,  method="fdr")
    
    contrast_res = tibble(id = row.names(PEM),
                          raw_lfc = effect_est,
                          raw_lfc_se = se_est,
                          shrunken_lfc = shrunkren_res$result$PosteriorMean,
                          shrunken_lfc_se = shrunkren_res$result$PosteriorSD,
                          p_value = shrunkren_res$result$lfsr)
    return(contrast_res)
  }
  
  
  if (length(parents) > 0){
    
    
    expressed_in_parents_mat = pnorm(estimate_matrix[,parents, drop=F] - log(abs_expr_thresh), sd = stderr_matrix[,parents, drop=F], lower.tail=FALSE)
    expressed_in_parents_mat = apply(expressed_in_parents_mat, 2, p.adjust, method="BH")
    
    expressed_in_parents_mat = expressed_in_parents_mat < sig_thresh
    expr_df$expressed_in_parents = Matrix::rowSums(expressed_in_parents_mat) > 0
    
    # higher_than_parents_stat = -t( sweep( t(estimate_matrix[,parents, drop=F]), 2, as.numeric(estimate_matrix[,cell_state]) , `-`) )
    # higher_than_parents_pval = pnorm(higher_than_parents_stat,
    #                                  sd = sqrt(sweep(t(stderr_matrix[,parents, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    
    
    
    # higher_than_parents_stat = new_estimate_matrix[,cell_state]
    # higher_than_parents_pval = p_value_matrix[, cell_state]
    res = contrast_helper(parents, children, estimate_matrix, stderr_matrix)
    higher_than_parents_stat = res[, "shrunken_lfc"]
    higher_than_parents_pval = res[,"p_value"]
    higher_than_parents_pval = apply(higher_than_parents_pval, 2, p.adjust, method="BH")
    
    higher_than_parents_mat = abs(higher_than_parents_stat) > log_fc_thresh & higher_than_parents_pval < sig_thresh
    expr_df$higher_than_parents = Matrix::rowSums(higher_than_parents_mat) > 0
    
    # lower_than_parents_pval = pnorm(-higher_than_parents_stat,
    #                                 sd = sqrt(sweep(t(stderr_matrix[,parents, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    # lower_than_parents_pval = p_value_matrix[,cell_state]
    lower_than_parents_pval = res[,"p_value"]
    lower_than_parents_pval = apply(lower_than_parents_pval, 2, p.adjust, method="BH")
    
    lower_than_parents_mat = abs(higher_than_parents_stat) > log_fc_thresh & lower_than_parents_pval < sig_thresh
    expr_df$lower_than_parents = Matrix::rowSums(lower_than_parents_mat) > 0
  }else{
    expr_df$expressed_in_parents = NA
    expr_df$expressed_in_siblings = NA
    expr_df$higher_than_parents = NA
    expr_df$lower_than_parents = NA
    expr_df$higher_than_all_siblings = NA
    expr_df$lower_than_all_siblings = NA
    expr_df$higher_than_siblings = NA
    expr_df$lower_than_siblings = NA
    expr_df$expressed_in_children = NA
    expr_df$higher_than_all_children = NA
    expr_df$lower_than_all_children = NA
    expr_df$higher_than_children = NA
    expr_df$lower_than_children = NA
  }
  
  if (length(siblings) > 0){
    
    res = contrast_helper(cell_state, siblings, estimate_matrix, stderr_matrix)
    
    expressed_in_siblings_mat = pnorm(estimate_matrix[,siblings, drop=F] - log(abs_expr_thresh), sd = stderr_matrix[,siblings, drop=F], lower.tail=FALSE)
    expressed_in_siblings_mat = apply(expressed_in_siblings_mat, 2, p.adjust, method="BH")
    
    expressed_in_siblings_mat = expressed_in_siblings_mat < sig_thresh
    expr_df$expressed_in_siblings = Matrix::rowSums(expressed_in_siblings_mat) > 0
    
    # higher_than_siblings_stat = -t(sweep(t(estimate_matrix[,siblings, drop=F]), 2, as.numeric(estimate_matrix[,cell_state]) , `-`))
    # higher_than_siblings_pval = pnorm(higher_than_siblings_stat,
    #                                   sd = sqrt(sweep(t(stderr_matrix[,siblings, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    
    higher_than_siblings_stat = res[, "shrunken_lfc"]
    higher_than_siblings_pval = res[,"p_value"]
    higher_than_siblings_pval = apply(higher_than_siblings_pval, 2, p.adjust, method="BH")
    
    higher_than_siblings_mat = abs(higher_than_siblings_stat) > log_fc_thresh & higher_than_siblings_pval < sig_thresh
    expr_df$higher_than_all_siblings = Matrix::rowSums(higher_than_siblings_mat) == ncol(higher_than_siblings_pval)
    expr_df$higher_than_siblings = Matrix::rowSums(higher_than_siblings_mat) > 0
    
    # lower_than_siblings_pval = pnorm(-higher_than_siblings_stat,
    #                                  sd = sqrt(sweep(t(stderr_matrix[,siblings, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    lower_than_siblings_pval = res[,"p_value"]
    lower_than_siblings_pval = apply(lower_than_siblings_pval, 2, p.adjust, method="BH")
    
    lower_than_siblings_mat = abs(higher_than_siblings_stat) > log_fc_thresh & lower_than_siblings_pval < sig_thresh
    expr_df$lower_than_all_siblings = Matrix::rowSums(lower_than_siblings_mat) == ncol(lower_than_siblings_mat)
    expr_df$lower_than_siblings = Matrix::rowSums(lower_than_siblings_mat) > 0
    
    
  }else{
    expr_df$expressed_in_siblings = NA
    expr_df$higher_than_all_siblings = NA
    expr_df$lower_than_all_siblings = NA
    expr_df$higher_than_siblings = NA
    expr_df$lower_than_siblings = NA
  }
  
  if (length(children) > 0){
    
    res = contrast_helper(cell_state, children, estimate_matrix, stderr_matrix)
    
    
    expressed_in_children_mat = pnorm(estimate_matrix[,children, drop=F] - log(abs_expr_thresh), sd = stderr_matrix[,children, drop=F], lower.tail=FALSE)
    expressed_in_children_mat = apply(expressed_in_children_mat, 2, p.adjust, method="BH")
    
    expressed_in_children_mat = expressed_in_children_mat < sig_thresh
    expr_df$expressed_in_children = Matrix::rowSums(expressed_in_children_mat) > 0
    
    # higher_than_children_stat = -t(sweep(t(estimate_matrix[,children, drop=F]), 2, as.numeric(estimate_matrix[,cell_state]) , `-`))
    # higher_than_children_pval = pnorm(higher_than_children_stat,
    #                                   sd = sqrt(sweep(t(stderr_matrix[,children, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    higher_than_children_stat = res[, "shrunken_lfc"]
    higher_than_children_pval = res[,"p_value"]
    higher_than_children_pval = apply(higher_than_children_pval, 2, p.adjust, method="BH")
    
    higher_than_children_mat = abs(higher_than_children_stat) > log_fc_thresh & higher_than_children_pval < sig_thresh
    expr_df$higher_than_all_children = Matrix::rowSums(higher_than_children_mat) == ncol(higher_than_children_pval)
    expr_df$higher_than_children = Matrix::rowSums(higher_than_children_mat) > 0
    
    # lower_than_children_pval = pnorm(-higher_than_children_stat,
    # sd = sqrt(sweep(t(stderr_matrix[,children, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    lower_than_children_pval = res[,"p_value"]
    # lower_than_children_pval = apply(lower_than_children_pval, 2, p.adjust, method="BH")
    
    lower_than_children_mat = abs(higher_than_children_stat) > log_fc_thresh & lower_than_children_pval < sig_thresh
    expr_df$lower_than_all_children = Matrix::rowSums(lower_than_children_mat) == ncol(lower_than_children_mat)
    expr_df$lower_than_children = Matrix::rowSums(lower_than_children_mat) > 0
    
    
  }else{
    expr_df$expressed_in_children = NA
    expr_df$higher_than_all_children = NA
    expr_df$lower_than_all_children = NA
    expr_df$higher_than_children = NA
    expr_df$lower_than_children = NA
  }
  
  expr_df = expr_df %>% tidyr::nest(data = !gene_id)
  
  message("      interpreting patterns")
  interpret_expression_pattern = function(pat_df){
    if (pat_df$expr_self){
      if (is.na(pat_df$expressed_in_parents)){
        # no parents, therefore no siblings
        #return ("Maintained")
        if (is.na(pat_df$expressed_in_children)){
          return("Maintained")
        } else {
          # no parent, but there are children
          if (pat_df$expressed_in_children == FALSE | pat_df$higher_than_all_children){
            # Higher than parent, and higher than children
            return("Precursor-specific")
          }
          else if (pat_df$higher_than_children){
            # no parent, higher than children
            return("Precursor-specific")
          }
          else if(pat_df$lower_than_children){
            # no parent, but lower than children
            return("Precursor-depleted")
          }
          else { # no parent same as children
            return("Maintained")
          }
        }
      }else if (pat_df$expressed_in_parents){
        # Expressed in self and parent
        if (is.na(pat_df$expressed_in_siblings)){
          # Expressed in self and parent and there are no siblings
          if (pat_df$higher_than_parents){
            if (is.na(pat_df$expressed_in_children)){
              return("Upregulated")
            } else {
              # there are children
              if (pat_df$expressed_in_children == FALSE | pat_df$higher_than_all_children){
                # Higher than parent, and higher than siblings
                return("Transiently upregulated")
              }
              else if(pat_df$lower_than_all_children){
                # lower than children
                return("Increasingly upregulated")
              }
              else { # same as children
                return("Upregulated")
              }
            }
          }
          else if(pat_df$lower_than_parents){
            if (is.na(pat_df$expressed_in_children)){
              return("Downregulated")
            } else {
              # there are children
              if (pat_df$lower_than_all_children){
                # Lower than parent, and lower than children
                return("Decreasingly downregulated")
              }
              else if(pat_df$lower_than_all_children){
                # lower than children
                return("Transiently downregulated")
              }
              else { # same as children
                return("Downregulated")
              }
            }
          }else{
            if (is.na(pat_df$expressed_in_children)){
              return("Maintained")
            } else {
              # same as parent, and there are children
              if (pat_df$expressed_in_children == FALSE | pat_df$higher_than_all_children){
                # Higher than parent, and higher than children
                return("Precursor-specific")
              }
              else if(pat_df$lower_than_all_children){
                # no parent, but lower than children
                return("Precursor-depleted")
              }
              else { # no parent same as children
                return("Maintained")
              }
            }
          }
        } else {
          # Expressed in self and parent and there are siblings
          if (pat_df$higher_than_parents){
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Higher than parent, and higher than siblings
              return("Specifically upregulated")
            }
            else if (pat_df$higher_than_siblings){
              # Higher than parent, and higher than siblings
              return("Selectively upregulated")
            }
            else if(pat_df$lower_than_siblings){
              # Higher than parent, but lower than siblings
              return("Upregulated")
            }
            else { # higher than parent, same as siblings
              return("Upregulated")
            }
          }
          else if(pat_df$lower_than_parents){
            if (pat_df$expressed_in_siblings == TRUE & pat_df$lower_than_all_siblings){
              # Lower than parent, and higher than siblings
              return("Specifically downregulated")
            }
            else if (pat_df$expressed_in_siblings == TRUE & pat_df$lower_than_siblings){
              # Lower than parent, and higher than some siblings
              return("Selectively downregulated")
            }
            else { # lower than parent, same as or higher than siblings
              return("Downregulated")
            }
          }
          else { # same as parent
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Same as parent, and higher than all siblings
              return("Specifically maintained")
            }
            else if (pat_df$higher_than_siblings){
              # Same as parent, and higher than some siblings
              return("Selectively maintained")
            }
            else if(pat_df$lower_than_all_siblings){
              # Same as parent, but lower than siblings
              return("Maintained")
            }
            else { # same as parent, same as siblings
              return("Maintained")
            }
          }
        }
        
      }else{
        # expressed in self but not in parent
        if (is.na(pat_df$expressed_in_siblings)){
          # Expressed in self, not in parent and there are no siblings
          if (pat_df$higher_than_parents)
            return("Activated")
          else if(pat_df$lower_than_parents)
            return("Downregulated") # shouldn't happen
          else
            return("Activated") # might happen if its above threshold but not significantly above parent (and parent is below thresh)
        } else {
          # expressed in self, not in parent, and there are siblings
          if (pat_df$higher_than_parents){ # Expressed in self and higher than parent and there are siblings
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Higher than parent, and higher than all siblings
              return("Specifically activated")
            } else if (pat_df$higher_than_siblings){
              # Higher than parent, and higher than some siblings
              return("Selectively activated")
            }
            else if(pat_df$lower_than_all_siblings){
              # Higher than parent (which is off), but lower than all siblings
              return("Activated")
            }
            if(pat_df$lower_than_siblings){
              # Higher than parent (which is off), but lower than some siblings
              return("Activated")
            }
            else { # Higher than parent (which is off), same as siblings
              return("Activated")
            }
          }
          else if(pat_df$lower_than_parents){
            # gene is expressed, lower in the parent (which is off)
            if (pat_df$higher_than_all_siblings){
              # Lower than parent, and higher than all siblings
              return("Absent") # shouldn't happen
            }
            else if (pat_df$higher_than_siblings){
              # Lower than parent, and higher than some siblings
              return("Absent") # shouldn't happen
            }
            else if(pat_df$lower_than_all_siblings){
              # Lower than parent and  lower than all siblings
              return("Absent")
            }
            else if(pat_df$lower_than_siblings){
              # Lower than parent and  lower than some siblings
              return("Absent")
            }
            else { # Lower than parent, same as siblings
              return("Absent")
            }
          }
          else { # same as parent (which is off)
            if (pat_df$higher_than_all_siblings){
              # Same as parent, and higher than all siblings
              return("Absent")
            }
            else if (pat_df$higher_than_siblings){
              # Same as parent, and higher than some siblings
              return("Absent")
            }
            else if(pat_df$lower_than_all_siblings){
              # Same as parent, but lower than all siblings
              return("Absent")
            }
            else if(pat_df$lower_than_siblings){
              # Same as parent, but lower than some siblings
              return("Absent")
            }
            else { # same as parent, same as siblings
              return("Absent")
            }
          }
        }
      }
      return ("Expressed")
    }else{
      # Not expressed in self
      if (is.na(pat_df$expressed_in_parents)){
        # no parents, therefore no siblings
        return ("Absent")
      }else if (pat_df$expressed_in_parents){
        # Not expressed in self, but expressed in parents
        if (is.na(pat_df$expressed_in_siblings)){
          # Not expressed in self, expressed parent and there are no siblings
          if(pat_df$lower_than_parents)
            return("Deactivated")
          else
            return("Absent") # shouldn't happen
        } else {
          # Not expressed in self, expressed in parent and there are siblings
          if(pat_df$lower_than_parents){
            # Lower than parent
            if(pat_df$lower_than_all_siblings){
              # Lower than parent and  lower than siblings
              return("Specifically deactivated")
            }
            else if(pat_df$lower_than_siblings){
              # Lower than parent and  lower than siblings
              return("Selectively deactivated")
            }
            return("Deactivated")
          }
          else {
            #Not expressed in self, not lower than parent
            return ("Absent")
          }
        }
      }else{
        # Not expressed in self or parents
        return ("Absent")
      }
      return ("Absent")
    }
    return ("Absent")
    #match_row = match(data.frame(t(pat_df)), data.frame(t(interp_table)))
    #interpetation[match_row]
  }
  #debug(interpret_expression_pattern)
  expr_df = expr_df %>% mutate(interpretation = purrr::map(.f = purrr::possibly(
    interpret_expression_pattern, NA_character_), .x = data))
  expr_df = expr_df %>% tidyr::unnest(interpretation)
  message("      completed ", cell_state)
  return(expr_df)
}
