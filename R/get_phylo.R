get_edge_list_by_clade <- function (epidemic, split_by=NULL) {
  epidemic <- data.frame(epidemic)
  if(is.null(split_by)) split_by <- factor(epidemic$parents)
  edges <- lapply(split(epidemic, split_by), function (x) {
    parent_id <- x[1, 1]
    if (parent_id == -1) {
      return(NULL)
    }
    new_order <- order(x$infection_time, decreasing=TRUE)
    id_new_order <- as.numeric(rownames(x)[new_order])-1
    parent_inf_time <- as.numeric(epidemic[parent_id+1, 2])
    n <- (nrow(x)+1)# Number of lineages in current clade (number of children + founding lineage)
    # Assumes that the parent lineage is sampled after giving birth to the last offspring
    edge_matrix <- cbind(c(-1:-n, -2:-n),
                         c(-2:-n, parent_id, id_new_order[(2:n)-1]))
    num_vertices <- nrow(x)*2+1
    vertex.times <- cbind(c(-1:-n, parent_id, id_new_order),
                          c(parent_inf_time,
                            rev(x$infection_time[new_order]),
                            epidemic$recovery_time[parent_id+1],
                            rev(x$recovery_time[new_order])))
    edge_lengths <- apply(edge_matrix, 1, function (pair) {
      time.diff <- vertex.times[vertex.times[, 1] == pair[2], 2] -
        vertex.times[vertex.times[, 1] == pair[1], 2]
      if (time.diff < 0) browser()
      return(time.diff)
    })
    if (any(is.na(edge_lengths))) browser()
    if (parent_id==0) {
      edge_matrix <- edge_matrix[-1, ]
      edge_matrix <- apply(edge_matrix, 2, function (COL) {
        sapply(COL, function (ELE) ifelse(ELE>=0, ELE, ELE+1))
      })
      edge_lengths <- edge_lengths[-1]
    }
    return(list(edge_matrix=edge_matrix, edge_lengths=edge_lengths))
  })
  return (edges)
}

summarise_clade_edges <- function (by_clade) {
  do.call(rbind, lapply(seq_along(by_clade), function (i) {
    if (is.null(by_clade[[i]])) return(NULL)
    MATRIX <- by_clade[[i]]$edge_matrix
    matrix(sapply(c(MATRIX), function(x) {
      if (x < 0) return(paste0(i, "_", (-1*x)))
      return(x)
    }), ncol=2)
  }))
}

adjust_nodes <- function (parents, edge_mat) {
  do.call(rbind, lapply(parents, function (parent_id) {
    id <- which(edge_mat[, 2] == as.character(parent_id)) # row number of original lineage to parent
    ptr <- max(id)  # row number of final lineage to parent
    # Find row number of duplicate node
    ptr2 <- which(edge_mat[, 2]==edge_mat[ptr, 1])
    while (length(ptr2) > 0) {
      ptr <- ptr2
      ptr2 <- which(edge_mat[, 2]==edge_mat[ptr, 1])
    }
    return(c(min(id), # edge to delete
             ptr, # node to collapse
             edge_mat[min(id), 1]) # rename node to this
    )
  }))
}

create_tree_from_edges <- function (ntip, edge_mat, edge_lengths, replace_nodes, tip.labels=NULL) {
  nnode <- ntip - 1
  tr <- rtree(ntip)
  tr$edge <- edge_mat
  tr$edge[as.numeric(replace_nodes[, 2]), 1] <- replace_nodes[, 3]
  tr$edge <- tr$edge[-as.numeric(replace_nodes[, 1]), ]
  tr$edge.length <- as.numeric(edge_lengths[-as.numeric(replace_nodes[, 1])])
  node_map <- matrix(unique(grep("_", tr$edge[, 1], value=TRUE)))
  node_map <- cbind(node_map, 1:length(node_map)+ntip)
  tr$edge <- matrix(sapply(c(tr$edge), function (x) {
    if (grepl("_", x)) return(as.numeric(node_map[node_map[, 1]==x, 2]))
    return(as.numeric(x))
  }), ncol=2)
  tr$edge[, 2] <- sapply(tr$edge[, 2], function (x) ifelse(x<ntip, x+1, x))
  if (!is.null(tip.labels)) tr$tip.label <- tip.labels
  return(tr)
}

#' Generate the phylogenetic tree for an epidemic based on who infected whom
#'
#' @param epidemic A matrix detailing who infected whom at what time
#'
#' @return
#' @export
#'
#' @examples
get_phylo <- function (epidemic) {
  ##### DO NOT PLOT THE RESULTING TREE IF MORE THAN ONE ROOT
  epidemic <- data.frame(epidemic)
  clades <- split(epidemic, factor(epidemic$parents))
  # Divide edge list into clades (tree formed by a parent and its offspring)
  by_clade <- get_edge_list_by_clade(epidemic)
  edge_mat <- summarise_clade_edges(by_clade)
  edge_lengths <- unlist(lapply(by_clade, function (x) x$edge_lengths))
  replace_nodes <- adjust_nodes(unique(epidemic$parents)[-1:-2], edge_mat)
  ntip <- nrow(epidemic)
  tree <- create_tree_from_edges(ntip, edge_mat, edge_lengths, replace_nodes,
                                 paste0(1:nrow(epidemic)-1, "_", epidemic$recovery_time))
  return(tree)
}

