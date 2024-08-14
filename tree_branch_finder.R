library("phytools")


find_branch_number <- function(start_node, end_node, edge_matrix) {
  
  edge <- which(edge_matrix[,1] == start_node & edge_matrix[,2] == end_node)
  
    return(edge)
  
}

