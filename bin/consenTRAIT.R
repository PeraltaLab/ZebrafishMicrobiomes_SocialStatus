################################################################################
#                                                                              #
#   consenTRAIT function                                                       #
#                                                                              #
#   Originally by: Adam Martiny                                                #
#   Minor Modification by: Mario Muscarella                                    #
#                                                                              #
#   Last Update: 2015-02-18                                                    #
#                                                                              #
#   Citation: Martiny et al. 2013. Phylogenetic conservatism of functional     #
#               traits in microorganisms. ISME J. 7(4): 830–838                #
#                                                                              #
################################################################################

require("data.table")||install.packages("data.table");require("data.table")
require("adephylo")||install.packages("adephylo");require("adephylo")
require("ape")||install.packages("ape");require("ape")

consenTRAIT <- function(data = "", phy = "", cutoff = 0.9){
  if (class(phy) != "phylo"){
    stop("tree must be a phylo class object")
  }
  if(is.rooted(phy) == FALSE){
    stop("tree must be rooted")
  }
  tree <- phy
  table <- data
  sub.trees <- subtrees(tree)
  tip.dists <- distRoot(tree, method = "patristic")
  node.dists <- as.matrix(dist.nodes(tree))[,length(tree$tip.label) + 1]
  
  # Output File
  ConsenTrait <- as.data.frame(matrix(NA, nrow = dim(table)[2], ncol = 7))
  colnames(ConsenTrait) <- c("Trait", "Clusters", "Cluster Size", 
                             "Distance", "CI_95", "Distance.L", "CI_95.L")
  
  #Starting ConsenTrait Algorithm
  cutoff = cutoff
  
  for(i in 1:dim(table)[2]){
    print(paste("Analyzing Trait ", i, " of ", dim(table)[2], ": ",
                colnames(table)[i], sep = ""), quote = F)
    
    # Define Inputs
    trait <- colnames(table)[i]
    genomes <- rownames(table)
    temp <- table[, i]
    
    # Initiate Temp Results
    origins <- vector(mode = "character", length = 0)
    positives <- vector(mode = "character", length = 0)
    cluster_size <- vector(mode = "numeric", length = 0)
    cluster_dist <- vector(mode = "list", length = 0)
    
    # Loop through all subtrees and determine if any subtrees have > 90% positives
    for (j in 1:length(sub.trees)){
      tip_names <- sub.trees[[j]]$tip.label
      tree_name <- sub.trees[[j]]$name
      if(mean(temp[which(genomes %in% tip_names)]) > cutoff ) {
        match_test <- match(tip_names, positives)
        if(all(is.na(match_test))){
          positives <- c(positives, tip_names)
          origins <- c(origins, tree_name)
          cluster_size <- c(cluster_size, length(sub.trees[[j]]$tip.label))
          temp_dist <- as.vector(tip.dists[which(tree$tip.label %in% tip_names)] - 
                                   node.dists[length(tree$tip.label) + j])
          cluster_dist <- c(cluster_dist, list(temp_dist))
        }
      }
    }
    cluster_dist_mean <- unlist(lapply(cluster_dist, mean))
    cluster_number <- length(cluster_size)
    cluster_size <- mean(cluster_size)
    cluster_mean_dist <- mean(cluster_dist_mean)
    cluster_ci_dist <- ci(cluster_dist_mean)
    cluster_mean_dist_log <- mean(log10(cluster_dist_mean))
    cluster_ci_dist_log <- ci(log10(cluster_dist_mean))
    
    # Save to output
    ConsenTrait[i, 1] <- trait
    ConsenTrait[i, 2] <- cluster_number
    ConsenTrait[i, 3] <- cluster_size
    ConsenTrait[i, 4] <- cluster_mean_dist
    ConsenTrait[i, 5] <- cluster_ci_dist
    ConsenTrait[i, 6] <- cluster_mean_dist_log
    ConsenTrait[i, 7] <- cluster_ci_dist_log
  }
  return(ConsenTrait)
}


  