# the problem of this function right now is that 
# it fails in merging together clusters sharing elements that are not contiguous
connect_KNN <- function(KNN){
    
    null_idxs <- which(unlist(lapply(KNN, is.null)))
    KNN[null_idxs] <- null_idxs

    i <- 1
    j <- i + 1
    k <- i

    cluster <- KNN[[k]]
    cluster_name <- i

    KNN.connected <- list()
    while (j < length(KNN)){
        
        while(j <= length(KNN) && length(intersect(cluster, KNN[[j]])) > 0){
            cluster <- sort(unique(unlist(KNN[k:j])))
            i <- i + 1
            j <- i + 1
        }

        KNN.connected[[cluster_name]] <-  cluster
        
        k <- j
        i <- k
        j <- i + 1

        cluster <- KNN[[min(k, length(KNN))]]
        cluster_name <- cluster_name + 1
        
    }

    KNN.connected[[cluster_name]] <-  cluster
    
    return(KNN.connected)
}