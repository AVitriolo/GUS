connect_KNN_igraph <- function(KNN) {

  null_idxs <- which(vapply(KNN, is.null, logical(1)))
  KNN[null_idxs] <- null_idxs

  n <- length(KNN)

  # Original checks ALL pairs (i,j) -- O(n^2) intersect() calls, which is
  # what's costing ~25s per grid point. But KNN[[i]] only ever contains
  # spatial neighbors (indices within a small offset of i, per
  # find_neighbors()'s K), so KNN[[i]] and KNN[[j]] can only overlap if i
  # and j are within that same small offset of each other -- anything
  # farther apart is guaranteed disjoint. WINDOW is computed directly from
  # the actual data (max index offset seen in any KNN[[i]]) rather than
  # hardcoded, so this stays correct even if K changes upstream.
  offsets <- unlist(lapply(seq_len(n), function(i) KNN[[i]] - i))
  max_offset <- if (length(offsets) > 0) max(abs(offsets)) else 0
  WINDOW <- max(1, 2 * max_offset)

  edges <- integer(0)
  for (i in seq_len(n - 1)) {
    j_end <- min(i + WINDOW, n)
    for (j in (i + 1):j_end) {
      if (length(intersect(KNN[[i]], KNN[[j]])) > 0) {
        edges <- c(edges, i, j)
      }
    }
  }

  g <- igraph::graph(edges = edges, directed = FALSE, n = n)
  comps <- igraph::components(g)

  clusters_idxs <- split(x = seq_along(KNN), f = comps$membership)

  KNN.connected <- lapply(clusters_idxs, function(idx){sort(unique(unlist(KNN[idx])))})

  return(list(
    clusters   = KNN.connected,
    graph      = g,
    membership = comps$membership
  ))
}