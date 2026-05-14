find_neighbors <- function(data, idx, K = 1, max_distance = 1000, return_self = T){

	if (!("chr" %in% colnames(data))){
		data$chr <- as.integer(gsub(pattern = "_[0-9]+", replacement = "", x = data$CpGID))
	}
	
	if (!("start" %in% colnames(data))){
		data$start <- as.integer(gsub(pattern = "chr[0-9]+_", replacement = "", x = data$CpGID))
	}

    n <- nrow(data)

	if (!(return_self)){neighbors <- c()}
	else if ((return_self)){neighbors <- c(idx)}    

    for (k in 1:K){

        prev_row <- data[max(0, (idx-k)),]
		curr_row <- data[idx,]
		next_row <- data[(idx+k),]
	
		prev_chr <- prev_row$chr
		curr_chr <- curr_row$chr
		next_chr <- next_row$chr
	
	if(length(prev_chr) == 0){prev_chr <- "chrNA"; prev_bin_number <- NA}
	if(is.na(next_chr)){next_chr <- "chrNA"; next_bin_number <- NA}
	
	if((idx > 1) & (prev_chr == curr_chr)){
	  neighbors <- c(neighbors, as.integer(idx-k))
	}
	
	if((idx < n) & (next_chr == curr_chr)){
	  neighbors <- c(neighbors, as.integer(idx+k))
	}

	# add filtering for max_distance
	# compute distance of all neighborhood from the central CpG
	# filter out if dist > max_distance
	if ((!return_self)){

	  CpGIDs.neighborhood <- data[c(neighbors), "CpGID"]	  
	  CpGID.current <- data[c(idx), "CpGID"]

	  starts.neighborhood <- as.integer(unlist(lapply(X = CpGIDs.neighborhood, FUN = function(x){strsplit(x = x, split = "_")[[1]][2]})))
	  starts.current <- as.integer(unlist(lapply(X = CpGID.current, FUN = function(x){strsplit(x = x, split = "_")[[1]][2]})))
	  abs_distances <- abs(starts.current - starts.neighborhood)
	  to.remove_idxs <- which(abs_distances > max_distance)

	  if (length(to.remove_idxs) > 0){
		neighbors <- neighbors[-to.remove_idxs]
	  	}

	  }

    }
	
    return(sort(neighbors)) 

}