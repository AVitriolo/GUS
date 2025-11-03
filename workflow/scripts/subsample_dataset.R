source("workflow/scripts/helpers/compute_corr.R")

options(scipen=999)                                                                                       # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                                        # read args

input_path                               <- args$input_path                                               # take input path
output_path 	                         <- args$output_path                                              # take output path

CpGs <- paste0("CpG_",1:1000)                                                                             # simulate CpGs

vec_sampled_avg_corr <- c()                                                                               # initialize avg correlations vector
vec_sampled_dfs <- c()                                                                                    # initialize sampled_dfs vector
abs_corr_threshold <- 0.5                                                                                 # set avg correlation threshold

times_sampled <- rep(x = 0, length(CpGs)); names(times_sampled) <- CpGs                                   # initialize number of CpGs sampled
expected_sampling_prob <- rep((1 / length(CpGs)),length(CpGs))                                            # expected in uniform distribution
B <- 1000                                                                                                 # montecarlo iterations
significance_thr <- 0.05                                                                                  # significance threshold
pval <- 0                                                                                                 # initialize p-value

pseudocount <- 0.0001                                                                                     # pseudocount for weight calculation
weight <- 1 / (times_sampled + pseudocount)                                                               # compute weight
random_seed <- 123                                                                                        # set random seed

lower_bound_CpG <- 25                                                                                     # min times a CpG has to be sampled at all
min_num_sampled <- 25                                                                                     # min number of times at least min_sampled_CpGs have to be sampled
min_sampled_CpGs <- length(CpGs) / 2                                                                      # min number of CpGs to be sampled min_num_sampled times

cond <- !((pval >= significance_thr) & 
            (sum((times_sampled >= min_num_sampled)) >= min_sampled_CpGs) &
                (min(times_sampled) >= lower_bound_CpG)
                )                                                                                         # condition to break the loop: distribution must be uniform with 
                                                                                                          # at least min_sampled_CpGs to be sampled min_num_sampled times
                                                                                                          # and the least sampled CpG must be sampled at least 

set.seed(random_seed)

while(cond)                                                                                               # while
{

sampled_CpGs <- sample(x = CpGs, size = 100, replace = F, prob = weight)                                  # random sample, without replacement

df_id <- paste0(sort(sampled_CpGs), collapse = "_")

sampled_df <- original_df[,sampled_CpGs]
vec_sampled_dfs[df_id] <- sampled_df

sampled_avg_corr <- compute_corr(input_df = sampled_df)
vec_sampled_avg_corr[df_id] <- sampled_avg_corr

times_sampled[sampled_CpGs] <- times_sampled[sampled_CpGs] + 1                                            # update times_sampled
weight <- 1 / (times_sampled + pseudocount)                                                               # update weights
gof_test <- chisq.test(x = times_sampled, p = expected_sampling_prob, simulate.p.value = TRUE, B = B)     # compute goodness of fit test
pval <- gof_test$p.value                                                                                  # take p.value

cond <- !((pval >= significance_thr) &                                                                    # update condition
            (sum((times_sampled >= min_num_sampled)) >= min_sampled_CpGs) &
                (min(times_sampled) >= lower_bound_CpG)
                )           
                                                                              
}

selected_dfs_mask <- names(vec_sampled_avg_corr[vec_sampled_avg_corr < abs_corr_threshold])
selected_dfs <- vec_sampled_dfs[selected_dfs_mask]

for(idx in seq(selected_dfs)){

    output_df <- selected_dfs[idx]
    output_file <- paste0(output_path, "/", selected_dfs_mask[idx])

    write.table(output_df,
            file=output_file,
            row.names=T,
            col.names=T,
            quote=F,
            sep="\t")
            
}