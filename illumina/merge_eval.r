

# get output from merge_and_map.r
load("counts.rda")


# Look at correlations between all pairs of samples
#   The point of this plot is that the highest correlations should be found between related samples of the same replicates or same selections, not between 
#   unrelated samples (black circles in plot).

for (lib in names(bio_rep)) {
    print(sprintf("Correlations between all samples in lib %s",lib))
    cor_mat = cor(counts[[lib]][,names(samples[[1]])], method="pearson", use="complete.obs")
    print(cor_mat)


    # look at correlations between replicates of the same samples vs all pairs
    rep_mat = matrix(0, nrow=nrow(cor_mat), ncol=ncol(cor_mat), dimnames=dimnames(cor_mat))
    sample_mat = matrix(0, nrow=nrow(cor_mat), ncol=ncol(cor_mat), dimnames=dimnames(cor_mat))
    sample_idx = 0
    replicate_idx = 0
    replicate_idx_map = list()
    for (sample in names(bio_rep[[lib]])) {
        sample_idx = sample_idx +1
        for (br1 in bio_rep[[lib]][[sample]]) {
	    # mark all replicates of same sample
            sel_rep1 = sprintf("%s-%s", sample, br1)
            for (br2 in bio_rep[[lib]][[sample]]) {
                sel_rep2 = sprintf("%s-%s", sample, br2)
                rep_mat[sel_rep1,sel_rep2] = sample_idx
            }
	    # maintain a map from replicate string to an integer index
	    if (! br1 %in% names(replicate_idx_map)) {	    
	        replicate_idx = replicate_idx +1
		replicate_idx_map[[br1]] = replicate_idx
	    }
	    # mark all samples of same replica
	    for (sample2 in names(bio_rep[[lib]])) {
	        sel_rep2 = sprintf("%s-%s", sample2, br1)
		sample_mat[sel_rep1,sel_rep2] = replicate_idx_map[[br1]]
	    }
        }
    }
    cors = cor_mat[upper.tri(cor_mat)]
    sample_idx = rep_mat[upper.tri(rep_mat)]
    rep_idx = sample_mat[upper.tri(sample_mat)]

    quartz(width=8, height=8)
    
    s = sprintf("Library %s sample correlations", lib)
    plot(seq_along(cors), cors, col=sample_idx+1, pch=1+rep_idx, ylim=c(0.2,1.1), xlab="Sample pairs", ylab="Pearson correlation", main=s, lwd=2)
    n_col = max(sample_idx)+1
    legend("top", c("Non-replicate pair",sprintf("Sample '%s' rep. pair",names(bio_rep[[lib]]))), pch=16, col=seq(n_col), ncol=min(3,ceiling(n_col/2.0)))
    n_pch = max(rep_idx)+1
    legend("bottom", c("Different rep. pair",sprintf("Rep. %s sample",names(replicate_idx_map))), pch=seq(n_pch), col=1, lwd=2, lty=NA, ncol=min(3,ceiling(n_pch/2.0)))
    quartz.save(sprintf("all_pair_cor_%s.png",lib), type="png")
    
}
