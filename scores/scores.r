options(width=160)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    counts_file = "../call/counts_prkn_bc.csv"
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript scores.r  <counts.csv>")
    quit(save="no")
} else {
    counts_file = args[1]
}
print(sprintf("Input counts file: %s",counts_file))

# Settings
settings = list(file=counts_file)

# Samples
settings$bio_reps = list()
settings$bio_reps[["ao"]] = sprintf("%d", seq(3))
settings$bio_reps[["dmso"]] = sprintf("%d", seq(3))
settings$n_facs_bins = 4
settings$idx_facs_bins = seq(settings$n_facs_bins)
settings$pop_facs_bins = list()
settings$pop_facs_bins[["ao"]] = c(0.15, 0.30, 0.35, 0.20)
settings$pop_facs_bins[["dmso"]] = c(0.15, 0.30, 0.35, 0.20)

# Required number of reads in all FACS bins needed to trust PSI of a variant
settings$threshold_counts_per_rep = 50

# Require this number of replicate measurements to keep score
settings$min_rep = 2

stopifnot( length(settings$idx_facs_bins) == settings$n_facs_bins )
stopifnot(all( sapply(settings$pop_facs_bins, length) == settings$n_facs_bins ))

# names of coulmns with counts
replicas = c()
cn_counts = c()
for (cell in names(settings$bio_reps)) {
    for (bio_rep in settings$bio_reps[[cell]]) {
        replicas = c(replicas, sprintf("%s_%s",cell,bio_rep))
        for (facs_bin in seq(settings$n_facs_bins)) {
            cn_counts =c(cn_counts, sprintf("%s_bin%d.%s",cell,facs_bin,bio_rep))
        }
    }
}

# needed for making output files
wt = list()
wt[["name"]] = "PRKN"
# wt[["dna"]] = "ATGATCGTGTTTGTGCGCTTCAACAGCTCCCACGGCTTTCCAGTGGAGGTGGACTCTGATACCAGCATCTTCCAGCTGAAGGAGGTGGTGGCAAAGAGGCAGGGAGTGCCAGCCGACCAGCTGAGAGTGATCTTCGCCGGCAAGGAGCTGAGAAACGATTGGACAGTGCAGAATTGCGACCTGGATCAGCAGTCCATCGTGCACATCGTGCAGCGGCCCTGGAGAAAGGGACAGGAGATGAACGCAACCGGCGGCGACGATCCAAGGAATGCAGCAGGAGGATGTGAGAGGGAGCCTCAGTCTCTGACCAGAGTGGACCTGTCTAGCTCCGTGCTGCCTGGCGATAGCGTGGGCCTGGCCGTGATCCTGCACACCGACTCCAGGAAGGATTCTCCACCTGCAGGCAGCCCAGCAGGCCGGAGCATCTATAACTCCTTTTACGTGTATTGCAAGGGACCATGTCAGCGGGTGCAGCCTGGCAAGCTGAGAGTGCAGTGCAGCACCTGTAGGCAGGCCACACTGACCCTGACACAGGGCCCTTCCTGCTGGGACGATGTGCTGATCCCAAATAGAATGTCTGGCGAGTGCCAGAGCCCTCACTGTCCAGGCACATCCGCCGAGTTCTTTTTCAAGTGTGGCGCCCACCCCACCTCCGACAAGGAGACATCTGTGGCCCTGCACCTGATCGCCACCAACTCTCGGAATATCACATGCATCACCTGTACAGACGTGAGATCTCCTGTGCTGGTGTTTCAGTGCAACAGCAGGCACGTGATCTGCCTGGATTGTTTCCACCTGTATTGCGTGACCAGGCTGAATGACCGCCAGTTTGTGCACGATCCACAGCTGGGATACTCCCTGCCATGCGTGGCAGGATGTCCCAACTCTCTGATCAAGGAGCTGCACCACTTCAGGATCCTGGGCGAGGAGCAGTACAATCGCTATCAGCAGTACGGAGCAGAGGAGTGCGTGCTGCAGATGGGAGGCGTGCTGTGCCCACGCCCCGGCTGTGGAGCAGGCCTGCTGCCTGAGCCAGACCAGAGGAAGGTGACCTGCGAGGGAGGAAACGGCCTGGGATGTGGATTTGCCTTCTGCCGGGAGTGTAAGGAGGCCTATCACGAGGGCGAGTGCAGCGCCGTGTTTGAGGCATCCGGAACCACAACCCAGGCATACAGGGTGGATGAGAGAGCAGCAGAGCAGGCAAGATGGGAGGCAGCATCCAAGGAGACAATCAAGAAGACAACCAAGCCCTGCCCTAGGTGTCACGTGCCTGTGGAGAAGAATGGCGGCTGCATGCACATGAAGTGTCCACAGCCCCAGTGCAGACTGGAGTGGTGTTGGAACTGCGGCTGTGAGTGGAATAGGGTGTGCATGGGCGACCACTGGTTCGATGTG"
wt[["aa"]] = "MIVFVRFNSSHGFPVEVDSDTSIFQLKEVVAKRQGVPADQLRVIFAGKELRNDWTVQNCDLDQQSIVHIVQRPWRKGQEMNATGGDDPRNAAGGCEREPQSLTRVDLSSSVLPGDSVGLAVILHTDSRKDSPPAGSPAGRSIYNSFYVYCKGPCQRVQPGKLRVQCSTCRQATLTLTQGPSCWDDVLIPNRMSGECQSPHCPGTSAEFFFKCGAHPTSDKETSVALHLIATNSRNITCITCTDVRSPVLVFQCNSRHVICLDCFHLYCVTRLNDRQFVHDPQLGYSLPCVAGCPNSLIKELHHFRILGEEQYNRYQQYGAEECVLQMGGVLCPRPGCGAGLLPEPDQRKVTCEGGNGLGCGFAFCRECKEAYHEGECSAVFEASGTTTQAYRVDERAAEQARWEAASKETIKKTTKPCPRCHVPVEKNGGCMHMKCPQPQCRLEWCWNCGCEWNRVCMGDHWFDV"


##
##  Build a data frame of coutns per amino acid variant
##
# consider keeping syn. WT here but still merging different DNA variants of the same aa_variants
counts = read.csv(counts_file)
stopifnot(all( c("subst_aa","subst_dna",cn_counts) %in% colnames(counts) ))

# consider single amino acid variants and wild-types (DNA level WT, i.e. excluding synonymous WT)
# ic_use = which(counts$n_subst_aa==1 | counts$n_subst_dna==0)

# First, find all single amino acid variants
ic_aa1 = which(counts$n_subst_aa==1)
print(sprintf("Aggregating read counts for %d of %d (%.1f%%) barcodes of single amino acid variants",
              length(ic_aa1), nrow(counts), length(ic_aa1)/nrow(counts)*100))
for (cn in cn_counts) {
    print(sprintf("   using %8d of %8d (%.1f%%) read counts for %s", sum(counts[ic_aa1,cn]), sum(counts[,cn]), sum(counts[ic_aa1,cn])/sum(counts[,cn])*100, cn))
}

# aggregate barcode counts per amino acid substitution
df = counts[ic_aa1,c("subst_aa",cn_counts)]
agg = aggregate(df[,cn_counts], by=list(df$subst_aa), sum)
aa1 = data.frame(subst=agg[,1], wt=NA, resi=NA, mut=NA)
aa1[,cn_counts] = agg[,cn_counts]

# extract amino acids and residue number from substitution string
nc = nchar(aa1$subst)
aa1$wt = substr(aa1$subst,1,1)
aa1$resi = as.numeric( substr(aa1$subst,2,nc-1) )
aa1$mut = substr(aa1$subst,nc,nc)

# re-oreder variants
aa1 = aa1[order(aa1$resi,aa1$mut),]
rownames(aa1) = NULL

# Second, synonymous wild-type variants with a single DNA substitutions
ic_syn1 = which(counts$n_subst_aa == 0 & counts$n_subst_dna == 1)
print(sprintf("Aggregating read counts for %d of %d (%.1f%%) barcodes of single nucl synonymous wild types",
              length(ic_syn1), nrow(counts), length(ic_syn1)/nrow(counts)*100))
for (cn in cn_counts) {
    print(sprintf("   using %8d of %8d (%.1f%%) read counts for %s", sum(counts[ic_syn1,cn]), sum(counts[,cn]), sum(counts[ic_syn1,cn])/sum(counts[,cn])*100, cn))
}

# aggregate barcode counts per nucleotide substitution
df = counts[ic_syn1,c("subst_dna",cn_counts)]

# calculate amino acid numbering and make substitution string
df$resi = (as.numeric( substr(df$subst_dna,2,nchar(df$subst_dna)-1) )-1) %/% 3 + 1
df$wt = strsplit(wt$aa,"")[[1]][df$resi]
df$subst_aa = paste0(df$wt, df$resi, "=")

agg = aggregate(df[,cn_counts], by=list(df$subst_aa), sum)
syn1 = data.frame(subst=agg[,1], wt=NA, resi=NA, mut=NA)
syn1[,cn_counts] = agg[,cn_counts]

# extract amino acids and residue number from substitution string
nc = nchar(syn1$subst)
syn1$wt = substr(syn1$subst,1,1)
syn1$resi = as.numeric( substr(syn1$subst,2,nc-1) )
syn1$mut = substr(syn1$subst,nc,nc)

# re-oreder variants
syn1 = syn1[order(syn1$resi),]
rownames(syn1) = NULL

# Third, merge barcodes for wildtype
ic_wt = which(counts$n_subst_aa == 0 & counts$n_subst_dna == 0)
print(sprintf("Aggregating read counts for %d of %d (%.1f%%) barcodes of wild type (DNA level)",
              length(ic_wt), nrow(counts), length(ic_wt)/nrow(counts)*100))
for (cn in cn_counts) {
    print(sprintf("   using %8d of %8d (%.1f%%) read counts for %s", sum(counts[ic_wt,cn]), sum(counts[,cn]), sum(counts[ic_wt,cn])/sum(counts[,cn])*100, cn))
}
wtrow = data.frame(subst="WT", wt="", resi=0, mut="")
wtrow[,cn_counts] = apply(counts[ic_wt,cn_counts], MARGIN=2, sum)

# Merge into a single data frame
raw = rbind(wtrow, aa1, syn1)

# All substitutions should be unique
stopifnot( nrow(raw) == length(unique(raw$subst)) )

# Each barcode should only be used once
ic = c(ic_wt, ic_aa1, ic_syn1)
stopifnot( length(ic) == length(unique(ic)) )

print(sprintf("Collected data for %d substitutions using %d barcodes (%.1f%%)", nrow(raw), length(ic), length(ic)/nrow(counts)*100))


##
##  Calculate average-bin scores
##

# Calculate RPM with psudo counts 
settings$pseudocounts = 1
calc_rpm = function(v) { v/sum(v)*10^6 }
raw_rpm = data.frame(subst=raw$subst)
for (cn in cn_counts) {
    raw_rpm[,cn] = calc_rpm(raw[,cn] + settings$pseudocounts)
}

# Function to calculate PSI based on rows in a data frame
protein_stability_index = function(df, name, indices=NA, populations=NA) {
    # Each FACS gate should have an index
    if (all(is.na(indices))) {
        indices = seq(ncol(df))
    } else {
        stopifnot( length(indices) == ncol(df) )
    }

    calc_psi = function(rpm,idx) {
        rpm_sum = sum(rpm)
	psi = sum( rpm/rpm_sum * idx )
	return( c(rpm_sum, psi) )
    }
    calc_psi_pop = function(rpm,idx,pop) {
        rpm_sum = sum(rpm * pop)
	psi = sum( rpm/rpm_sum * idx * pop )
	return( c(rpm_sum, psi) )
    }

    # Each FACS gate should have a population
    if (all(is.na(populations))) {
        # populations = rep(1/ncol(df), times=ncol(df))
	print(sprintf("    Calc. PSI for %s using indices %s",name,paste0(indices,collapse=",")))
	psi = t(apply(df, MARGIN=1, calc_psi, idx=indices))
    } else {
        stopifnot( length(populations) == ncol(df) )
	print(sprintf("    Calc. PSI for %s using indices %s and populations %s",name,paste0(indices,collapse=","),paste0(populations,collapse=",")))
	psi = t(apply(df, MARGIN=1, calc_psi_pop, idx=indices, pop=populations))
    }
    
    # put in a data frame and give column name
    ret_df = data.frame(psi)
    colnames(ret_df) = paste0(name,c("_rpm_sum","_psi"))
    return(ret_df)
}


# Calculate PSI
raw_psi = data.frame(subst=raw$subst)
for (cell in names(settings$bio_reps)) {
    for (bio_rep in settings$bio_reps[[cell]]) {
        replica = sprintf("%s_%s",cell,bio_rep)
        cn_gates = sprintf("%s_bin%d.%s",cell,seq(settings$n_facs_bins),bio_rep)
        cn_reads = sprintf("%s_read_sum",replica)
        print(sprintf("Calculate PSI for %s",replica))
    
        # Total number of reads per variant over all gates
        raw_psi[,cn_reads] = apply(raw[,cn_gates], MARGIN=1, sum)

        # PSI per variant
        raw_psi = cbind(raw_psi, protein_stability_index(raw_rpm[,cn_gates],  name = replica,
	                                                 indices = settings$idx_facs_bins,
							 populations = settings$pop_facs_bins[[cell]]))

        i_rm = which(raw_psi[,cn_reads] < settings$threshold_counts_per_rep)
	if (length(i_rm) > 0) {
            print(sprintf("    Removing %d PSI scores with less than %d total reads",length(i_rm),settings$threshold_counts_per_rep))
	    print(raw_psi[i_rm,sprintf("%s_%s",replica,c("read_sum","psi"))])
            raw_psi[i_rm,sprintf("%s_psi",replica)] = NA
	}
    }
}

# plot(raw_psi$bio1_facs1_psi, raw_psi$bio1_facs2_psi)
# i = which(raw_psi$bio1_facs1_read_sum < 100 | raw_psi$bio1_facs2_read_sum < 100)
# points(raw_psi[i,"bio1_facs1_psi"], raw_psi[i,"bio1_facs2_psi"], pch=20, col=2)


# Calculate a mean of replica
for (cell in names(settings$bio_reps)) {
    replica_cns = sprintf("%s_%s_psi",cell,settings$bio_reps[[cell]])
    raw_psi[,paste0(cell,"_mean")] = apply(raw_psi[,replica_cns], MARGIN=1, mean, na.rm=T)
    raw_psi[,paste0(cell,"_sd")]   = apply(raw_psi[,replica_cns], MARGIN=1, sd, na.rm=T)

    # calc mean of cell representation
    rpm_cns = sprintf("%s_%s_rpm_sum",cell,settings$bio_reps[[cell]])
    raw_psi[,paste0(cell,"_cpm")] = apply(raw_psi[,rpm_cns], MARGIN=1, mean, na.rm=T)
}

i_rm_all = c()

# Remove scores with less than requested replicates measured
for (treat in names(settings$bio_reps)) {
    cn = sprintf("%s_rep",treat)
    cn_reps = sprintf("%s_%s_psi",treat,settings$bio_reps[[treat]])
    raw_psi[,cn] = apply(raw_psi[,cn_reps], MARGIN=1, function(v){ sum(! is.na(v)) })
    print(sprintf("Replicate distribution for %s scores",treat))
    print(table(raw_psi[,cn]))
    i_rm = which(raw_psi[,cn] < settings$min_rep)
    print(sprintf("Removing %d mean PSI scores with less than %d measurements", length(i_rm), settings$min_rep))
    raw_psi[i_rm,sprintf("%s_mean",treat)] = NA
    raw_psi[i_rm,sprintf("%s_sd",treat)] = NA
    i_rm_all = c(i_rm_all,i_rm)
}

i_rm_all = unique(i_rm_all)
i_rm_all = i_rm_all[order(i_rm_all)]
print(sprintf("All %d removed measurements", length(i_rm_all)))
cns = c("subst","ao_1_read_sum","ao_1_psi","ao_2_read_sum","ao_2_psi","ao_3_read_sum","ao_3_psi",
        "dmso_1_read_sum","dmso_1_psi","dmso_2_read_sum","dmso_2_psi","dmso_3_read_sum","dmso_3_psi",
	"ao_mean","ao_sd","dmso_mean","dmso_sd","ao_rep","dmso_rep")
print(raw_psi[i_rm_all,cns])


##
## Normalize to FACS distributions
##
# per library, plot replica PSI distributions, weighted distributions, per-bin average SD and the FACS transformed distributions
# Map each sorting replica to a FACS distribution that represents the library
load("facs.rda")      # get facs_events, facs_table, and facs_set

# select a FACS data set per treatment
facs_table$rep = NA
facs_table[which(facs_table$name == "P+K 22"),  "rep"] = "ao_mean"
facs_table[which(facs_table$name == "P+K D"),  "rep"]  = "dmso_mean"

# Plotting
quartz(width=8, height=6)
rand_unif = runif(10^6, min=0.0, max=1.0)
par(ask=T)

for (treat in c("ao","dmso")) {
    cn_ag = sprintf("%s_mean", treat)
    fi = which(facs_table$rep==cn_ag)
    stopifnot(length(fi)==1)

    xf = log10( facs_events[[facs_table[fi,"name"]]][,facs_set$cn_ratio] )
    pm = quantile(xf, c(0.0001,0.9999))
    xf = xf[which(pm[1] < xf & xf < pm[2])]
    hf = density(xf, bw="nrd", n=1024)

    # Function that transform uniformly distributed numbers in [0,1] to FACS distribution
    facs_quant = splinefun(cumsum(hf$y)/sum(hf$y), hf$x, method="hyman")

    # Score distribution
    xs = raw_psi[,cn_ag]
    hs = density(xs, bw="nrd", n=1024, na.rm=T)
	
    # Function that transform degron-score-distributed numbers to uniform distribution in [0,1]
    score_cdf = splinefun(hs$x, cumsum(hs$y)/sum(hs$y), method="hyman")

    # Function that transforms score distributed numbers to FACS distributed numbers
    # this produces negative numbers, can I avoid that?
    score2facs = function(x){ facs_quant(score_cdf(x)) }
    score2facs_deriv = function(x){ facs_quant(score_cdf(x), deriv=1) * score_cdf(x, deriv=1) }

    cn_facs = sprintf("%s_fnorm", treat)
    raw_psi[,cn_facs] = score2facs(raw_psi[,cn_ag])
    cn_facs_sd = sprintf("%s_fnorm_sd", treat)
    raw_psi[,cn_facs_sd] = raw_psi[,sprintf("%s_sd", treat)] * score2facs_deriv(raw_psi[,cn_ag])
    
    # Repeat with a score histogram that is weighted by the representation among cells
    cn_cells = sprintf("%s_cpm", treat)
    is = which(! is.na(xs))
    hw = density(xs[is], bw="ucv", weights=raw_psi[is,cn_cells]/sum(raw_psi[is,cn_cells]), n=1024)
    score_w_cdf = splinefun(hw$x, cumsum(hw$y)/sum(hw$y))
    scorew2facs = function(x){ facs_quant(score_w_cdf(x)) }
    cn_wfacs = sprintf("%s_fwnorm", treat)
    raw_psi[,cn_wfacs] = scorew2facs(raw_psi[,cn_ag])

    # Plot
    # cn_rep = sub("_psi","",cn)
    hff = hist(xf, breaks=200, plot=F)
    hr = hist(facs_quant(rand_unif), breaks=200, plot=F)
    hsf = hist(raw_psi[,cn_facs], breaks=200, plot=F)
    hsfw = hist(raw_psi[,cn_wfacs], breaks=200, plot=F)
    # hsfw = whist(raw_psi[,cn_wfacs], raw_psi[,cn_cells], breaks=hsf$breaks, plot=F)
    # log10 xlim=c(-1.5,0.0) non-log xlim=c(0.0,0.6)
    plot(0,0,col=0, xlim=pm, ylim=c(0.0,10), xlab="Score", ylab="Density",
         main=sprintf("%s, %d scores, %d breaks",cn_ag,sum(! is.na(raw_psi[,cn_ag])),settings$fnorm_hs_nbreaks))
    lines(hsf$mids, hsf$density, col=1, lwd=2)
    lines(hr$mids, hr$density, col=3)
    lines(hff$mids, hff$density, col=2)
    lines(hsfw$mids, hsfw$density,col=4)
    lines(hf$x, hf$y, col=5)
    hsf_lab = sprintf("Transformed score distribution (%.2f - %.2f)", hsf$breaks[1], hsf$breaks[length(hsf$breaks)])
    hsfw_lab = sprintf("Transf. score distribution weighted (%.2f - %.2f)", hsfw$breaks[1], hsfw$breaks[length(hsfw$breaks)])
    legend("top", c(hsf_lab,hsfw_lab,"Target FACS distribution","Transformed unif. random numbers","Kerneled target"),
           col=c(1,4,2,3,5), lty=1)
    # legend("top", c(hsf_lab,"Target FACS distribution","Transformed unif. random numbers"),
    #        col=c(1,2,3), lty=1, lwd=c(2,1,1))
}


##
## Derived scores and normalization
##

iwt = which(raw_psi$subst == "WT")
iC431A = which(raw_psi$subst == "C431A")
stopifnot(all(raw$subst == raw_psi$subst))
inons = which(raw$mut == "*")
isyn = which(raw$mut == "=")

# Scale such that mean nonsence score equals C431A score
dmso_factor = (raw_psi[iC431A,"ao_mean"] - mean(raw_psi[inons,"ao_mean"], na.rm=T)) / (raw_psi[iC431A,"dmso_mean"] - mean(raw_psi[inons,"dmso_mean"], na.rm=T))

# Calculate activity score
activity_ag = raw_psi$ao_mean - dmso_factor*raw_psi$dmso_mean
activity_ag_sd = sqrt(raw_psi$ao_sd^2 + (dmso_factor*raw_psi$dmso_sd)^2)

# Min-max normalization
scale = 1/(activity_ag[iwt] - activity_ag[iC431A])
raw_psi$activity = (activity_ag - activity_ag[iC431A]) * scale
raw_psi$activity_sd = activity_ag_sd * scale
print(sprintf("Calc. activity score using a DMSO factor of %.4f and min-max normalized to WT (%.4f) and C431A (%.4f) with a scaling of %.4f",
              dmso_factor, activity_ag[iwt], activity_ag[iC431A], scale))

print(sprintf("Activity score uncertainty from %d WT syn %.3f and %d nonsense %.3f",
      sum(! is.na(raw_psi[isyn,"activity"])), sd(raw_psi[isyn,"activity"], na.rm=T),
      sum(! is.na(raw_psi[inons,"activity"])), sd(raw_psi[inons,"activity"], na.rm=T)))
n = sum(! is.na(raw_psi$activity))
# stopifnot(all( is.na(raw_psi$activity) == is.na(raw_psi$activity_sd) ))
print(sprintf("Summary of %d (%.1f%% incl. *=) score std. dev. (signal=1)", n, n*100.0/(nchar(wt[["aa"]])*21) ))
print(summary(raw_psi$activity_sd))

# Keima difference score
activity_fn = raw_psi$ao_fnorm - raw_psi$dmso_fnorm
activity_fn_sd = sqrt(raw_psi$ao_fnorm_sd^2 + raw_psi$dmso_fnorm_sd^2)

# Min-max normalization of FACS-normalized scores
scale_fn = 1/(activity_fn[iwt] - activity_fn[iC431A])
raw_psi$activity_fnorm = (activity_fn - activity_fn[iC431A]) * scale_fn
raw_psi$activity_fnorm_sd = activity_fn_sd * scale_fn
print(sprintf("Calc. fluorescence-normalized activity score min-max normalized to WT (%.4f) and C431A (%.4f) with a scaling of %.4f",
              activity_fn[iwt], activity_fn[iC431A], scale_fn))

# # Two-point normalization of FACS-normalized scores
# #   dy/dx = (y1-y0)/(x1-x0);  y = (x-x0) * dy/dx + y0; with y the norm. values, x the original values and points (x0,y0) andd (x1,y1)
# ltp_wt = log10(1.0)
# ltp_c431a = log10(0.296)
# scale_fn = (ltp_wt - ltp_c431a) / (activity_fn[iwt] - activity_fn[iC431A])
# raw_psi$activity_fnorm = (activity_fn - activity_fn[iC431A]) * scale_fn + ltp_c431a
# raw_psi$activity_fnorm_sd = activity_fn_sd * scale_fn
# print(sprintf("Calc. fluorescence-normalized activity score two-point normalized with WT (%.4f) value 1.0 and C431A (%.4f) LTP value %.4f with a scaling of %.4f",
#               activity_fn[iwt], activity_fn[iC431A], ltp_c431a, scale_fn))


print(sprintf("FACS-normalized score uncertainty from %d WT syn %.3f and %d nonsense %.3f",
      sum(! is.na(raw_psi[isyn,"activity_fnorm"])), sd(raw_psi[isyn,"activity_fnorm"], na.rm=T),
      sum(! is.na(raw_psi[inons,"activity_fnorm"])), sd(raw_psi[inons,"activity_fnorm"], na.rm=T)))
n = sum(! is.na(raw_psi$activity_fnorm))
# stopifnot(all( is.na(raw_psi$activity_fnorm) == is.na(raw_psi$activity_fnorm_sd) ))
print(sprintf("Summary of %d (%.1f%% incl. *=) FACS-normalized score std. dev. (signal=1)", n, n*100.0/(nchar(wt[["aa"]])*21) ))
print(summary(raw_psi$activity_fnorm_sd))


##
##  Per residue data
##

# coverage, which single amino acid subst have calculated scores
aa_one = strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

# only count mutations that have both conditions measured
i_aa1 = which(raw$mut %in% aa_one & ! is.na(raw_psi$activity))
agg = aggregate(raw[i_aa1,"mut"], by=list(raw[i_aa1,"resi"]), paste0, collapse="")

residue = data.frame(resi=seq(nchar(wt[["aa"]])), wt=strsplit(wt[["aa"]],"")[[1]])
residue[,"mut"] = agg[match(residue$resi,agg[,1]),2]
residue[which(is.na(residue$mut)),"mut"] = ""
residue$coverage = sapply(residue$mut, nchar)

# Median scores
agg = aggregate(raw_psi[i_aa1,"activity"], by=list(raw[i_aa1,"resi"]), median, na.rm=T)
residue$med_score = agg[match(residue$resi,agg[,1]),2]

# agg = aggregate(raw_psi[i_aa1,"dmso_mean"], by=list(raw[i_aa1,"resi"]), median, na.rm=T)
# residue$med_dmso = agg[match(residue$resi,agg[,1]),2]

# agg = aggregate(raw_psi[i_aa1,"ao_mean"], by=list(raw[i_aa1,"resi"]), median, na.rm=T)
# residue$med_ao = agg[match(residue$resi,agg[,1]),2]

# Nonsense scores
in_ctrl = which(raw$mut == "*" & ! is.na(raw_psi$activity))
in_ctrl_match = in_ctrl[match(residue$resi,raw[in_ctrl,"resi"])]
residue$nons_score = raw_psi[in_ctrl_match,"activity"]

# in_ctrl = which(raw$mut == "*" & ! is.na(raw_psi$dmso_norm))
# in_ctrl_match = in_ctrl[match(residue$resi,raw[in_ctrl,"resi"])]
# residue$nons_dmso = raw_psi[in_ctrl_match,"dmso_norm"]

# in_treat = which(raw$mut == "*" & ! is.na(raw_psi$ao_norm))
# in_treat_match = in_treat[match(residue$resi,raw[in_treat,"resi"])]
# residue$nons_ao = raw_psi[in_treat_match,"ao_norm"]

# Synonymous WT scores
is_ctrl = which(raw$mut == "=" & ! is.na(raw_psi$activity))
is_ctrl_match = is_ctrl[match(residue$resi,raw[is_ctrl,"resi"])]
residue$syn_score = raw_psi[is_ctrl_match,"activity"]

# is_ctrl = which(raw$mut == "=" & ! is.na(raw_psi$dmso_norm))
# is_ctrl_match = is_ctrl[match(residue$resi,raw[is_ctrl,"resi"])]
# residue$syn_dmso = raw_psi[is_ctrl_match,"dmso_norm"]

# is_treat = which(raw$mut == "=" & ! is.na(raw_psi$ao_norm))
# is_treat_match = is_treat[match(residue$resi,raw[is_treat,"resi"])]
# residue$syn_ao = raw_psi[is_treat_match,"ao_norm"]

# report coverage
print(sprintf("Data has %d activity measurements,  %d of %d (%.1f%%) single variants",
              sum(! is.na(raw_psi$delta_psi)),
              length(i_aa1), 19*nchar(wt[["aa"]]), 100*length(i_aa1)/(19*nchar(wt[["aa"]])) ))
# print(sprintf("  For ctrl cells,  %d (%.1f%%) nonsense, %d (%.1f%%) synonymous and WT",
# 	      length(in_ctrl), 100*length(in_ctrl)/nchar(wt[["aa"]]),
# 	      length(is_ctrl), 100*length(is_ctrl)/nchar(wt[["aa"]]) ))
# print(sprintf("  For treated cells,  %d (%.1f%%) nonsense, %d (%.1f%%) synonymous and WT",
# 	      length(in_treat), 100*length(in_treat)/nchar(wt[["aa"]]),
# 	      length(is_treat), 100*length(is_treat)/nchar(wt[["aa"]]) ))


#
# Dump 
#

# full data set
write.csv(raw_psi, file="prkn_activity_full.csv", row.names=F, quote=F)

# Dump CSV with per-residues data
write.csv(residue, file="prkn_activity_residues.csv", row.names=F, quote=F)

# new data frame without missing scores
i = which(! is.na(raw_psi$activity))
prkn_func = raw_psi[i,c("subst","ao_mean","ao_sd","dmso_mean","dmso_sd","activity","activity_sd","activity_fnorm","activity_fnorm_sd")]
colnames(prkn_func) = c("var","ao","ao_sd","dmso","dmso_sd","activity","activity_sd","act_keima","act_keima_sd")
nc = nchar(prkn_func$var)
pf_resi = as.numeric(substr(prkn_func$var,2,nc-1))
pf_resi[which(prkn_func$var=="WT")] = 0
pf_mut = substr(prkn_func$var, nc, nc)
prkn_func = prkn_func[order(pf_resi,pf_mut),]
print(sprintf("Final score set has %d variants, incl. %d WT, %d nonsense and %d synonymous (of %d residues)",
              length(i), sum(prkn_func$var=="WT"), sum(raw[i,"mut"] == "*"), sum(raw[i,"mut"] == "="), nrow(residue)))

# Dump CSV with substitutions
write.csv(prkn_func, file="prkn_activity.csv", row.names=F, quote=F)

# save everything in R data file
save(raw, raw_psi, residue, prkn_func, wt, settings, scale, cn_counts, replicas, file="prkn_activity.rda")

# dump prism file
f = file("prism_mave_prkn_func.txt", "wt")
write("# --------------------", f)
write("# version: 1", f)
write("# protein:", f)
write("#     name: PRKN", f)
write(sprintf("#     sequence: %s",wt[["aa"]]), f)
write("#     organism: Homo sapiens (Human)", f)
write("#     uniprot: O60260", f)
write("# mave:", f)
write("#     organism: Homo sapiens (Human)", f)
write("#     cloning: Landing pad", f)
write("#     expression: Overexpression", f)
write("#     technology: Sort-seq", f)
write("#     doi: unpublished", f)
write("#     year: 2025", f)
write("# variants:", f)
write(sprintf("#     number: %d",nrow(prkn_func)), f)
write(sprintf("#     coverage: %.2f",nrow(prkn_func)/(20*nchar(wt[["aa"]])+1) ), f) # including nons. and wt
write(sprintf("#     depth: %.2f",mean(nchar(residue[which(residue$coverage > 0),"mut"]))), f)
write("#     width: single mutants", f)
write("# columns:", f)
write("#     dmso_psi: PSI score of non-treated cells", f)
write("#     dmso_sd: Score standard deviation between replicates", f)
write("#     ao_psi: PSI score of AO treated cells", f)
write("#     ao_sd: Score standard deviation between replicates", f)
write("#     score: Functional score is ao_psi-dmso_psi normalized to WT (score one) and C431A (score zero)", f)
write("#     score_sd: Error of functional score is propagated from ao_sd and dmso_sd", f)
write("#     keima: Functional score is ao_psi-dmso_psi normalized to WT (score one) and C431A (score zero)", f)
write("#     keima_sd: Error of functional score is propagated from ao_sd and dmso_sd", f)
write("# --------------------", f)
write("#", f)
write("# Unpublished version of October 2025 - kristoffer.johansson@bio.ku.dk", f)
write("# ", f)
write(sprintf("%5s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s", "var", "dmso_psi", "dmso_sd", "ao_psi", "ao_sd", "score", "score_sd","keima","keima_sd"), f)
write(sprintf("%5s  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f",
              prkn_func$var, prkn_func$dmso, prkn_func$dmso_sd, prkn_func$ao, prkn_func$ao_sd,
	      prkn_func$activity, prkn_func$activity_sd, prkn_func$activity_fnorm, prkn_func$activity_fnorm_sd), f)
close(f)

print("")
print(settings)


