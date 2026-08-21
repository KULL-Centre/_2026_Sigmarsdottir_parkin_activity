options(width=160)

aa_one = strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

prkn = read.csv("prkn.csv")

nc = nchar(prkn$var)
prkn$aa = substr(prkn$var, 1, 1)
prkn$resi = as.numeric(substr(prkn$var, 2, nc-1))
prkn$mut = substr(prkn$var, nc, nc)

iwt = which(prkn$var == "WT")
prkn[iwt,"aa"] = ""
prkn[iwt,"resi"] = 0
prkn[iwt,"mut"] = ""

isyn = which(prkn$mut == "=")
inons = which(prkn$mut == "*")
i1 = which(prkn$mut %in% aa_one)

iR42P = which(prkn$var == "R42P")
iC431A = which(prkn$var == "C431A")


# #
# # Distributions and scatter of average gate scores of individual AO and DMSO assays
# #

# plot_distributions = function(ht1, htn, hts, hc1, hcn, hcs, ...) {
#     plot(0,0,col=0, ylab="Counts", ...)
#     lines(hc1$mids, hc1$counts/10, col=1, lty=2)
#     lines(hcn$mids, hcn$counts, col=2, lty=2)
#     lines(hcs$mids, hcs$counts, col=4, lty=2)
#     lines(ht1$mids, ht1$counts/10, col=1, lty=1)
#     lines(htn$mids, htn$counts, col=2, lty=1)
#     lines(hts$mids, hts$counts, col=4, lty=1)
#     legend("topright", c("DMSO ctrl","AO treated"), lty=c(2,1), lwd=1, col=1)
#     legend("topleft", c("Single var. /10","WT syn.","Nonsense"), lty=1, col=c(1,4,2))
# }

# quartz(width=12, height=4)
# par(mfrow=c(1,3), bg="white")

# breaks = 40
# hc1 = hist(prkn[i1,"dmso"], breaks=breaks, plot=F)
# ht1 = hist(prkn[i1,"ao"],   breaks=breaks, plot=F)
# hcs = hist(prkn[isyn,"dmso"], breaks=breaks, plot=F)
# hts = hist(prkn[isyn,"ao"],   breaks=breaks, plot=F)
# hcn = hist(prkn[inons,"dmso"], breaks=breaks, plot=F)
# htn = hist(prkn[inons,"ao"],   breaks=breaks, plot=F)
# plot_distributions(ht1, htn, hts, hc1, hcn, hcs, xlim=c(1,4), xlab="Avg. gate", ylim=c(0,160))

# hc1 = hist(prkn[i1,"funcf_dmso"], breaks=breaks, plot=F)
# ht1 = hist(prkn[i1,"funcf_ao"],   breaks=breaks, plot=F)
# hcs = hist(prkn[isyn,"funcf_dmso"], breaks=breaks, plot=F)
# hts = hist(prkn[isyn,"funcf_ao"],   breaks=breaks, plot=F)
# hcn = hist(prkn[inons,"funcf_dmso"], breaks=breaks, plot=F)
# htn = hist(prkn[inons,"funcf_ao"],   breaks=breaks, plot=F)
# plot_distributions(ht1, htn, hts, hc1, hcn, hcs, xlim=c(0,1), xlab="Keima", ylim=c(0,160))

# hc1 = hist(log10(prkn[i1,"funcf_dmso"]), breaks=breaks, plot=F)
# ht1 = hist(log10(prkn[i1,"funcf_ao"]),   breaks=breaks, plot=F)
# hcs = hist(log10(prkn[isyn,"funcf_dmso"]), breaks=breaks, plot=F)
# hts = hist(log10(prkn[isyn,"funcf_ao"]),   breaks=breaks, plot=F)
# hcn = hist(log10(prkn[inons,"funcf_dmso"]), breaks=breaks, plot=F)
# htn = hist(log10(prkn[inons,"funcf_ao"]),   breaks=breaks, plot=F)
# plot_distributions(ht1, htn, hts, hc1, hcn, hcs, xlim=c(-1.4,0), xlab="log Keima", ylim=c(0,200))

# quartz.save("prkn_func_distributions.png", type="png")

# # scatter plot moved to score analyses
# plot_scatter = function(x, y, i_low, i_high, ...) {
#     plot(x, y, pch=16, cex=.3, ...)
#     points(x[ila],    y[ila],    pch=16, cex=0.4, col=5)
#     points(x[iha],    y[iha],    pch=16, cex=0.4, col=8)    
#     points(x[inons],  y[inons],  pch=16, cex=0.4, col=2)
#     points(x[isyn],   y[isyn],   pch=16, cex=0.4, col=4)
#     points(x[iwt],    y[iwt],    pch=16, cex=1.5, col=3)
#     points(x[iR42P],  y[iR42P],  pch=16, cex=1.5, col=7)
#     points(x[iC431A], y[iC431A], pch=16, cex=1.5, col=6)
#     legend("top", c("WT","C431A","R42P","Single var.","WT syn.","Nonsense","Low-abun","High-abun"), pch=16, col=c(3,6,7,1,4,2,5,8), ncol=4)
# }

# quartz(width=12, height=4)
# par(mfrow=c(1,3), bg="white")
# ila = which(prkn$abun < 0.3 & prkn$mut %in% aa_one)
# iha = which(prkn$abun > 0.7 & prkn$mut %in% aa_one)
# plot_scatter(prkn$ao, prkn$dmso, ila, iha, xlab="AO treated [Avg. gate]", ylab="DMSO Ctrl [Avg. gate]", ylim=c(1.6,4.0))
# plot_scatter(prkn$funcf_ao, prkn$funcf_dmso, ila, iha, xlab="AO treated [Keima acedic/neutral]", ylab="DMSO Ctrl [Keima acedic/neutral]", ylim=c(0.0,1.2))
# plot_scatter(prkn$funcf_ao, prkn$funcf_dmso, ila, iha, xlab="AO treated [log Keima acedic/neutral]", ylab="DMSO Ctrl [log Keima acedic/neutral]", log="xy", xlim=c(0.05,1.0), ylim=c(0.03,2.0))
# quartz.save("prkn_func_scatter.png", type="png")

# quartz(width=6, height=6)
# par(bg="white")
# plot_scatter(prkn$ao, prkn$dmso, ila, iha, xlab="AO treated [Avg. gate]", ylab="DMSO Ctrl [Avg. gate]", ylim=c(1.6,4.0))
# abline(c(0.5,1))
# abline(c(1.5,1/1.77), lty=2)
# legend("bottomright", c("DMSO factor 1.0","DMSO factor 1.77"), lty=c(1,2), col=1)
# quartz.save("prkn_func_scatter_lin.png", type="png")


#
# Activity score by combining average gate scores
#

plot_activity_distribution = function(act_score, xlim=NA, ylim=NA, ...) {
    breaks = 25
    stopifnot( length(act_score) == nrow(prkn) )
    hd1 = hist(act_score[i1], breaks=breaks, plot=F)
    hds = hist(act_score[isyn], breaks=breaks, plot=F)
    hdn = hist(act_score[inons], breaks=breaks, plot=F)
    if (any(is.na(xlim))) { xlim = c(min(act_score,na.rm=T), max(act_score,na.rm=T)) }
    if (any(is.na(ylim))) { ylim = c(0,max(hd1$counts/10)) }
    plot(0,0,col=0, xlim=xlim, ylim=ylim, xlab="Activity score", ylab="Counts", ...)
    lines(hd1$mids, hd1$counts/10, col=1, lwd=1)
    lines(hdn$mids, hdn$counts, col=2, lwd=1)
    lines(hds$mids, hds$counts, col=4, lwd=1)
    points(act_score[iwt], 0.0, pch=16, col=3, cex=2)
    points(act_score[iC431A], 0.0, pch=16, col=6, cex=2)
    points(act_score[iR42P], 0.0, pch=16, col=7, cex=2)
    legend("topleft", c("Single var. /10","WT syn.","Nonsense","WT","C431A","R42P"), pch=c(NA,NA,NA,16,16,16),
           lty=c(1,1,1,NA,NA,NA), col=c(1,4,2,3,6,7))
}

quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")

# ao as score 1?
act_score1 = prkn$ao - prkn$dmso
plot_activity_distribution(act_score1, main="delta-AG w DMSO factor 1.0")

dmso_factor = (prkn[iC431A,"ao"] - mean(prkn[inons,"ao"], na.rm=T)) / (prkn[iC431A,"dmso"] - mean(prkn[inons,"dmso"], na.rm=T))
act_score2 = prkn$ao - dmso_factor*prkn$dmso
plot_activity_distribution(act_score2, main=sprintf("delta-AG w DMSO factor %.3f",dmso_factor))

# act_score3 = prkn$funcf_ao - prkn$funcf_dmso
# plot_activity_distribution(act_score3, main="delta-Keima")

# quartz.save("prkn_ag_lin_score.png", type="png")

# quartz(width=6, height=6)
# act_score4 = prkn$ao - 4.0*prkn$dmso
act_score3 = prkn$act_keima
plot_activity_distribution(act_score3, main="delta Keima")
quartz.save("prkn_ag_lin.png", type="png")


plot_act_fitness = function(x, y, ...) {
    rp = cor(x, y, method="pearson", use="complete.obs")
    plot(x, y, pch=16, cex=.3, main=sprintf("Pearson %.2f",rp), ...)
    points(x[inons],  y[inons],  pch=16, cex=0.4, col=2)
    points(x[isyn],   y[isyn],   pch=16, cex=0.4, col=4)
    points(x[iwt],    y[iwt],    pch=16, cex=1.5, col=3)
    points(x[iR42P],  y[iR42P],  pch=16, cex=1.5, col=7)
    points(x[iC431A], y[iC431A], pch=16, cex=1.5, col=6)
    legend("top", c("Single var.","C431A","R42P"), pch=16, col=c(1,6,7), ncol=3)
}

quartz(width=9, height=6)
par(mfrow=c(2,3), bg="white")
plot_act_fitness(prkn$peve, act_score1, xlab="popEVE", ylab="Activity AO-DMSO [AG]", ylim=c(-1.6,1.5))
plot_act_fitness(prkn$peve, act_score2, xlab="popEVE", ylab="Activity AO-1.77*DMSO [AG]", ylim=c(-4.5,0.4))
plot_act_fitness(prkn$peve, act_score3, xlab="popEVE", ylab="Activity AO-DMSO [log Keima]", ylim=c(-1,4))
plot_act_fitness(prkn$am_score, act_score1, xlab="AlphaMissense", ylab="Activity AO-DMSO [AG]", ylim=c(-1.6,1.5))
plot_act_fitness(prkn$am_score, act_score2, xlab="AlphaMissense", ylab="Activity AO-1.77*DMSO [AG]", ylim=c(-4.5,0.4))
plot_act_fitness(prkn$am_score, act_score3, xlab="AlphaMissense", ylab="Activity AO-DMSO [log Keima]", ylim=c(-1,4))
quartz.save("prkn_fitness_act.png", type="png")

quartz(width=9, height=6)
par(mfrow=c(2,3), bg="white")
plot_act_fitness(prkn$peve,     act_score1+prkn$abun, xlab="popEVE", ylab="Activity+Abundance AO-DMSO [AG]",
                 ylim=c(-1.6,2.5))
plot_act_fitness(prkn$peve,     act_score2+prkn$abun, xlab="popEVE", ylab="Activity+Abundance AO-1.77*DMSO [AG]",
                 ylim=c(-4.5,1.0))
plot_act_fitness(prkn$peve,     act_score3+prkn$abun, xlab="popEVE", ylab="Activity+Abundance AO-DMSO [log Keima]",
                 ylim=c(-1,4))
plot_act_fitness(prkn$am_score, act_score1+prkn$abun, xlab="AlphaMissense", ylab="Activity+Abundance AO-DMSO [AG]",
                 ylim=c(-1.6,2.5))
plot_act_fitness(prkn$am_score, act_score2+prkn$abun, xlab="AlphaMissense", ylab="Activity+Abundance AO-1.77*DMSO [AG]",
                 ylim=c(-4.2,1.1))
plot_act_fitness(prkn$am_score, act_score3+prkn$abun, xlab="AlphaMissense", ylab="Activity+Abundance AO-DMSO [log Keima]",
                 ylim=c(-1,4))
quartz.save("prkn_fitness_act_abun.png", type="png")

 
# quartz(width=16, height=4)
# par(mfrow=c(1,4), bg="white")
# plot_act_fitness(prkn$abun, prkn$dmso, xlab="Abundance [AG]", ylab="DMSO [AG]")
# plot_act_fitness(prkn$abun, prkn$ao, xlab="Abundance [AG]", ylab="AO [AG]")
# plot_act_fitness(prkn$abun, prkn$ao-prkn$dmso, xlab="Abundance [AG]", ylab="AO-DMSO [AG]")
# plot_act_fitness(prkn$abun, prkn$ao-1.77*prkn$dmso, xlab="Abundance [AG]", ylab="AO-1.77*DMSO [AG]")
# quartz.save("prkn_act_abun.png", type="png")


quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")

# i = which((! is.na(prkn$peve)) & (! is.na(prkn$ao)) & (! is.na(prkn$dmso)) & (! is.na(prkn$abun)))
fit = lm(peve ~ ao + dmso + abun, data=prkn)
summary(fit)
coef_dmso = coef(fit)["dmso"]/coef(fit)["ao"]
coef_abun = coef(fit)["abun"]/coef(fit)["ao"]
print(sprintf("Opt. popEVE gives DMSO factor %.2f, Abun. factor %.2f",coef_dmso,coef_abun))
fitness_pe = prkn$ao + coef_dmso*prkn$dmso + coef_abun*prkn$abun
plot_act_fitness(prkn$peve, fitness_pe, ylim=c(0.3,4.5), xlab="popEVE", ylab="Fitness AO-0.41*DMSO+1.72*Abund [AG]")
legend("bottomright", c(sprintf("DMSO %.2f",coef_dmso), sprintf("Abun %.2f",coef_abun)), pch=NA)

# Include nonsense variants with C431A score
peve_nons = prkn$peve
peve_nons[inons] = peve_nons[iC431A]
fit = lm(peve_nons ~ prkn$ao + prkn$dmso + prkn$abun)
summary(fit)
coef_dmso = coef(fit)["prkn$dmso"]/coef(fit)["prkn$ao"]
coef_abun = coef(fit)["prkn$abun"]/coef(fit)["prkn$ao"]
fitness_pen = prkn$ao + coef_dmso*prkn$dmso + coef_abun*prkn$abun
plot_act_fitness(peve_nons, fitness_pen, xlab="popEVE", ylim=c(-3.5,1.8), ylab="Fitness AO-1.53*DMSO+1.46*Abund [AG]")
legend("bottomright", c(sprintf("DMSO %.2f",coef_dmso), sprintf("Abun %.2f",coef_abun)), pch=NA)

fit = lm(peve_nons ~ ao + dmso, data=prkn)
coef_dmso = coef(fit)["dmso"]/coef(fit)["ao"]
rp = cor(peve_nons, prkn$ao+coef_dmso*prkn$dmso, use="complete.obs", method="pearson")
print(sprintf("Opt. popEVE incl nons. wo. abun. gives DMSO factor %.2f and pearson %.2f",coef_dmso,rp))
print(summary(fit))
fitness_pen = prkn$ao + coef_dmso*prkn$dmso
plot_act_fitness(peve_nons, fitness_pen, xlab="popEVE", ylab="Fitness AO-3.95*DMSO [AG]")
legend("bottomright", sprintf("DMSO %.2f",coef_dmso), pch=NA)

quartz.save("prkn_fitness_pe_opt.png", type="png")



quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")

fit = lm(am_score ~ ao + dmso + abun, data=prkn)
summary(fit)
coef_dmso = coef(fit)["dmso"]/coef(fit)["ao"]
coef_abun = coef(fit)["abun"]/coef(fit)["ao"]
print(sprintf("Opt. AlphaM. gives DMSO factor %.2f, Abun. factor %.2f",coef_dmso,coef_abun))
fitness_pe = prkn$ao + coef_dmso*prkn$dmso + coef_abun*prkn$abun
plot_act_fitness(prkn$am_score, fitness_pe, ylim=c(1.5,6), xlab="AlphaMissense", ylab="Fitness AO-DMSO+Abund [AG]")
legend("bottomleft", c(sprintf("DMSO %.2f",coef_dmso), sprintf("Abun %.2f",coef_abun)), pch=NA)

# Include nonsense variants with C431A score
am_nons = prkn$am_score
am_nons[inons] = am_nons[iC431A]
fit = lm(am_nons ~ prkn$ao + prkn$dmso + prkn$abun)
summary(fit)
coef_dmso = coef(fit)["prkn$dmso"]/coef(fit)["prkn$ao"]
coef_abun = coef(fit)["prkn$abun"]/coef(fit)["prkn$ao"]
fitness = prkn$ao + coef_dmso*prkn$dmso + coef_abun*prkn$abun
plot_act_fitness(am_nons, fitness, ylim=c(0,5), xlab="AlphaMissense", ylab="Fitness [AG]")
legend("bottomleft", c(sprintf("DMSO %.2f",coef_dmso), sprintf("Abun %.2f",coef_abun)), pch=NA)

fit = lm(am_nons ~ ao + dmso, data=prkn)
coef_dmso = coef(fit)["dmso"]/coef(fit)["ao"]
print(sprintf("Opt. AlphaM incl nons. wo. abun. gives DMSO factor %.2f",coef_dmso))
print(summary(fit))
fitness = prkn$ao + coef_dmso*prkn$dmso
plot_act_fitness(am_nons, fitness, xlab="AlphaMissense", ylab="Fitness [AG]")
legend("bottomleft", sprintf("DMSO %.2f",coef_dmso), pch=NA)

quartz.save("prkn_fitness_am_opt.png", type="png")


ib = which(prkn$clinvar %in% c("likely_benign","benign"))
ip = which(prkn$clinvar %in% c("likely_pathogenic","pathogenic"))
ig = which(is.na(prkn$clinvar) & ! is.na(prkn$gnomad41))
plot_act_fitness = function(x, y, ...) {
    rp = cor(x, y, method="pearson", use="complete.obs")
    plot(x, y, pch=16, cex=.3, main=sprintf("Pearson %.2f",rp), ...)
    # points(x[inons],  y[inons],  pch=16, cex=0.4, col=2)
    # points(x[isyn],   y[isyn],   pch=16, cex=0.4, col=4)
    points(x[ig],     y[ig],     pch=16, col=3, cex=.6)
    points(x[ib],     y[ib],     pch=16, col=4, cex=1.5)
    points(x[ip],     y[ip],     pch=16, col=2, cex=1.5)
    # points(x[iwt],    y[iwt],    pch=16, cex=1.5, col=3)
    # points(x[iR42P],  y[iR42P],  pch=16, cex=1.5, col=7)
    # points(x[iC431A], y[iC431A], pch=16, cex=1.5, col=6)
    legend("top", c("Single var.","GnomAD","Pathogenic","Benign"), pch=16, col=c(1,3,2,4), ncol=2)
}


fit = lm(peve_nons ~ prkn$activity + prkn$abun)
summary(fit)
coef_abun_pen = coef(fit)["prkn$abun"]/coef(fit)["prkn$activity"]
fitness_pen = prkn$activity + coef_abun_pen*prkn$abun

fit = lm(am_nons ~ prkn$activity + prkn$abun)
summary(fit)
coef_abun_amn = coef(fit)["prkn$abun"]/coef(fit)["prkn$activity"]
fitness_amn = prkn$activity + coef_abun_amn*prkn$abun

quartz(width=10, height=5)
par(mfrow=c(1,2), bg="white")
plot_act_fitness(peve_nons, fitness_pen, xlab="popEVE", ylab=sprintf("Activity + %.2f*abundance [AG]",coef_abun_pen), ylim=c(-1.2,4.5))
plot_act_fitness(am_nons, fitness_amn, xlab="AlphaMissense", ylab=sprintf("Activity + %.2f*abundance [AG]",coef_abun_amn), ylim=c(-1.1,5.5))
quartz.save("prkn_fitness_act_abun.png", type="png")


# fit = glm(lab ~ activity + abun, data=prkn, family=binomial)
# summary(fit)
# coef_abun_lr = coef(fit)["abun"]/coef(fit)["activity"]
# fitness_lr = prkn$activity + coef_abun_lr*prkn$abun
# plot_act_fitness(peve_nons, fitness_lr, xlab="popEVE", ylab=sprintf("Activity + %.2f*abundance [AG]",coef_abun_lr), ylim=c(-1.2,3))


##
## Thresholds for activity
##

plot_threshold = function(cn, ...) {
    plot_activity_distribution(prkn[,cn], ...)
    cut_lof = quantile(prkn[inons,cn], probs=0.95, na.rm=T)
    cut_wtl = quantile(prkn[isyn,cn], probs=0.05, na.rm=T)
    abline(v=c(cut_lof, cut_wtl), lty=2)
    
    ib = which(prkn$clinvar %in% c("likely_benign","benign"))
    ip = which(prkn$clinvar %in% c("likely_pathogenic","pathogenic"))
    points(prkn[ib,cn], 100+rnorm(length(ib),0,10), pch=16, col=4)
    points(prkn[ip,cn], 100+rnorm(length(ip),0,10), pch=16, col=2)
    
    ig = which(is.na(prkn$clinvar) & ! is.na(prkn$gnomad41))
    points(prkn[ig,cn], 50+rnorm(length(ig),0,10), pch=16, col=8, cex=.6)
    legend("topright", c("Pathogenic","Benign","GnomAD",sprintf("LOF %.2f",cut_lof),sprintf("WTL %.2f",cut_wtl)), col=c(2,4,8,1,1),
           pch=c(16,16,16,NA,NA), lty=c(NA,NA,NA,2,2))
}

quartz(width=12, height=6)
par(mfrow=c(1,2), bg="white")
plot_threshold("activity", xlim=c(-0.6,1.8), ylim=c(0,210), main="Average gate")
plot_threshold("act_keima", xlim=c(-1.2,2.3), ylim=c(0,210), main="log Keima")
quartz.save("prkn_activity_thresholds.png", type="png")

# Add functional categories
cut_lof = 0.36
cut_wtl = 0.73
prkn$cat = "Intermediate"
prkn[which(is.na(prkn$activity)),"cat"] = NA
prkn[which(prkn$activity < cut_lof),"cat"] = "LOF"
prkn[which(prkn$activity > cut_wtl),"cat"] = "WTL"
prkn[which(prkn$activity < cut_lof & prkn$activity+prkn$activity_sd > cut_lof),"cat"] = "possibleLOF"
prkn[which(prkn$activity > cut_wtl & prkn$activity-prkn$activity_sd < cut_wtl),"cat"] = "possibleWTL"
prkn$cat = as.factor(prkn$cat)
print(table(prkn$cat))

i1_meas = setdiff(i1, which(is.na(prkn$activity)))
print(sprintf("Single variant percents of %d", length(i1_meas)))
print(table(prkn[i1_meas,"cat"])*100.0/length(i1_meas))

# Add functional categories, Keima based
cutf_lof = 0.52
cutf_wtl = 0.59
prkn$catf = "Intermediate"
prkn[which(is.na(prkn$act_keima)),"catf"] = NA
prkn[which(prkn$act_keima < cutf_lof),"catf"] = "LOF"
prkn[which(prkn$act_keima > cutf_wtl),"catf"] = "WTL"
prkn[which(prkn$act_keima < cutf_lof & prkn$act_keima+prkn$act_keima_sd > cutf_lof),"catf"] = "possibleLOF"
prkn[which(prkn$act_keima > cutf_wtl & prkn$act_keima-prkn$act_keima_sd < cutf_wtl),"catf"] = "possibleWTL"
prkn$catf = as.factor(prkn$catf)
print(table(prkn$catf))

i1_meas = setdiff(i1, which(is.na(prkn$act_keima)))
print(sprintf("Single variant percents of %d", length(i1_meas)))
print(table(prkn[i1_meas,"catf"])*100.0/length(i1_meas))


#
# ROC analysis
#

eval_roc = function(y_true, y_pred, order="from-right", n_thresholds=100, plot=F) {
    # y_true should be two-level factors or bool
    # With order "from-right", y_pred should be likelihood of second level (or TRUE), e.g. type="response"
    stopifnot( order %in% c("from-right","from-left"))
    if (! is.factor(y_true) ) { y_true = as.factor(y_true); stopifnot( length(levels(y_true)) == 2 ) }
    if (any(is.na(y_true)) | any(is.na(y_pred))) { print("ERROR: Data contain NA values"); return(NA) }
    delta_score = max(y_pred) - min(y_pred)
    threshold = seq( floor(min(y_pred)*n_thresholds/delta_score), ceiling(max(y_pred)*n_thresholds/delta_score) ) *delta_score/n_thresholds
    true_mask = y_true == levels(y_true)[2]
    n_pos = sum(true_mask)
    n_neg = length(y_true) - n_pos
    if (order == "from-left") {
        tp = sapply(threshold, function(t){ sum(y_pred <= t & true_mask) })
        fp = sapply(threshold, function(t){ sum(y_pred <= t & ! true_mask) })        
    } else {    # from-right
        threshold = rev(threshold)
        tp = sapply(threshold, function(t){ sum(y_pred > t & true_mask) })
        fp = sapply(threshold, function(t){ sum(y_pred > t & ! true_mask) })
    }
    tpr = tp / n_pos
    fpr = fp / n_neg
    # Area under the curve
    auc = sum( (c(fpr,1)-c(0,fpr))*(c(tpr,1)+c(0,tpr)) )/2
    if (plot) {
        plot(fpr, tpr, type="l", xlab="False positive rate", ylab="True positive rate", main=sprintf("AUC: %.3f",auc))
	it = which.max( tpr - fpr )
        points(fpr[it], tpr[it], pch=16)
        text(fpr[it], tpr[it], sprintf("Threshold %.3f (%.2f,%.2f)", threshold[it], fpr[it], tpr[it]), pos=4, xpd=T)
	legend("bottomright", c(sprintf("Predicting label '%s'",levels(y_true)[2]), sprintf("Predicting '%s'",order),
	                        sprintf("Positive %d, negative %d",n_pos,n_neg)), pch=NA)
    }
    df = data.frame(threshold=threshold, rev_thres=rev(threshold), tpr=tpr, fpr=fpr, specificity=1-fpr, 
                    ppv = tp / (tp+fp), npv = (n_neg-fp) / (n_neg-fp + n_pos-tp), n=tp+fp)
    return(list(auc=auc, pred=y_pred, label=y_true, roc=df))
}

plot_roc = function(scores, score_names, score_colors=NA, ...) {
    if (any(is.na(score_colors))) { score_colors=seq_along(scores) }
    plot(0,0,col=0, xlim=c(0,1), ylim=c(0,1), xlab="False positive rate", ylab="True positive rate", ...)
    legends = c()
    # plot in reverse order to have first scores on top of other curves
    for (ic in rev(seq_along(scores))) {
        stopifnot(length(scores[[ic]]) == length(prkn$lab))
	n_lab = sum(! is.na(prkn$lab))
        i = which((! is.na(scores[[ic]])) & (! is.na(prkn$lab)) )
	roc = eval_roc(prkn[i,"lab"], scores[[ic]][i], plot=F)
	print(sprintf("AUC of %s on %d of %d labels: %.2f",score_names[ic],length(i),n_lab,roc$auc))
        lines(roc$roc[,"fpr"], roc$roc[,"tpr"], col=score_colors[ic])
	legends = c(legends, sprintf("%s AUC %.3f",score_names[ic],roc$auc))
    }
    legend("bottomright", rev(legends), lty=1, col=score_colors)
}

prkn$lab = NA

print("Clinvar lables")
prkn[which(prkn$clinvar %in% c("likely_benign","benign")),"lab"] = "wt-like"
prkn[which(prkn$clinvar %in% c("likely_pathogenic","pathogenic")),"lab"] = "lof"
table(prkn$lab)

quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")
scores      = list(prkn$activity, prkn$act_keima, prkn$ao, act_score2, prkn$abun, prkn$peve, -prkn$am_score)
score_names = c("AO-1.77*DMSO", "Keima", "AO", "AO-DMSO", "Abundance", "popEVE", "AlphaMissense")
plot_roc(scores, score_names, main="ClinVar")

# Gnomad recordings that are not in clinvar assumed to be benign
print("GnamAD lables added")
prkn[which(is.na(prkn$clinvar) & ! is.na(prkn$gnomad41)),"lab"] = "wt-like"
table(prkn$lab)

scores      = list(prkn$activity, prkn$act_keima, prkn$ao, act_score2, prkn$abun, prkn$peve, -prkn$am_score)
score_names = c("AO-1.77*DMSO", "Keima", "AO", "AO-DMSO", "Abundance", "popEVE", "AlphaMissense")
plot_roc(scores, score_names, main="ClinVar and GnomAD")

print("Synonymous and nonsense lables added")
prkn[which(prkn$mut=="="),"lab"] = "wt-like"
prkn[which(prkn$mut=="*"),"lab"] = "lof"
table(prkn$lab)

scores      = list(prkn$activity, prkn$act_keima, prkn$ao, act_score2)
score_names = c("AO-1.77*DMSO", "Keima", "AO", "AO-DMSO")
plot_roc(scores, score_names, main="ClinVar, GnomAD, WT syn. and Nons.")

quartz.save("prkn_roc.png", type="png")


# fit to labels
print("Logistic regression of all labels")
prkn$lab = as.factor(prkn$lab)
lr_fit = glm(lab ~ ao + dmso + abun, data=prkn, family=binomial)
coef_dmso = coef(lr_fit)["dmso"]/coef(lr_fit)["ao"]
coef_abun = coef(lr_fit)["abun"]/coef(lr_fit)["ao"]
quartz(width=5, height=5)
roc = eval_roc(lr_fit$y, lr_fit$fitted.values, plot=T)
legend("bottomright", c(sprintf("Score coef %.2f",coef(lr_fit)["ao"]),sprintf("DMSO coef %.2f",coef_dmso),sprintf("Abundance coef %.2f",coef_abun)), pch=NA)
print(sprintf("Activity and abundance optimized to all labels AUC %.2f. lab ~ sigmoid(%.2f + %.2f*(AO + %.2f*DMSO + %.2f*abund )) (logistic regression)",
              roc$auc, coef(lr_fit)["(Intercept)"], coef(lr_fit)["ao"], coef_dmso, coef_abun))
print(summary(lr_fit))
quartz.save("prkn_roc_act_abun_logi.png", type="png")


# ROC with thresholds marked
cut_lof = 0.36
cut_wtl = 0.73
cutf_lof = 0.52
cutf_wtl = 0.59

quartz(width=10, height=5)
par(mfrow=c(1,2), bg="white")

i = which((! is.na(prkn$activity)) & (! is.na(prkn$lab)) )
lab_fac = factor(prkn[i,"lab"], levels=c("wt-like","lof"))
roc = eval_roc(lab_fac, prkn[i,"activity"], "from-left", plot=T)
it_lof = which.min(abs( roc$roc[,"threshold"]-cut_lof ))
points(roc$roc[it_lof,"fpr"], roc$roc[it_lof,"tpr"], pch=17)
# text(  roc$roc[it_lof,"tpr"], roc$roc[it_lof,"fpr"], sprintf("LOF %.2f (%.3f,%.3f)", cut_lof, fpr[it_lof], tpr[it_lof]), pos=4, xpd=T)
text( 0.1, 0.9, sprintf("LOF threshold %.2f (%.3f,%.3f)", cut_lof, roc$roc[it_lof,"fpr"], roc$roc[it_lof,"tpr"]), pos=4, xpd=T)

i = which((! is.na(prkn$act_keima)) & (! is.na(prkn$lab)) )
lab_fac = factor(prkn[i,"lab"], levels=c("wt-like","lof"))
rocf = eval_roc(lab_fac, prkn[i,"act_keima"], "from-left", plot=T)
itf_lof = which.min(abs( rocf$roc[,"threshold"]-cutf_lof ))
points(rocf$roc[itf_lof,"fpr"], rocf$roc[itf_lof,"tpr"], pch=17)
# text(  rocf$roc[itf_lof,"fpr"], rocf$roc[itf_lof,"tpr"], sprintf("LOF %.2f (%.3f,%.3f)", cut_lof, fpr[itf_lof], tpr[itf_lof]), pos=4, xpd=T)
text( 0.1, 0.8, sprintf("LOF threshold %.2f (%.3f,%.3f)", cutf_lof, rocf$roc[itf_lof,"fpr"], rocf$roc[itf_lof,"tpr"]), pos=4, xpd=T)

quartz.save("prkn_roc_lof_threshold.png", type="png")


quartz(width=10, height=5)
par(mfrow=c(1,2), bg="white")

i = which((! is.na(prkn$activity)) & (! is.na(prkn$lab)) )
lab_fac = factor(prkn[i,"lab"], levels=c("lof","wt-like"))
roc = eval_roc(lab_fac, prkn[i,"activity"], "from-right", plot=T)
it_wtl = which.min(abs( roc$roc[,"threshold"]-cut_wtl ))
points(roc$roc[it_wtl,"fpr"], roc$roc[it_wtl,"tpr"], pch=17)
text(  roc$roc[it_wtl,"fpr"], roc$roc[it_wtl,"tpr"], sprintf("WTL %.2f (%.3f,%.3f)",cut_wtl, roc$roc[it_wtl,"fpr"], roc$roc[it_wtl,"tpr"]), pos=4, xpd=T)

i = which((! is.na(prkn$act_keima)) & (! is.na(prkn$lab)) )
lab_fac = factor(prkn[i,"lab"], levels=c("lof","wt-like"))
rocf = eval_roc(lab_fac, prkn[i,"act_keima"], "from-right", plot=T)
itf_wtl = which.min(abs( rocf$roc[,"threshold"]-cutf_wtl ))
points(rocf$roc[itf_wtl,"fpr"], rocf$roc[itf_wtl,"tpr"], pch=17)
text(rocf$roc[itf_wtl,"fpr"], rocf$roc[itf_wtl,"tpr"],
     sprintf("WTL %.2f (%.3f,%.3f)",cut_wtl, rocf$roc[itf_wtl,"fpr"], rocf$roc[itf_wtl,"tpr"]), pos=4, xpd=T)

quartz.save("prkn_roc_wtl_threshold.png", type="png")


##
## Low throughput validation experiments
##

plot_ltp = function(x, y, dx=NA, dy=NA, ...) {
    rp = cor(x, y, method="pearson", use="complete.obs")
    plot(x, y, pch=16, cex=.8, main=sprintf("Pearson %.2f",rp), ...)
    if (! all(is.na(dx))) { arrows(x0=x-dx, x1=x+dx, y0=y, code=3, angle=90, length=.03) }
    if (! all(is.na(dy))) { arrows(x0=x, y0=y-dy, y1=y+dy, code=3, angle=90, length=.03) }
    points(x[iwt],    y[iwt],    pch=16, cex=1.5, col=3)
    points(x[iR42P],  y[iR42P],  pch=16, cex=1.5, col=7)
    points(x[iC431A], y[iC431A], pch=16, cex=1.5, col=6)
    legend("topleft", c("Single var.","WT","C431A","R42P"), pch=16, col=c(1,3,6,7), ncol=1)
}

quartz(width=10, height=5)
par(mfrow=c(1,2), bg="white")
plot_ltp(prkn$activity, prkn$ltp, prkn$activity_sd, prkn$ltp_sd+1e-3, xlab="HTP activity [AG]", ylab="LTP activity [Keima]", xlim=c(-0.3,1.3))
plot_ltp(prkn$act_keima, prkn$ltp, prkn$act_keima_sd, prkn$ltp_sd+1e-3, xlab="HTP activity [log Keima]", ylab="LTP Activity [Keima]", xlim=c(-0.3,1.3))
quartz.save("prkn_validation_ltp.png", type="png")


##
## Yi et al 2019
##

quartz(width=10, height=5)
par(mfrow=c(1,2), bg="white")
plot_ltp(prkn$activity, prkn$yi19, prkn$activity_sd, xlab="HTP activity [AG]", ylab="Yi19 [Fluorescence]", xlim=c(-0.9,1.5))
plot_ltp(prkn$act_keima, prkn$yi19, prkn$act_keima_sd, xlab="HTP activity [log Keima]", ylab="Yi19 [Fluorescence]", xlim=c(-0.9,1.5))
quartz.save("prkn_validation_yi19.png", type="png")

quartz(width=10, height=5)
par(mfrow=c(1,2), bg="white")
plot_ltp(prkn$activity, log10(prkn$yi19), prkn$activity_sd, xlab="HTP activity [AG]", ylab="Yi19 [log Fluorescence]", xlim=c(-0.9,1.5))
plot_ltp(prkn$act_keima, log10(prkn$yi19), prkn$act_keima_sd, xlab="HTP activity [log Keima]", ylab="Yi19 [log Fluorescence]", xlim=c(-0.9,1.5))
quartz.save("prkn_validation_yi19_log.png", type="png")

# quartz(width=5, height=5)
# plot_ltp(10^prkn$act_keima, prkn$yi19, xlab="HTP activity [Keima]", ylab="Yi19 [Fluorescence]", xlim=c(0,35))


##
## Functional siets
##

plot_scatter_ctrl = function(x, y, i_low, i_high, ...) {
    plot(x, y, pch=16, cex=.3, ...)
    points(x[inons],  y[inons],  pch=16, cex=0.4, col=2)
    # points(x[isyn],   y[isyn],   pch=16, cex=0.4, col=4)
    points(x[iwt],    y[iwt],    pch=16, cex=1.5, col=3)
    points(x[iR42P],  y[iR42P],  pch=16, cex=1.5, col=7)
    points(x[iC431A], y[iC431A], pch=16, cex=1.5, col=6)
    legend("topright", c("WT","C431A","R42P","Single var.","Nons."), pch=16, col=c(3,6,7,1,2), ncol=3)
}

quartz(width=10, height=5)
par(mfrow=c(1,2), bg="white")
plot_scatter_ctrl(prkn$abun, prkn$activity,  ylab="HTP activity [AG]",        xlab="Abundance [AG]", ylim=c(-1.5,3))
plot_scatter_ctrl(prkn$abun, prkn$act_keima, ylab="HTP activity [log Keima]", xlab="Abundance [AG]", ylim=c(-1.5,4))
quartz.save("prkn_act_abun_ctrl.png", type="png")


ib = which(prkn$clinvar %in% c("likely_benign","benign"))
ip = which(prkn$clinvar %in% c("likely_pathogenic","pathogenic"))

plot_scatter_clin = function(x, y, i_low, i_high, ...) {
    plot(x, y, pch=16, cex=.3, ...)
    points(x[ib], y[ib], pch=16, cex=1, col=4)
    points(x[ip], y[ip], pch=16, cex=1, col=2)
    points(x[iwt],    y[iwt],    pch=16, cex=1.5, col=3)
    points(x[iR42P],  y[iR42P],  pch=16, cex=1.5, col=7)
    points(x[iC431A], y[iC431A], pch=16, cex=1.5, col=6)
    legend("topright", c("WT","C431A","R42P","Single var.","Benign","Pathogenic"), pch=16, col=c(3,6,7,1,4,2), ncol=3)
}

quartz(width=10, height=5)
par(mfrow=c(1,2), bg="white")
plot_scatter_clin(prkn$abun, prkn$activity,  ylab="HTP activity [AG]",        xlab="Abundance [AG]", ylim=c(-1.5,3))
plot_scatter_clin(prkn$abun, prkn$act_keima, ylab="HTP activity [log Keima]", xlab="Abundance [AG]", ylim=c(-1.5,4))
quartz.save("prkn_act_abun_clin.png", type="png")


##
## Controls
##

# add activity scores here
print("Frequent pathogenic variants")
cns = c("var","gnomad21","gnomad41","clinvar","activity","activity_sd","act_keima","act_keima_sd","abun","abun_sd","rasa","gemme","rosetta",
        "am_score","am_class","peve")
print(prkn[which(prkn$var %in% c("R275W","G430D")),cns])


##
## Per residue
##

prknr = read.csv("../scores/prkn_activity_residues.csv")
prkn$act_med = prknr[match(prkn$resi,prknr$resi),"med_score"]

# make sure to use the same data in the following
ii = which(! is.na(prkn[i1,"activity"]))
i = i1[ii]

# The explained variance is R^2 and for a linear model with intersept this is rp^2
rp = cor(prkn[i,"activity"], prkn[i,"act_med"], method="pearson")

# a model of variant score based on position medians
fit = lm(prkn[i,"activity"] ~ prkn[i,"act_med"])
act_pred = coef(fit)[1] + coef(fit)[2] * prkn[i,"act_med"]
# Fraction of unexplained variance is residual variance divided by original variance
fuv = var(prkn[i,"activity"]-act_pred)/var(prkn[i,"activity"])

print(sprintf("Fraction of unexplained variance %.1f%% (ratio of residual to total var) and explained %.1f%% (Pearson %.2f squared). Sum %.3f%%",
              fuv*100, rp*rp*100, rp, fuv*100+rp*rp*100))

quartz(width=8,height=4)
par(mfrow=c(1,2), mar=c(5,4,2,2)+.1)
plot(prkn[i,"act_med"], prkn[i,"activity"])
abline(coef(fit), col=2)
legend("topleft", sprintf("Pearson %.2f",rp), pch=NA)

ii = which(prkn[i1,"mut"]=="A" & ! is.na(prkn[i1,"activity"]))
i = i1[ii]
plot(prkn[i,"rasa"], prkn[i,"act_med"])
rp = cor(prkn[i,"rasa"], prkn[i,"act_med"], method="pearson")
rs = cor(prkn[i,"rasa"], prkn[i,"act_med"], method="spearman")
legend("bottomright", c(sprintf("Pearson %.2f",rp),sprintf("Spearman %.2f",rs)), pch=NA)

