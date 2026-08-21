options(width=160)

aa_one = strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

##
## Data
##
load("prkn_activity.rda")
pf = raw_psi
colnames(pf) = c("var",colnames(pf)[2:ncol(pf)])

nc = nchar(pf$var)
pf$aa = substr(pf$var, 1, 1)
pf$resi = as.numeric(substr(pf$var, 2, nc-1))
pf$mut = substr(pf$var, nc, nc)

iwt = which(pf$var == "WT")
pf[iwt,"aa"] = ""
pf[iwt,"resi"] = 0
pf[iwt,"mut"] = ""

isyn = which(pf$mut == "=")
inons = which(pf$mut == "*")
i1 = which(pf$mut %in% aa_one)

# Control variant, low abundance, residual function when overexpressed
iR42P = which(pf$var == "R42P")

# Control catalytic site variant, high abundance, no activity
iC431A = which(pf$var == "C431A")

# Abundance data
pa = read.csv2("../../project_vamp/call_jan2022/vamp_parkin.csv")
ia2f = match(pf$var,pa$var)
pf$abun = as.numeric( pa[ia2f,"vamp_score"] )
pf$abun_sd = as.numeric( pa[ia2f,"vamp_std"] )


##
## Scores and distributions
##
breaks = 40
hc1 = hist(pf[i1,"dmso_mean"], breaks=breaks, plot=F)
ht1 = hist(pf[i1,"ao_mean"],   breaks=breaks, plot=F)
hcs = hist(pf[isyn,"dmso_mean"], breaks=breaks, plot=F)
hts = hist(pf[isyn,"ao_mean"],   breaks=breaks, plot=F)
hcn = hist(pf[inons,"dmso_mean"], breaks=breaks, plot=F)
htn = hist(pf[inons,"ao_mean"],   breaks=breaks, plot=F)

quartz(width=12, height=6)
par(mfrow=c(1,2), bg="white")
plot(0,0,col=0, xlim=c(1,4), ylim=c(0,160), xlab="PSI", ylab="Counts")
lines(hc1$mids, hc1$counts/10, col=1, lty=2)
lines(hcn$mids, hcn$counts, col=2, lty=2)
lines(hcs$mids, hcs$counts, col=4, lty=2)
lines(ht1$mids, ht1$counts/10, col=1, lty=1)
lines(htn$mids, htn$counts, col=2, lty=1)
lines(hts$mids, hts$counts, col=4, lty=1)
legend("topright", c("DMSO ctrl","AO treated"), lty=c(2,1), lwd=1, col=1)
legend("topleft", c("Single var. /10","WT syn.","Nonsense"), lty=1, col=c(1,4,2))

plot(pf$ao_mean, pf$dmso_mean, pch=16, cex=.3, xlab="AO treated", ylab="DMSO Ctrl", ylim=c(1.6,4.0))

ila = which(pf$abun < 0.3 & pf$mut %in% aa_one)
points(pf[ila,"ao_mean"], pf[ila,"dmso_mean"], pch=16, cex=.4, col=5)
iha = which(pf$abun > 0.7 & pf$mut %in% aa_one)
points(pf[iha,"ao_mean"], pf[iha,"dmso_mean"], pch=16, cex=.4, col=8)

points(pf[inons,"ao_mean"], pf[inons,"dmso_mean"], pch=16, cex=.4, col=2)
points(pf[isyn,"ao_mean"], pf[isyn,"dmso_mean"], pch=16, cex=.4, col=4)
points(pf[iwt,"ao_mean"], pf[iwt,"dmso_mean"], pch=16, cex=1.5, col=3)
points(pf[iR42P,"ao_mean"], pf[iR42P,"dmso_mean"], pch=16, cex=1.5, col=7)
points(pf[iC431A,"ao_mean"], pf[iC431A,"dmso_mean"], pch=16, cex=1.5, col=6)
legend("top", c("WT","C431A","R42P","Single var.","WT syn.","Nonsense","Low-abun","High-abun"), pch=16, col=c(3,6,7,1,4,2,5,8), ncol=4)

quartz.save("prkn_ao_dmso.png", type="png")


# Keima AO, DMSO scatter plot
quartz(width=12, height=6)
par(mfrow=c(1,2), bg="white")
plot(10^pf$ao_fnorm, 10^pf$dmso_fnorm, pch=16, cex=.3, xlab="AO Keima", ylab="DMSO Keima")
points(10^pf[ila,"ao_fnorm"],    10^pf[ila,"dmso_fnorm"],    pch=16, cex=.4, col=5)
points(10^pf[iha,"ao_fnorm"],    10^pf[iha,"dmso_fnorm"],    pch=16, cex=.4, col=8)
points(10^pf[inons,"ao_fnorm"],  10^pf[inons,"dmso_fnorm"],  pch=16, cex=.4, col=2)
points(10^pf[isyn,"ao_fnorm"],   10^pf[isyn,"dmso_fnorm"],   pch=16, cex=.4, col=4)
points(10^pf[iwt,"ao_fnorm"],    10^pf[iwt,"dmso_fnorm"],    pch=16, cex=1.5, col=3)
points(10^pf[iR42P,"ao_fnorm"],  10^pf[iR42P,"dmso_fnorm"],  pch=16, cex=1.5, col=7)
points(10^pf[iC431A,"ao_fnorm"], 10^pf[iC431A,"dmso_fnorm"], pch=16, cex=1.5, col=6)
# legend("top", c("WT","C431A","R42P","Single var.","WT syn.","Nonsense","Low-abun","High-abun"), pch=16, col=c(3,6,7,1,4,2,5,8), ncol=4)

plot(pf$ao_fnorm, pf$dmso_fnorm, pch=16, cex=.3, xlim=c(-1.3,0.4), ylim=c(-1.7,0.4), xlab="AO Keima, log", ylab="DMSO Keima, log")
points(pf[ila,"ao_fnorm"], pf[ila,"dmso_fnorm"], pch=16, cex=.4, col=5)
points(pf[iha,"ao_fnorm"], pf[iha,"dmso_fnorm"], pch=16, cex=.4, col=8)
points(pf[inons,"ao_fnorm"], pf[inons,"dmso_fnorm"], pch=16, cex=.4, col=2)
points(pf[isyn,"ao_fnorm"], pf[isyn,"dmso_fnorm"], pch=16, cex=.4, col=4)
points(pf[iwt,"ao_fnorm"], pf[iwt,"dmso_fnorm"], pch=16, cex=1.5, col=3)
points(pf[iR42P,"ao_fnorm"], pf[iR42P,"dmso_fnorm"], pch=16, cex=1.5, col=7)
points(pf[iC431A,"ao_fnorm"], pf[iC431A,"dmso_fnorm"], pch=16, cex=1.5, col=6)
legend("bottom", c("WT","C431A","R42P","Single var.","WT syn.","Nonsense","Low-abun","High-abun"), pch=16, col=c(3,6,7,1,4,2,5,8), ncol=4)
quartz.save("prkn_fnorm_scatter.png", type="png")


##
## Delta-PSI 
##

plot_distr = function(x, xlim=NA, ylim=NA, ylab="Density", ...) {
    xmin = min(x, na.rm=T); xmax = max(x, na.rm=T)
    if (all(is.na(xlim))) { xlim = c(xmin,xmax)}
    breaks = seq(xmin-1e-3, xmax+1e-3, length.out=50)
    h1 = hist(x[i1],   breaks=breaks, plot=F)
    hs = hist(x[isyn],   breaks=breaks, plot=F)
    hn = hist(x[inons],   breaks=breaks, plot=F)
    if (all(is.na(ylim))) { ylim = c(0, max(c(h1$density,hn$density,hs$density))*1.2)}
    plot(0,0,col=0, xlim=xlim, ylim=ylim, ylab=ylab, ...)
    lines(h1$mids, h1$density, col=1)
    lines(hn$mids, hn$density, col=2)
    lines(hs$mids, hs$density, col=4)
    points(x[iR42P], 0, pch=16, col=7, cex=2)
    points(x[iwt], 0, pch=16, col=4, cex=2)
    points(x[iC431A], 0, pch=16, col=6, cex=2)
    legend("topleft", c("Lib. single","Lib. WT syn.","Lib. Nons.","Lib. WT","Lib. C431A","Lib. R42P"), lty=c(1,1,1,NA,NA,NA),
           pch=c(NA,NA,NA,16,16,16), col=c(1,4,2,4,6,7), ncol=2)
}


# quartz(width=12, height=6)
# par(mfrow=c(1,2), bg="white")

# dmso_factor = 1.0
# # dmso_factor = 4.0
# delta_psi = pf$ao_mean - dmso_factor*pf$dmso_mean
# plot_distr(delta_psi, xlab="AO - a*DMSO [avg. gate]", main=sprintf("DMSO factor a=%.2f",dmso_factor))

# dmso_factor = (pf[iC431A,"ao_mean"] - mean(pf[inons,"ao_mean"])) / (pf[iC431A,"dmso_mean"] - mean(pf[inons,"dmso_mean"]))
# delta_psi = pf$ao_mean - dmso_factor*pf$dmso_mean
# plot_distr(delta_psi, xlab="AO - a*DMSO [avg. gate]", main=sprintf("DMSO factor a=%.2f",dmso_factor))

# quartz.save("prkn_delta_score.png", type="png")

quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")

plot_distr(pf$ao_mean, xlab="AO [avg. gate]", main="AO")

dmso_factor = 1.0
# dmso_factor = 4.0
delta_psi = pf$ao_mean - dmso_factor*pf$dmso_mean
plot_distr(delta_psi, xlab="AO - a*DMSO [avg. gate]", main=sprintf("DMSO factor a=%.2f",dmso_factor))

dmso_factor = (pf[iC431A,"ao_mean"] - mean(pf[inons,"ao_mean"])) / (pf[iC431A,"dmso_mean"] - mean(pf[inons,"dmso_mean"]))
delta_psi = pf$ao_mean - dmso_factor*pf$dmso_mean
plot_distr(delta_psi, xlab="AO - a*DMSO [avg. gate]", main=sprintf("DMSO factor a=%.2f",dmso_factor))

quartz.save("prkn_delta_score.png", type="png")


##
## Scatter plots AO vs DMSO
##
# from anno/analyses that no longer has individual AO and DMSO transforms - I use prkn_func_scatter_lin.png!
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


##
## Nonsense variants in AO Keima is bimodal?!
##
# Robust regression (towards outliers) in R: MASS::rlm
# Passing Bablok and Deming regression (loss symm. in y and x) in R: (MethComp and) mcr packages
#   https://www.r-bloggers.com/2015/09/deming-and-passing-bablok-regression-in-r/
library("mcr")

add_reg = function(x,y) {
    fit1 = lm(y ~ x)
    abline(coef(fit1), lty=3)
    fit2 = lm(x ~ y)
    abline(c(-coef(fit2)[1]/coef(fit2)[2], 1/coef(fit2)[2] ), lty=2)
    # we are measuring the same thing (inactivity of nonsence) in two different ways (conditions)
    fit3 = mcreg(x, y, method.reg="Deming")
    # fit3 = mcreg(x, y, method.reg="PaBa") # robust version based on Kendall rank-correlation
    abline(coef(fit3)[1:2], lty=1)
    legend("topleft", c(sprintf("Fit DMSO~AO 1/slope = %.2f", 1/coef(fit1)[2]),
                        sprintf("Fit AO~DMSO slope = %.2f", coef(fit2)[2]),
 			sprintf("Daming fit slope = %.2f", 1/coef(fit3)[2])), lty=c(3,2,1))
}

in_high_ao = which(pf$mut=='*' & pf$ao_fnorm > (pf[iwt,"ao_fnorm"]+0.05))
in_low_ao  = which(pf$mut=='*' & pf$ao_fnorm < (pf[iwt,"ao_fnorm"]-0.05))

# pymol select
paste(pf[in_high_ao,"resi"], collapse="+")
paste(pf[in_low_ao,"resi"], collapse="+")

quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")

plot(pf$abun, pf$ao_mean, xlab="Abundance", ylab="AO avg. gate")
points(pf[in_high_ao,"abun"], pf[in_high_ao,"ao_mean"], pch=16, col=4)
points(pf[in_low_ao,"abun"],  pf[in_low_ao,"ao_mean"],  pch=16, col=2)

rp = cor(pf[inons,"ao_mean"], pf[inons,"dmso_mean"], method="pearson", use="complete.obs")
plot(pf$ao_mean, pf$dmso_mean, pch=16, cex=.3, xlab="AO avg. gate", ylab="DMSO avg. gate", main=sprintf("Nonsense Pearson %.2f",rp))
points(pf[in_high_ao,"ao_mean"], pf[in_high_ao, "dmso_mean"], pch=16, cex=.6, col=4)
points(pf[in_low_ao,"ao_mean"], pf[in_low_ao,"dmso_mean"], pch=16, cex=.6, col=2)
points(pf[isyn,"ao_mean"], pf[isyn,"dmso_mean"], pch=16, cex=.4, col=5)
points(pf[iwt,"ao_mean"], pf[iwt,"dmso_mean"], pch=16, cex=1.5, col=3)
points(pf[iR42P,"ao_mean"], pf[iR42P,"dmso_mean"], pch=16, cex=1.5, col=7)
points(pf[iC431A,"ao_mean"], pf[iC431A,"dmso_mean"], pch=16, cex=1.5, col=6)
legend("bottom", c("WT","C431A","R42P","WT syn.","Nons low AO","Nons high AO"), pch=16, col=c(3,6,7,5,2,4), ncol=3)

add_reg(pf[inons,"ao_mean"], pf[inons,"dmso_mean"])

rp = cor(pf[inons,"ao_fnorm"], pf[inons,"dmso_fnorm"], method="pearson", use="complete.obs")
plot(pf$ao_fnorm, pf$dmso_fnorm, pch=16, cex=.3, xlim=c(-1.3,0.4), ylim=c(-1.7,0.4), xlab="AO Keima, log", ylab="DMSO Keima, log",
     main=sprintf("Nonsense Pearson %.2f",rp))
points(pf[in_high_ao,"ao_fnorm"], pf[in_high_ao,"dmso_fnorm"], pch=16, cex=.6, col=4)
points(pf[in_low_ao,"ao_fnorm"], pf[in_low_ao,"dmso_fnorm"], pch=16, cex=.6, col=2)
points(pf[isyn,"ao_fnorm"], pf[isyn,"dmso_fnorm"], pch=16, cex=.4, col=5)
points(pf[iwt,"ao_fnorm"], pf[iwt,"dmso_fnorm"], pch=16, cex=1.5, col=3)
points(pf[iR42P,"ao_fnorm"], pf[iR42P,"dmso_fnorm"], pch=16, cex=1.5, col=7)
points(pf[iC431A,"ao_fnorm"], pf[iC431A,"dmso_fnorm"], pch=16, cex=1.5, col=6)
legend("bottom", c("WT","C431A","R42P","WT syn.","Nons low AO","Nons high AO"), pch=16, col=c(3,6,7,5,2,4), ncol=3)

add_reg(pf[inons,"ao_fnorm"], pf[inons,"dmso_fnorm"])

quartz.save("nons_ao.png", type="png")


##
## Abundance
##

plot_abun_scat = function(x, y, ...) {
    # abundance data does not have WT syn data
    rp = cor(x, y, method="pearson", use="complete.obs")
    plot(x, y, pch=16, cex=.3, main=sprintf("Pearson %.2f",rp), ...)
    points(x[inons],  y[inons],  pch=16, cex=0.4, col=2)
    points(x[iwt],    y[iwt],    pch=16, cex=1.5, col=3)
    points(x[iR42P],  y[iR42P],  pch=16, cex=1.5, col=7)
    points(x[iC431A], y[iC431A], pch=16, cex=1.5, col=6)
    legend("top", c("Single var.","Nonsense","WT","C431A","R42P"), pch=16, col=c(1,2,3,6,7), ncol=3)    
}

quartz(width=16, height=4)
par(mfrow=c(1,4), bg="white")
plot_abun_scat(pf$abun, pf$dmso_mean, xlab="Abundance [AG]", ylab="DMSO [AG]")
plot_abun_scat(pf$abun, pf$ao_mean, xlab="Abundance [AG]", ylab="AO [AG]")
plot_abun_scat(pf$abun, pf$activity, xlab="Abundance [AG]", ylab="Activity [AG]")
plot_abun_scat(pf$abun, pf$activity_fnorm, xlab="Abundance [AG]", ylab="Activity [Keima,  log]")
quartz.save("prkn_act_abun.png", type="png")


##
## Representation among cells
##

quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")

breaks = seq(0,300,10)
ha = hist(pf[2:nrow(pf),"ao_cpm"], breaks=breaks, plot=F)
hd = hist(pf[2:nrow(pf),"dmso_cpm"], breaks=breaks, plot=F)

ymax = max(c(ha$density,hd$density), na.rm=T)
plot(0,0,col=0, xlim=c(0,300), ylim=c(0,ymax), xlab="Representation among cells [CPM]", ylab="Density", main=sprintf("%d variants",nrow(pf)))
lines(ha$mids, ha$density, col=3)
lines(hd$mids, hd$density, col=4)
legend("topright", c(sprintf("AO, WT %.1f%%", pf[1,"ao_cpm"]*1e-4),
                     sprintf("DMSO, WT %.1f%%", pf[1,"dmso_cpm"]*1e-4)), lty=1, col=c(3,4))

rp = cor(pf[2:nrow(pf),"ao_cpm"], pf[2:nrow(pf),"dmso_cpm"], method="pearson", use="complete.obs")
plot(pf[2:nrow(pf),"ao_cpm"], pf[2:nrow(pf),"dmso_cpm"], xlab="AO [CPM]", ylab="DMSO [CPM]", pch=16, cex=.3, main=sprintf("Pearson %.2f",rp))
abline(c(0,1), col=2)

rp = cor(pf[2:nrow(pf),"ao_cpm"], pf[2:nrow(pf),"abun"], method="pearson", use="complete.obs")
plot(pf[2:nrow(pf),"ao_cpm"], pf[2:nrow(pf),"abun"], ylim=c(0,1.6), xlab="AO [CPM]", ylab="Abundance score", pch=16, cex=.3, main=sprintf("Pearson %.2f",rp))
points(pf[inons,"ao_cpm"], pf[inons,"abun"], pch=16, cex=.4, col=2)
points(pf[iR42P,"ao_cpm"], pf[iR42P,"abun"], pch=16, cex=1.5, col=7)
points(pf[iC431A,"ao_cpm"], pf[iC431A,"abun"], pch=16, cex=1.5, col=6)
legend("top", c("C431A","R42P","Nons."), pch=16, col=c(6,7,2), ncol=3)

quartz.save("cells.png", type="png")


##
## Uncertainty
##

# Relation between score and uncertainty
quartz(width=8, height=8)
par(mfrow=c(2,2), bg="white")
plot(pf$ao_mean,        pf$ao_sd,             pch=16, cex=.5, main="AO", xlab="Score", ylab="Std. dev")
plot(pf$dmso_mean,      pf$dmso_sd,           pch=16, cex=.5, main="DMSO", xlab="Score", ylab="Std. dev")
plot(pf$activity,       pf$activity_sd,       pch=16, cex=.5, main="Activity", xlab="Score", ylab="Std. dev")
plot(pf$activity_fnorm, pf$activity_fnorm_sd, pch=16, cex=.5, main="Keima activity", xlab="Score", ylab="Std. dev")
quartz.save("score_std.png", type="png")

# Relation between activity uncertainty and number of replicates
quartz(width=8, height=8)
par(mfrow=c(2,2), bg="white")
plot(raw_psi$ao_sd,       rnorm(nrow(raw_psi),0,0.1)+raw_psi$ao_rep,   ylim=c(0.5,3.5), pch=16, cex=.3,
     xlab="Std. dev", ylab="Number of AO replicas", main="AO")
plot(raw_psi$dmso_sd,     rnorm(nrow(raw_psi),0,0.1)+raw_psi$dmso_rep, ylim=c(0.5,3.5), pch=16, cex=.3,
     xlab="Std. dev", ylab="Number of DMSO replicas", main="DMSO")
plot(raw_psi$activity_sd, rnorm(nrow(raw_psi),0,0.1)+raw_psi$ao_rep,   ylim=c(0.5,3.5), pch=16, cex=.3,
     xlab="Std. dev", ylab="Number of AO replicas", main="Activity")
plot(raw_psi$activity_sd, rnorm(nrow(raw_psi),0,0.1)+raw_psi$dmso_rep, ylim=c(0.5,3.5), pch=16, cex=.3,
     xlab="Std. dev", ylab="Number of DMSO replicas", main="Activity")
quartz.save("sd_rep.png", type="png")

# Mean and min number of read counts among replicates
ao_mean   = apply(raw_psi[,c("ao_1_read_sum","ao_2_read_sum","ao_3_read_sum")],       MARGIN=1, mean, na.rm=T)
dmso_mean = apply(raw_psi[,c("dmso_1_read_sum","dmso_2_read_sum","dmso_3_read_sum")], MARGIN=1, mean, na.rm=T)
ao_min    = apply(raw_psi[,c("ao_1_read_sum","ao_2_read_sum","ao_3_read_sum")],       MARGIN=1, min, na.rm=T)
dmso_min  = apply(raw_psi[,c("dmso_1_read_sum","dmso_2_read_sum","dmso_3_read_sum")], MARGIN=1, min, na.rm=T)

# Test read count threshold to plot
cut = 100

# Relation between replicate read counts and std. dev. of activity scores
quartz(width=8, height=8)
par(mfrow=c(2,2), bg="white")
plot(raw_psi$activity_sd, ao_mean+1,   pch=16, cex=.3, log="y", main="AO mean", xlab="Activity std. dev.", ylab="Replica read counts +1")
abline(h=c(cut,settings$threshold_counts_per_rep), lty=c(3,2))
legend("topright", c(sprintf("Test thres. %d",cut),sprintf("Applied thres. %d",settings$threshold_counts_per_rep)), lty=c(3,2))
plot(raw_psi$activity_sd, dmso_mean+1, pch=16, cex=.3, log="y", main="DMSO mean", xlab="Activity std. dev.", ylab="Replica read counts +1")
abline(h=c(cut,settings$threshold_counts_per_rep), lty=c(3,2))
legend("topright", c(sprintf("Test thres. %d",cut),sprintf("Applied thres. %d",settings$threshold_counts_per_rep)), lty=c(3,2))
plot(raw_psi$activity_sd, ao_min+1,    pch=16, cex=.3, log="y", main="AO min", xlab="Activity std. dev.", ylab="Replica read counts +1")
abline(h=c(cut,settings$threshold_counts_per_rep), lty=c(3,2))
legend("topright", c(sprintf("Test thres. %d",cut),sprintf("Applied thres. %d",settings$threshold_counts_per_rep)), lty=c(3,2))
plot(raw_psi$activity_sd, dmso_min+1,  pch=16, cex=.3, log="y", main="DMSO min", xlab="Activity std. dev.", ylab="Replica read counts +1")
abline(h=c(cut,settings$threshold_counts_per_rep), lty=c(3,2))
legend("topright", c(sprintf("Test thres. %d",cut),sprintf("Applied thres. %d",settings$threshold_counts_per_rep)), lty=c(3,2))
quartz.save("sd_act_counts.png", type="png")

# Relation between replicate read counts and std. dev. of AO and DMSO scores
quartz(width=8, height=8)
par(mfrow=c(2,2), bg="white")
plot(raw_psi$ao_sd,   ao_mean+1,   pch=16, cex=.3, log="y", main="AO mean", xlab="AO std. dev.", ylab="Replica read counts +1")
abline(h=c(cut,settings$threshold_counts_per_rep), lty=c(3,2))
legend("topright", c(sprintf("Test thres. %d",cut),sprintf("Applied thres. %d",settings$threshold_counts_per_rep)), lty=c(3,2))
plot(raw_psi$dmso_sd, dmso_mean+1, pch=16, cex=.3, log="y", main="DMSO mean", xlab="DMSO std. dev.", ylab="Replica read counts +1")
abline(h=c(cut,settings$threshold_counts_per_rep), lty=c(3,2))
legend("topright", c(sprintf("Test thres. %d",cut),sprintf("Applied thres. %d",settings$threshold_counts_per_rep)), lty=c(3,2))
plot(raw_psi$ao_sd,   ao_min+1,    pch=16, cex=.3, log="y", main="AO min", xlab="AO std. dev.", ylab="Replica read counts +1")
abline(h=c(cut,settings$threshold_counts_per_rep), lty=c(3,2))
legend("topright", c(sprintf("Test thres. %d",cut),sprintf("Applied thres. %d",settings$threshold_counts_per_rep)), lty=c(3,2))
plot(raw_psi$dmso_sd, dmso_min+1,  pch=16, cex=.3, log="y", main="DMSO min", xlab="DMSO std. dev.", ylab="Replica read counts +1")
abline(h=c(cut,settings$threshold_counts_per_rep), lty=c(3,2))
legend("topright", c(sprintf("Test thres. %d",cut),sprintf("Applied thres. %d",settings$threshold_counts_per_rep)), lty=c(3,2))
quartz.save("sd_ao_dmso_counts.png", type="png")

# Relations between activity avg. gate and Keima transform with uncertainties as error bars
quartz(width=6, height=6)
plot(0,0,col=0, xlim=c(1,4), ylim=c(-1.8,1.5), xlab="Average gate", ylab="Keima, log10", main="AO scores")
arrows(x0=pf$ao_mean, y0=pf$ao_fnorm-pf$ao_fnorm_sd, y1=pf$ao_fnorm+pf$ao_fnorm_sd, code=3, angle=90, length=0)
arrows(x0=pf$ao_mean-pf$ao_sd, x1=pf$ao_mean+pf$ao_sd, y0=pf$ao_fnorm, code=3, angle=90, length=0)
points(pf$ao_mean, pf$ao_fnorm, pch=16, cex=.5, col=2)
points(pf[iwt,"ao_mean"], pf[iwt,"ao_fnorm"], pch=16, col=3, cex=1.5)
points(pf[iC431A,"ao_mean"], pf[iC431A,"ao_fnorm"], pch=16, col=6, cex=1.5)
points(pf[iR42P,"ao_mean"], pf[iR42P,"ao_fnorm"], pch=16, col=7, cex=1.5)
legend("topleft", c("WT","C431A","R42P"), pch=16, col=c(3,6,7), ncol=1)
quartz.save("prkn_ao_error_prop.png", type="png")

quartz(width=6, height=6)
plot(0,0,col=0, xlim=c(1,4), ylim=c(-1.8,1.5), xlab="Average gate", ylab="Keima, log10", main="DMSO scores")
arrows(x0=pf$dmso_mean, y0=pf$dmso_fnorm-pf$dmso_fnorm_sd, y1=pf$dmso_fnorm+pf$dmso_fnorm_sd, code=3, angle=90, length=0)
arrows(x0=pf$dmso_mean-pf$dmso_sd, x1=pf$dmso_mean+pf$dmso_sd, y0=pf$dmso_fnorm, code=3, angle=90, length=0)
points(pf$dmso_mean, pf$dmso_fnorm, pch=16, cex=.5, col=2)
points(pf[iwt,"dmso_mean"], pf[iwt,"dmso_fnorm"], pch=16, col=3, cex=1.5)
points(pf[iC431A,"dmso_mean"], pf[iC431A,"dmso_fnorm"], pch=16, col=6, cex=1.5)
points(pf[iR42P,"dmso_mean"], pf[iR42P,"dmso_fnorm"], pch=16, col=7, cex=1.5)
legend("topleft", c("WT","C431A","R42P"), pch=16, col=c(3,6,7), ncol=1)
quartz.save("prkn_dmso_error_prop.png", type="png")

quartz(width=6, height=6)
plot(0,0,col=0, xlim=c(-2,3), ylim=c(-3,4), xlab="Average gate", ylab="Keima, log10", main="Activity scores")
arrows(x0=pf$activity, y0=pf$activity_fnorm-pf$activity_fnorm_sd, y1=pf$activity_fnorm+pf$activity_fnorm_sd, code=3, angle=90, length=0)
arrows(x0=pf$activity-pf$activity_sd, x1=pf$activity+pf$activity_sd, y0=pf$activity_fnorm, code=3, angle=90, length=0)
points(pf$activity, pf$activity_fnorm, pch=16, cex=.3, col=2)
points(pf[iwt,"activity"], pf[iwt,"activity_fnorm"], pch=16, col=3, cex=1.5)
points(pf[iC431A,"activity"], pf[iC431A,"activity_fnorm"], pch=16, col=6, cex=1.5)
points(pf[iR42P,"activity"], pf[iR42P,"activity_fnorm"], pch=16, col=7, cex=1.5)
legend("topleft", c("WT","C431A","R42P"), pch=16, col=c(3,6,7), ncol=1)
quartz.save("prkn_activity_error_prop.png", type="png")


plot_elbow = function(y, signal=NA, ...) {
    print(signal)
    y = y[which(! is.na(y))]
    y = y[order(y, decreasing=T)]
    n = length(y)
    
    x1=1;      y1=y[x1]
    x2=n;      y2=y[x2]
    dx=x2-x1;  dy = y2-y1
    l = 1.0/sqrt(dy*dy+dx*dx)
    d = sapply(seq(x1,x2), function(x0){ abs(x0*dy - y[x0]*dx + x2*y1 - y2*x1) *l })
    
    # # normalized such that distances are the same in both dimensions - makes no difference?!
    # xn = seq(0,n-1)/(n-1)
    # yn = (y-y[n])/(y[1]-y[n])
    # # distance from point x0,y0 to line through x1=0,y1=1 and x2=1,y2=0
    # d = apply(cbind(xn,yn), MARGIN=1, function(v){ abs(v[1] + v[2] - 1) / sqrt(2) })
    
    xm = which.max(d)
    ym = y[xm]
    plot(y, type="l", xlab="Index", ylab="Uncertainty", ...)
    points(xm, ym, pch=16)
    text(xm, ym, sprintf("(%.0f,%.2f)",xm,ym), pos=4)
    if (is.na(signal)) {
        legend("topright", sprintf("Discard %d (%.1f%%)",xm,xm*100.0/n), pch=NA)
    } else {
        legend("topright", c(sprintf("Discard %d (%.1f%%)",xm,xm*100.0/n),sprintf("Signal to noise %.2f",signal/ym)), pch=NA)
    }
    return( c(xm,ym) )
}

quartz(width=10, height=5)
par(mfrow=c(1,2), bg="white")
elbow_ao   = plot_elbow(raw_psi$ao_sd, main="AO", signal=raw_psi[iwt,"ao_mean"]-raw_psi[iC431A,"ao_mean"], ylim=c(0,1.7))
elbow_dmso = plot_elbow(raw_psi$dmso_sd, main="DMSO", signal=raw_psi[iwt,"dmso_mean"]-raw_psi[iC431A,"dmso_mean"], ylim=c(0,1.7))
quartz.save("elbow_ao_dmso.png", type="png")

quartz(width=10, height=5)
par(mfrow=c(1,2), bg="white")
elbow_act  = plot_elbow(raw_psi$activity_sd, main="Activity", signal=raw_psi[iwt,"activity"]-raw_psi[iC431A,"activity"])
elbow_fact  = plot_elbow(raw_psi$activity_fnorm_sd, main="Activity, Keima", signal=raw_psi[iwt,"activity_fnorm"]-raw_psi[iC431A,"activity_fnorm"])
quartz.save("elbow_activity.png", type="png")

ao_sig = raw_psi[iwt,"ao_mean"] - raw_psi[iC431A,"ao_mean"]
print(sprintf("AO score S2N based on WT syn. %.1f, based on nons. %.1f, and based on elbow analysis %.1f",
              ao_sig/sd(raw_psi[isyn,"ao_mean"],na.rm=T), ao_sig/sd(raw_psi[inons,"ao_mean"],na.rm=T), 1.0/elbow_ao[2]))
print(sprintf("Activity score S2N based on WT syn. %.1f, based on nons. %.1f, and based on elbow analysis %.1f",
              1.0/sd(raw_psi[isyn,"activity"],na.rm=T), 1.0/sd(raw_psi[inons,"activity"],na.rm=T), 1.0/elbow_act[2]))
print(sprintf("Activity Keima score S2N based on WT syn. %.1f, based on nons. %.1f, and based on elbow analysis %.1f",
              1.0/sd(raw_psi[isyn,"activity_fnorm"],na.rm=T), 1.0/sd(raw_psi[inons,"activity_fnorm"],na.rm=T), 1.0/elbow_fact[2]))


##
## FACS distributions
##

load("facs/facs.rda")

breaks = seq(-1.8, 0.0, 0.02)
# breaks = seq(0,1.0,0.02)

xt_fc = facs_events[["C431A 22"]][,"ratio.acitic.neutral....PE.CF594.A.BV605.A"]
xc_fc = facs_events[["C431A D"]][,"ratio.acitic.neutral....PE.CF594.A.BV605.A"]
xt_fw = facs_events[["WT 22"]][,"ratio.acitic.neutral....PE.CF594.A.BV605.A"]
xc_fw = facs_events[["WT D"]][,"ratio.acitic.neutral....PE.CF594.A.BV605.A"]
xt_fl = facs_events[["P+K 22"]][,"ratio.acitic.neutral....PE.CF594.A.BV605.A"]
xc_fl = facs_events[["P+K D"]][,"ratio.acitic.neutral....PE.CF594.A.BV605.A"]
x_all = c(xt_fc,xc_fc,xt_fw,xc_fw,xt_fl,xc_fl)

low_cut = 1e-2
high_cut = quantile(x_all, 0.999)
breaks = seq(log10(low_cut), log10(high_cut), length.out=100)
print(sprintf("Plotting fluorescence from %.2g to %.2g discarding %d low and %d high-fluorescence events",
              low_cut, high_cut, sum(x_all<low_cut), sum(x_all>high_cut)))

xt_fc = log10( xt_fc[low_cut<xt_fc & xt_fc<high_cut] )
ht_fc = hist(xt_fc, breaks=breaks, plot=F)

xc_fc = log10( xc_fc[low_cut<xc_fc & xc_fc<high_cut] )
hc_fc = hist(xc_fc, breaks=breaks, plot=F)

xt_fw = log10( xt_fw[low_cut<xt_fw & xt_fw<high_cut] )
ht_fw = hist(xt_fw, breaks=breaks, plot=F)

xc_fw = log10( xc_fw[low_cut<xc_fw & xc_fw<high_cut] )
hc_fw = hist(xc_fw, breaks=breaks, plot=F)

xt_fl = log10( xt_fl[low_cut<xt_fl & xt_fl<high_cut] )
ht_fl = hist(xt_fl, breaks=breaks, plot=F)

xc_fl = log10( xc_fl[low_cut<xc_fl & xc_fl<high_cut] )
hc_fl = hist(xc_fl, breaks=breaks, plot=F)

# just considering FACS distributions as gaussians, what is signal-to-noise ratio
quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")
x = seq(-2,1,0.001)

plot(0,0,col=0, xlim=c(-1.5,0), ylim=c(0,6), xlab="Keima ratio, log10", ylab="Density", main="AO")
lines(ht_fl$mids, ht_fl$density, col=1)
lines(ht_fw$mids, ht_fw$density, col=5)
lines(ht_fc$mids, ht_fc$density, col=6)
lines(x, dnorm(x, mean(xt_fl), sd(xt_fl)), col=1, lty=2)
lines(x, dnorm(x, mean(xt_fw), sd(xt_fw)), col=5, lty=2)
lines(x, dnorm(x, mean(xt_fc), sd(xt_fc)), col=6, lty=2)
legend("topleft", c("Lib.","WT","C431A"), lty=1, col=c(1,5,6))
legend("topright", c("FACS","Gaus. fit"), lty=c(1,2), col=1)

plot(0,0,col=0, xlim=c(-1.5,0.0), ylim=c(0,7), xlab="Keima ratio, log10", ylab="Density", main="DMSO")
lines(hc_fl$mids, hc_fl$density, col=1)
lines(hc_fw$mids, hc_fw$density, col=5)
lines(hc_fc$mids, hc_fc$density, col=6)
lines(x, dnorm(x, mean(xc_fl), sd(xc_fl)), col=1, lty=2)
lines(x, dnorm(x, mean(xc_fw), sd(xc_fw)), col=5, lty=2)
lines(x, dnorm(x, mean(xc_fc), sd(xc_fc)), col=6, lty=2)
legend("topleft", c("Lib.","WT","C431A"), lty=1, col=c(1,5,6))
legend("topright", c("FACS","Gaus. fit"), lty=c(1,2), col=1)

plot(0,0,col=0, xlim=c(-0.4,0.9), ylim=c(0,3.5), xlab="delta-Keima ratio, log10", ylab="Density", main="AO-DMSO")
lines(x, dnorm(x, mean(xt_fl)-mean(xc_fl), sqrt(var(xt_fl)+var(xc_fl))), col=1)
delta_wt = mean(xt_fw)-mean(xc_fw)
delta_wt_sd = sqrt(var(xt_fw)+var(xc_fw))
delta_c431a = mean(xt_fc)-mean(xc_fc)
lines(x, dnorm(x, delta_wt, delta_wt_sd), col=5)
lines(x, dnorm(x, delta_c431a, sqrt(var(xt_fc)+var(xc_fc))), col=6)
legend("topleft", c("Library","WT","C421A"), lty=1, col=c(1,5,6))
legend("topright", sprintf("Singal-to-noise %.2f", (delta_wt-delta_c431a)/delta_wt_sd), pch=NA)

quartz.save("facs_gauss.png", type="png")


quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")
plot_distr(pf$ao_mean, xlab="Average gate", main="AO")
plot_distr(pf$dmso_mean, xlab="Average gate", , main="DMSO")
plot_distr(pf$activity, xlab="Activity score", main="Activity")
quartz.save("average_gate_distributions.png", type="png")


quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")

plot_distr(pf$ao_fnorm, xlab="Keima ratio, log10", xlim=c(-1.2,0), ylim=c(0,8), main="AO Keima transform")
lines(ht_fw$mids, ht_fw$density, col=5)
lines(ht_fc$mids, ht_fc$density, col=6)
legend("topright", c("FACS WT","FACS C431A"), lty=c(1,1), col=c(5,6))

plot_distr(pf$dmso_fnorm, xlab="Keima ratio, log10", xlim=c(-1.2,0), ylim=c(0,8), main="DMSO Keima transform")
lines(hc_fw$mids, hc_fw$density, col=5)
lines(hc_fc$mids, hc_fc$density, col=6)
legend("topright", c("FACS WT","FACS C431A"), lty=c(1,1), col=c(5,6))

delta_fnorm = pf$ao_fnorm - pf$dmso_fnorm
plot_distr(delta_fnorm, xlab="delta-of-log10 Keima ratio", xlim=c(-0.6,1.0), ylim=c(0,6), main="AO-DMSO Keima transform")
x = seq(-1,1,0.001)
lines(x, dnorm(x, mean(xt_fw)-mean(xc_fw), sqrt(var(xt_fw)+var(xc_fw))), col=5)
lines(x, dnorm(x, mean(xt_fc)-mean(xc_fc), sqrt(var(xt_fc)+var(xc_fc))), col=6)
legend("topright", c("FACS WT","FACS C431A"), lty=c(1,1), col=c(5,6))

quartz.save("facs_cyto_with_transforms.png", type="png")


quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")

plot_distr(pf$ao_fwnorm, xlab="Keima ratio, log10", xlim=c(-1.2,0), ylim=c(0,8), main="AO Keima transform, weighted")
lines(ht_fw$mids, ht_fw$density, col=5)
lines(ht_fc$mids, ht_fc$density, col=6)
legend("topright", c("FACS WT","FACS C431A"), lty=c(1,1), col=c(5,6))

plot_distr(pf$dmso_fwnorm, xlab="Keima ratio, log10", xlim=c(-1.2,0), ylim=c(0,8), main="DMSO Keima transform, weighted")
lines(hc_fw$mids, hc_fw$density, col=5)
lines(hc_fc$mids, hc_fc$density, col=6)
legend("topright", c("FACS WT","FACS C431A"), lty=c(1,1), col=c(5,6))

delta_fwnorm = pf$ao_fwnorm - pf$dmso_fwnorm
plot_distr(delta_fwnorm, xlab="delta-of-log10 Keima ratio", xlim=c(-0.6,1.0), ylim=c(0,6), main="AO-DMSO Keima transform, weighted")
x = seq(-1,1,0.001)
lines(x, dnorm(x, mean(xt_fw)-mean(xc_fw), sqrt(var(xt_fw)+var(xc_fw))), col=5)
lines(x, dnorm(x, mean(xt_fc)-mean(xc_fc), sqrt(var(xt_fc)+var(xc_fc))), col=6)
legend("topright", c("FACS WT","FACS C431A"), lty=c(1,1), col=c(5,6))

quartz.save("facs_cyto_with_transforms_weighted.png", type="png")


quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")

plot_distr(pf$activity, xlab="Activity score", xlim=c(-1.0,1.8),  main="AO-1.77*DMSO avg. gate")

plot_distr(delta_fnorm, xlab="Activity score", xlim=c(-0.6,1.0), ylim=c(0,6), main="log10(AO)-log10(DMSO) Keima transform")
x = seq(-1,1,0.001)
lines(x, dnorm(x, mean(xt_fw)-mean(xc_fw), sqrt(var(xt_fw)+var(xc_fw))), col=5)
lines(x, dnorm(x, mean(xt_fc)-mean(xc_fc), sqrt(var(xt_fc)+var(xc_fc))), col=6)
legend("topright", c("FACS WT","FACS C431A"), lty=c(1,1), col=c(5,6))

plot_distr(delta_fwnorm, xlab="Activity score", xlim=c(-0.6,1.0), ylim=c(0,6), main="log10(AO)-log10(DMSO) Keima transform, weighted")
x = seq(-1,1,0.001)
lines(x, dnorm(x, mean(xt_fw)-mean(xc_fw), sqrt(var(xt_fw)+var(xc_fw))), col=5)
lines(x, dnorm(x, mean(xt_fc)-mean(xc_fc), sqrt(var(xt_fc)+var(xc_fc))), col=6)
legend("topright", c("FACS WT","FACS C431A"), lty=c(1,1), col=c(5,6))

quartz.save("facs_cyto_with_activity.png", type="png")


##
## Average gate replica correlations
##

plot_cor = function(df, pch=16, cex=.2, repnames=NA, n_tot_var=NA, ...) {
    text_rp = function(x, y, cex=1, ...) {
        rp = cor(x,y,method="pearson",use="complete.obs")
        x_mid = (max(x, na.rm=T) + min(x, na.rm=T)) /2.0
        y_mid = (max(y, na.rm=T) + min(y, na.rm=T)) /2.0
        text(x_mid, y_mid, sprintf("Pearson\n%.2f",rp), cex=2, ...)
    }
    if (any(is.na(repnames))) { repnames = colnames(df) } else { stopifnot(length(repnames)==ncol(df)) }
    rps = cor(df, method="pearson", use="complete.obs")
    rps = rps[upper.tri(rps)]
    if (ncol(df) == 2) {
        plot(df, xlab=repnames[1], ylab=repnames[2], pch=pch, cex=cex, ...)
    } else if (ncol(df) > 2) {
        sums = apply(df, MARGIN=2, sum)
	if (is.na(n_tot_var)) { n_tot_var = nrow(df) }
        # covers = apply(df, MARGIN=2, function(v) { sum(! is.na(v))*100.0/n_tot_var })
        # labels = sprintf("%s\n%.2f%% cover",repnames,covers)
        covers = apply(df, MARGIN=2, function(v) { sum(! is.na(v)) })
        labels = sprintf("%s\n%d scores",repnames,covers)
        plot(df, labels=labels, upper.panel=text_rp, pch=pch, cex=cex, ...)
	title(xlab="Average gate", ylab="Average gate")
    }
}

pdf("rep_cor_ao.pdf", width=7, height=7, pointsize=9)
plot_cor(raw_psi[,c("ao_1_psi","ao_2_psi","ao_3_psi")], repnames=c("AO rep. 1","AO rep. 2","AO rep. 3"))
dev.off()

pdf("rep_cor_dmso.pdf", width=7, height=7, pointsize=9)
plot_cor(raw_psi[,c("dmso_1_psi","dmso_2_psi","dmso_3_psi")], repnames=c("DMSO rep. 1","DMSO rep. 2","DMSO rep. 3"))
dev.off()

