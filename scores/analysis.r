options(width=160, digits=4)

# Passing Bablok and Deming regression (loss symm. in y and x)
library("mcr")

aa_one = strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

##
## Data
##
pf = read.csv("prkn_activity_full.csv")
colnames(pf) = c("var",colnames(pf)[2:ncol(pf)])

nc = nchar(pf$var)
pf$aa = substr(pf$var, 1, 1)
pf$resi = as.numeric(substr(pf$var, 2, nc-1))
pf$mut = substr(pf$var, nc, nc)

iwt = which(pf$var == "WT")
pf[iwt,"aa"] = ""
pf[iwt,"resi"] = 0
pf[iwt,"mut"] = ""

print("Score and std summary")
print(summary(pf$activity))
print(summary(pf$activity_sd))

isyn = which(pf$mut == "=")
inons = which(pf$mut == "*")
i1 = which(pf$mut %in% aa_one)
print(sprintf("Data has %d single var., %d nonsense, %d synonymous and %d WT", length(i1), length(inons), length(isyn), length(iwt)))

isyn_meas = which(pf$mut == "=" & ! is.na(pf$activity))
inons_meas = which(pf$mut == "*" & ! is.na(pf$activity))
i1_meas = which(pf$mut %in% aa_one & ! is.na(pf$activity))
print(sprintf("Without NA, data has %d single var., %d nonsense and %d synonymous", length(i1_meas), length(inons_meas), length(isyn_meas)))

cut = 0.614
isyn_filt = which(pf$mut == "=" & pf$activity_sd < cut)
inons_filt = which(pf$mut == "*" & pf$activity_sd < cut)
i1_filt = which(pf$mut %in% aa_one & pf$activity_sd < cut)
print(sprintf("After filtering std<%.3f, data has %d single var., %d nonsense and %d synonymous",
              cut, length(i1_filt), length(inons_filt), length(isyn_filt)))

# Control variant, low abundance, residual function when overexpressed
iR42P = which(pf$var == "R42P")

# Control catalytic site variant, high abundance, no activity
iC431A = which(pf$var == "C431A")


##
## DMSO scaling - supplementary Fig 6
##

plot_ao_dmso = function(x, y, dmso_fac, ...) {
    rp = cor(x[inons], y[inons], method="pearson", use="complete.obs")
    print(sprintf(sprintf("Nonsense Pearson %.2f",rp)))
    
    ymin = min(y,na.rm=T); ymax = max(y,na.rm=T); yd = ymax-ymin
    plot(x, y, pch=16, cex=.3, ylim=c(ymin-yd*0.28,ymax+yd*0.28), ...)
    points(x[inons],  y[inons],  pch=16, cex=.5,  col=2)
    points(x[isyn],   y[isyn],   pch=16, cex=.5,  col=4)
    points(x[iwt],    y[iwt],    pch=16, cex=1.5, col=3)
    points(x[iR42P],  y[iR42P],  pch=16, cex=1.5, col=7)
    points(x[iC431A], y[iC431A], pch=16, cex=1.5, col=6)
    
    dmso_inct = mean(y[inons]-x[inons]/dmso_fac, na.rm=T)
    abline(c(dmso_inct,1.0/dmso_fac), lty=1)
    
    fit_dam = mcreg(x[inons], y[inons], method.reg="Deming")
    abline(coef(fit_dam)[1:2], lty=2)
    s_dam = sprintf("Daming DMSO fac. %.2f", 1/coef(fit_dam)[2])
    
    legend("topleft", c("Single var.","WT syn.","Nons.","WT","C431A","R42P"), pch=16, col=c(1,4,2,3,6,7), pt.cex=c(.5,.8,.8,1.3,1.3,1.3), ncol=3)
    legend("bottomright", c(sprintf("DMSO factor %.2f",dmso_fac),s_dam), lty=c(1,2))
}    

# quartz(width=10, height=5)
# par(mfrow=c(1,2), mar=c(5,4,2,2)+.1, bg="white")
# plot_ao_dmso(pf$ao_mean,  pf$dmso_mean,  1.77, xlab="AO average bin", ylab="DMSO average bin")
# plot_ao_dmso(pf$ao_fnorm, pf$dmso_fnorm, 1.00, xlab="AO log10(mt-mKeima ratio)", ylab="DMSO log10(mt-mKeima ratio)")
rp_ab = cor(pf[inons,"ao_mean"],  pf[inons,"dmso_mean"], method="pearson", use="complete.obs")
rp_fn = cor(pf[inons,"ao_fnorm"],  pf[inons,"dmso_fnorm"], method="pearson", use="complete.obs")
print(sprintf("Pearson avg.bin %.2f and FACS-norm %.2f",rp_ab,rp_fn))

plot_distrib = function(x, xlim=NA, ylim=NA, ylab=NA, ...) {
    xmin = min(x,na.rm=T); xmax = max(x,na.rm=T)
    rf = 10
    breaks = seq(floor(xmin*rf)/rf, ceiling(xmax*rf)/rf, length.out=50)
    h1 = hist(x[i1],    breaks=breaks, plot=F)
    hs = hist(x[isyn],  breaks=breaks, plot=F)
    hn = hist(x[inons], breaks=breaks, plot=F)

    ylab = ifelse(is.na(ylab), "Density", ylab)
    if (any(is.na(xlim))) { xlim = c(xmin,xmax) }
    if (any(is.na(ylim))) { ylim = c(0,max(c(h1$density,hs$density,hn$density))) }
    
    plot(0,0,col=0, xlim=xlim, ylim=ylim, ylab=ylab, ...)
    lines(hn$mids, hn$density, col=2, lwd=1.5)
    lines(hs$mids, hs$density, col=4, lwd=1.5)
    lines(h1$mids, h1$density, col=1, lwd=1.5)
    points(x[iwt],    0, pch=16, cex=1.5, col=3)
    points(x[iR42P],  0, pch=16, cex=1.5, col=7)
    points(x[iC431A], 0, pch=16, cex=1.5, col=6)
    legend("topleft", c("Single var.","WT syn.","Nonsense","WT","C431A","R42P"),
           pch=c(NA,NA,NA,16,16,16), lty=c(1,1,1,NA,NA,NA), col=c(1,4,2,3,6,7), pt.cex=1.3, ncol=1)
}

# quartz(width=10, height=5)
# par(mfrow=c(1,2), mar=c(5,4,2,2)+.1, bg="white")
# plot_distrib(pf$activity, xlab="Activity score", xlim=c(-1.4,2))
# plot_distrib(pf$activity_fnorm, xlab="Activity log10(mt-mKeima ratio)")


ltp = read.csv("fig2c_ltp_validation.csv")
i_ltp2htp = match(pf$var,ltp$var)
pf$ltp    = ltp[i_ltp2htp,"Normalized_mean_LT"]
pf$ltp_sd = ltp[i_ltp2htp,"Normalized_sd_LT"]
# this should not be zero
pf[iwt,"ltp_sd"] = 1e-3

plot_ltp = function(x_all, dx_all, y_all, dy_all, mark_vars=c("WT","C431A","R42P"), mark_pos=c(2,2,2), xlim=NA, ylim=NA, ylab=NA, ...) {
    iltp = which(! is.na(pf$ltp))
    x = x_all[iltp]
    dx = dx_all[iltp] +0.001
    # y = pf[iltp,"ltp"]
    # dy = pf[iltp,"ltp_sd"] +0.001
    y = y_all[iltp]
    dy = dy_all[iltp] +0.001

    ylab = ifelse(is.na(ylab), "LTP mt-mKeima ratio", ylab)
    if (any(is.na(xlim))) { xlim = c(min(x-dx,na.rm=T), max(x+dx,na.rm=T)) }
    if (any(is.na(ylim))) { ylim = c(min(y-dy,na.rm=T), max(y+dy,na.rm=T)) }
    plot(0,0,col=0, xlim=xlim, ylim=ylim, ylab=ylab, pch=20, ...)
    arrows(x0=x-dx, y0=y, x1=x+dx, code=3, angle=90, length=.02, col="grey70")
    arrows(x0=x, y0=y-dy, y1=y+dy, code=3, angle=90, length=.02, col="grey70")
    points(x, y, pch=16, col=1)

    stopifnot(length(mark_pos)==length(mark_vars))
    names(mark_pos) = mark_vars
    for (mark_var in mark_vars) {
        ii = which(pf[iltp,"var"] == mark_var)
	stopifnot(length(ii) == 1)
        text(x[ii], y[ii], mark_var, pos=mark_pos[mark_var])
    }

    rp = cor(x, y, method="pearson", use="complete.obs")
    rs = cor(x, y, method="spearman", use="complete.obs")
    legend("topleft", c(sprintf("Pearson %.2f",rp), sprintf("Spearman %.2f",rs)), pch=NA)
    # legend("topleft", c(sprintf("Pearson %.2f",rp)), pch=NA)
}


quartz(width=10, height=10)
par(mfrow=c(2,2), mar=c(5,4,2,2)+.1, bg="white")
plot_ltp(pf$activity, pf$activity_sd, pf$ltp, pf$ltp_sd, xlab="HTP Activity score", mark_vars=c("WT","R42P"), mark_pos=c(2,2))
plot_ltp(pf$activity_fnorm, pf$activity_fnorm_sd, pf$ltp, pf$ltp_sd, xlab="HTP log10(mt-mKeima ratio)", mark_vars=c("WT","R42P","S246*"), mark_pos=c(2,2,2))

ltp_log = log10(pf$ltp)
ltp_log_sd = pf$ltp_sd / (pf$ltp*log(10))
plot_ltp(pf$activity_fnorm, pf$activity_fnorm_sd, ltp_log, ltp_log_sd, xlab="HTP log10(mt-mKeima ratio)", ylab="LTP log10(mt-mKeima ratio)", mark_vars=c("WT","R42P","S246*"), mark_pos=c(2,2,2))

act_fn_nonlog = 10^pf$activity_fnorm
act_fn_nonlog_sd = pf$activity_fnorm_sd * log(10) * act_fn_nonlog
plot_ltp(act_fn_nonlog, act_fn_nonlog_sd, pf$ltp, pf$ltp_sd, xlab="HTP mt-mKeima ratio", mark_vars=c("WT","R42P","S246*"), mark_pos=c(2,2,2))



pdf("dmso_scaling.pdf", width=7, height=4.66, pointsize=9)
par(mfrow=c(2,3), mar=c(5,4,2,2)+.1, bg="white")

plot_ao_dmso(pf$ao_mean,  pf$dmso_mean,  1.77, xlab="AO average bin", ylab="DMSO average bin")
text(0.8, 4.4, "A", xpd=T, font=2, cex=2)

plot_distrib(pf$activity, xlab="Activity score", xlim=c(-1.5,1.9))
text(-2.4, 3.01, "B", xpd=T, font=2, cex=2)

plot_ltp(pf$activity, pf$activity_sd, pf$ltp, pf$ltp_sd, xlab="HTP Activity score", mark_vars=c("WT","R42P"), mark_pos=c(2,2))
text(-0.95, 1.2, "C", xpd=T, font=2, cex=2)

plot_ao_dmso(pf$ao_fnorm, pf$dmso_fnorm, 1.00, xlab="AO log10(mt-mKeima ratio)", ylab="DMSO log10(mt-mKeima ratio)")
text(-1.61, 1.15, "D", xpd=T, font=2, cex=2)

plot_distrib(pf$activity_fnorm, xlab="Activity log10(mt-mKeima ratio)", xlim=c(-2.6,2.8))
text(-3.9, 1.74, "E", xpd=T, font=2, cex=2)

# plot_ltp(pf$activity_fnorm, pf$activity_fnorm_sd, pf$ltp, pf$ltp_sd, xlab="HTP log10(mt-mKeima ratio)", mark_vars=c("WT","R42P","S246*"), mark_pos=c(2,2,2))
plot_ltp(pf$activity_fnorm, pf$activity_fnorm_sd, ltp_log, ltp_log_sd, xlab="HTP log10(mt-mKeima ratio)", ylab="LTP log10(mt-mKeima ratio)", mark_vars=c("WT","R42P","S246*","C431A"), mark_pos=c(2,2,3,1))
text(-3.0, 0.15, "F", xpd=T, font=2, cex=2)

dev.off()



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
	title(xlab="Average bin", ylab="Average bin")
    }
}

pdf("rep_cor_ao.pdf", width=7, height=7, pointsize=9)
plot_cor(pf[,c("ao_1_psi","ao_2_psi","ao_3_psi")], repnames=c("AO rep. 1","AO rep. 2","AO rep. 3"))
dev.off()

pdf("rep_cor_dmso.pdf", width=7, height=7, pointsize=9)
plot_cor(pf[,c("dmso_1_psi","dmso_2_psi","dmso_3_psi")], repnames=c("DMSO rep. 1","DMSO rep. 2","DMSO rep. 3"))
dev.off()


##
## Per residue
##

pfr = read.csv("prkn_activity_residues.csv")
pf$act_med = pfr[match(pf$resi,pfr$resi),"med_score"]

# make sure to use the same data in the following
ii = which(! is.na(pf[i1,"activity"]))
i = i1[ii]

# The explained variance is R^2 and for a linear model with intersept this is rp^2
rp = cor(pf[i,"activity"], pf[i,"act_med"], method="pearson")

# a model of variant score based on position medians
fit = lm(pf[i,"activity"] ~ pf[i,"act_med"])
act_pred = coef(fit)[1] + coef(fit)[2] * pf[i,"act_med"]
# Fraction of unexplained variance is residual variance divided by original variance
fuv = var(pf[i,"activity"]-act_pred)/var(pf[i,"activity"])

print(sprintf("Fraction of unexplained variance %.1f%% (ratio of residual to total var) and explained %.1f%% (Pearson squared). Sum %.3f%%",
              fuv*100, rp*rp*100, fuv*100+rp*rp*100))

