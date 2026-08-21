options(width=160)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    infiles = c("Exp_5p5_C431A_AO.csv")
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript facs.r  <facs1.csv>  [facs2.csv  ...]")
    quit(save="no")
} else {
    # infiles = args[2:length(args)]
    infiles = args
}
# print(sprintf("Args: %s",paste0(args, collapse=" ")))


parse_facs_csv = function(filename, n_cols=20) {
    d = read.table(filename, sep=",", row.names=NULL, fill=T, col.names=sprintf("col%.02d",seq(n_cols)))

    # keep reading until an excess of columns are read (read.table(fill=T)) only looks at first 5 lines)
    while (! all(is.na(d[,n_cols]) | n_cols > 200) ) {
        n_cols = n_cols+10
        d = read.table(filename, sep=",", row.names=NULL, fill=T, col.names=sprintf("col%.02d",seq(n_cols)))
    }

    facs = list()
    
    # find tube name
    iname = which(d$col01 == "TUBE NAME")
    if (length(iname) < 1) { iname = which(grepl("TUBE", d$col01)) }
    if (length(iname) > 1) {
        print(sprintf("WARNING: Found %d TUBE NAME labels (using last): %s", length(iname), paste0(d[iname,"col02"],collapse=", ")))
	iname = iname[length(iname)]
    }
    if (length(iname) < 0) {
        print("WARNING: No TUBE NAME label found")
        facs$sample = "no_name"
    } else {
        facs$sample = d[iname,"col02"]
    }

    # find date
    idate = which(d$col01 == "$DATE")
    if (length(idate) < 1) { idate = which(grepl("DATE", d$col01)) }
    if (length(idate) > 1) {
        print(sprintf("WARNING: Found %d DATE labels (using last): %s", length(idate), paste0(d[idate,"col02"],collapse=", ")))
        idate = idate[length(idate)]
    }
    if (length(idate) < 0) {
        print("WARNING: No DATE label found")
      	facs$date = "no_date"
    } else {
        facs$date = d[idate,"col02"]
    }   
    

    # find the number of rows with settings assuming a header row. Look for the row from which everything in column 1 is numeric
    # print("Search end of settings lines")
    settings_rows = 0
    while (suppressWarnings( any(is.na(as.numeric(d[(settings_rows+2):nrow(d),"col01"]))) & settings_rows < 500)) {
        settings_rows = settings_rows +1
    }
    stopifnot(settings_rows < 500)

    # # print("Search number of settings cloumns")
    # settings_cols = n_cols
    # while (all(is.na(d[1:settings_rows,settings_cols])) & settings_cols > 0) {
    #     settings_cols = settings_cols -1
    # }
    # stopifnot(settings_cols > 0)
    # settings = d[1:settings_rows,1:settings_cols]
    
    facs$events = read.table(filename, sep=",", skip=settings_rows, header=T)
    return(facs)
}


##
## Settings
##
facs_set = list()
facs_set$infiles = infiles

facs_events = list()
facs_table = list()
for (infile in infiles) {
    print(sprintf("Reading %s",infile))
    facsfile = parse_facs_csv(infile)
    print(sprintf("    Found sample %s from %s with %d events and %d channels",facsfile$sample,facsfile$date,nrow(facsfile$events),ncol(facsfile$events)))
    name = facsfile$sample
    if ( name %in% names(facs_events)) {
        i = 1
	name = sprintf("%s_%02d", facsfile$sample, i)
	while (name %in% names(facs_events)) {
	    i = i+1
	    name = sprintf("%s_%02d", facsfile$sample, i)
	}
	print(sprintf("    Renaming %s to %s", facsfile$sample, name))
    }
    facs_table[[name]] = c(facsfile$sample, facsfile$date, nrow(facsfile$events), ncol(facsfile$events))
    facs_events[[name]] = facsfile$events
}

# channels we are using
# facs_set$cn_gfp = "GFP.A"
# facs_set$cn_cherry = "PE.Texas.Red.A"
# facs_set$cn_ratio = "Derived....GFP.A.PE.Texas.Red.A"

facs_set$cn_keima_acid = "PE.CF594.A"  # ex. 568nm;  em. 620nm
facs_set$cn_keima_neut = "BV605.A"     # ex. 440nm;  em. 620nm
facs_set$cn_ratio = "ratio.acitic.neutral....PE.CF594.A.BV605.A"

# only store used channels
cns = c(facs_set$cn_keima_acid, facs_set$cn_keima_neut, facs_set$cn_ratio)
for (sample in names(facs_events)) {
    stopifnot(all( cns %in% colnames(facs_events[[sample]]) ))
    facs_events[[sample]] = facs_events[[sample]][,cns]
}

# reformat as table
facs_table = data.frame(name = names(facs_table),
                        sample = sapply(facs_table, "[[", 1),
                        date = sapply(facs_table, "[[", 2),
                        events = sapply(facs_table, "[[", 3))
write.csv(facs_table, quote=F, row.names=F, file="facs_table.csv")


# plot all distributions
hl = list()
x_min = 0; x_max = 0; y_max = 0
for (sample in names(facs_events)) {
    h = hist(facs_events[[sample]][,facs_set$cn_ratio], breaks=200, plot=F)
    x_min = min(c(x_min, h$mids))
    x_max = max(c(x_max, h$mids))
    y_max = max(c(y_max, h$density))
    hl[[sample]] = h
}

quartz(width=15, height=6)
par(mfrow=c(1,2), bg="white")

# plot(0,0,col=0, xlim=c(x_min,x_max), ylim=c(0,y_max), xlab="GFP/mCherry", ylab="Density")
plot(0,0,col=0, xlim=c(-0.1,0.7), ylim=c(0,y_max), xlab="Keima acedic/neutral", ylab="Density")
for (i in seq_along(facs_events)) {
    sample = names(facs_events)[i]
    lines(hl[[sample]]$mids, hl[[sample]]$density, col=i, lty=i%/%8+1)
}
legend("topright", names(hl), col=seq(i), lty=seq(i)%/%8+1, ncol=3, cex=.8)

plot(0,0,col=0, xlim=c(-1.5,0.1), ylim=c(0,6), xlab="Keima acedic/neutral, log10", ylab="Density")
for (i in seq_along(facs_events)) {
    sample = names(facs_events)[i]
    ipos = which(facs_events[[sample]][,facs_set$cn_ratio] > 0.0)
    if (length(ipos) < nrow(facs_events[[sample]])) {
        print(sprintf("Not plotting %d events with non-positive fluorescense for %s",nrow(facs_events[[sample]])-length(ipos),sample))
    }
    h = hist(log10(facs_events[[sample]][ipos,facs_set$cn_ratio]), breaks=100, plot=F)
    lines(h$mids, h$density, col=i, lty=i%/%8+1)
}

quartz.save("all_distributions.png", type="png")


# # I guess these are the samples Vasileios used for plotting, see photo in slack DM from Aug 13, 2024
# facs_set$sel = list()
# facs_set$sel[["even1"]] = "E1"      # File 26-10-2022_E1_007.fcs    with 4.32 E5 events - alt E1b, E1b_02
# facs_set$sel[["even2"]] = "E2"      # File 07-12-2022_E2_007.fcs    with 4.51 E5 events - alt E2b_01, E2b_02
# facs_set$sel[["odd3"]]  = "O3b_01"  # File 26-10-2022_O3b_006.fcs   with 395055  events - alt O3b, O3_01
# facs_set$sel[["odd4"]]  = "O4"      # File 27-10-2022_O4_005.fcs    with 4.15 E5 events - alt O4b, O4_02, O4b_02
# facs_set$sel[["ct1"]]   = "CT1 I"   # File 16-12-2022_CT1 I_005.fcs with 4.67 E5 events - alt CT1 IIb

save(facs_events, facs_table, facs_set, file="facs.rda")

# # plot selected distributions
# quartz(width=8, height=6)
# plot(0,0,col=0, xlim=c(-0.1,1.1), ylim=c(0,6), xlab="Keims acedic/neutral", ylab="Density")
# i=0
# for (lib in names(facs_set$sel)) {
#     i=i+1
#     sample = facs_set$sel[[lib]]
#     lines(hl[[sample]]$mids, hl[[sample]]$density, col=i)
# }
# legend("topright", unlist(facs_set$sel), lty=1, col=seq(i))
# quartz.save("selected_distributions.png", type="png")
