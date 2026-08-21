options(width=160)

aa_one = strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
t2o = unlist(list("Ala"="A", "Cys"="C", "Asp"="D", "Glu"="E", "Phe"="F", "Gly"="G", "His"="H", "Ile"="I", "Lys"="K", "Leu"="L",
                  "Met"="M", "Asn"="N", "Pro"="P", "Gln"="Q", "Arg"="R", "Ser"="S", "Thr"="T", "Val"="V", "Trp"="W", "Tyr"="Y",
                  "Del"="del", "="="=", "Ter"="*"))

# Per-residue functional data
pfr = read.csv("../scores/prkn_activity_residues_251013.csv")

all_var = sprintf("%s%d%s", rep(pfr$wt,each=22), rep(pfr$resi,each=22), rep(c(aa_one,'*','='),nrow(pfr)))
prkn = data.frame(var=c("WT",all_var))
nc = nchar(prkn$var)
prkn$aa = substr(prkn$var, 1, 1)
prkn$resi = as.numeric(substr(prkn$var, 2, nc-1))
prkn$mut = substr(prkn$var, nc, nc)
prkn = prkn[-which(prkn$aa==prkn$mut),]

iwt = which(prkn$var == "WT")
prkn[iwt,"aa"] = ""
prkn[iwt,"resi"] = 0
prkn[iwt,"mut"] = ""

isyn = which(prkn$mut == "=")
inons = which(prkn$mut == "*")
i1 = which(prkn$mut %in% aa_one)

iR42P = which(prkn$var == "R42P")
iC431A = which(prkn$var == "C431A")

# Released activity scores
pf = read.csv("../scores/prkn_activity_251013.csv")
i = match(prkn$var,pf$var)
cns = colnames(pf)[2:ncol(pf)]
prkn[,cns] = pf[i,cns]

# abundance data with nonsense and the new categories, see abundance.r
pa = read.csv("prkn_abundance.csv")

i = match(prkn$var,pa$var)
prkn$abun = pa[i,"vamp_score"]
prkn$abun_sd = pa[i,"vamp_std"]

# annotations from clausen24
lene = read.csv("clausen24_source_data_s01.csv", na.strings=c("NaN","NA"))
i = match(prkn$var,lene$variant)

# check that missense abundance scores are the same
stopifnot(all( prkn[i1,"abun"] - lene[i[i1],"abundance_score"] < 1e-3, na.rm=T ))
stopifnot(all( prkn[i1,"abun_sd"] - lene[i[i1],"abundance_score_std"] < 1e-3, na.rm=T ))

# use GEMME, Rosetta, pLDDT, rASA and GnomeAD from parkin abundance paper
prkn$gemme = lene[i,"gemme_de"]
prkn$rosetta = lene[i,"normalized_rosetta_ddg"]
prkn$plddt = lene[i,"AF_confidence"]
prkn$rasa = lene[i,"rASA"]
prkn$gnomad21 = lene[i,"gnomAD_allele_frequency"]
prkn[which(prkn$mut=="="),c("plddt","rasa")] = NA

# annotations from cagiada & jonsson 2024
nico = read.csv("cagiada24_O60260-PRKN.csv")
i = match(prkn$var,nico$variant)

prkn$esm1b = nico[i,"esm_1b_masked_marginals"]
prkn$esmif = nico[i,"esm_if_masked_marginals"]
prkn$clinvar = nico[i,"clinvar_clinical_significance"]
# prkn[which(is.na(prkn$clinvar)),"clinvar"] = ""
# where clinvar is not available use NA consistantly (instead of empty str)
prkn[which(prkn$clinvar==""),"clinvar"] = NA
prkn[which(prkn$mut=="="),c("esm1b","esmif")] = NA

print("New clinvar (all)")
print(table(prkn$clinvar))
print("New clinvar (with abundance scores)")
print(table(prkn[which(! is.na(prkn$abun)),"clinvar"]))

print("Old clinvar (all)")
print(table(lene$disease))
print("Old clinvar (with abundance scores)")
print(table(lene[which(! is.na(lene$abundance_score)),"disease"]))


##
## gnomad 4.1
##
gnomad = read.csv("gnomAD_v4.1.0_ENSG00000185345_2025_09_16_10_26_39.csv")
ig_use = which(gnomad$VEP.Annotation=="missense_variant")
print(sprintf("Loaded gnomAD 4.1 with %d (%.1f%%) missense variants", length(ig_use), length(ig_use)/nrow(gnomad)*100))
print(sprintf("Of these are %d non-unique variants (at protein level)", length(gnomad[ig_use,"Protein.Consequence"]) - length(unique(gnomad[ig_use,"Protein.Consequence"]))))

nci = nchar(gnomad[ig_use,"Protein.Consequence"])
gnomad$aa3_wt = NA
gnomad[ig_use,"aa3_wt"] = substr(gnomad[ig_use,"Protein.Consequence"], 3, 5)
gnomad$aa3_mut = NA
gnomad[ig_use,"aa3_mut"] = substr(gnomad[ig_use,"Protein.Consequence"], nci-2, nci)
gnomad$resi = NA
gnomad[ig_use,"resi"] = as.numeric( substr(gnomad[ig_use,"Protein.Consequence"], 6, nci-3) )
gnomad$aa1_wt = t2o[gnomad$aa3_wt]
gnomad$aa1_mut = t2o[gnomad$aa3_mut]
gnomad$var = sprintf("%s%d%s", gnomad$aa1_wt, gnomad$resi, gnomad$aa1_mut)
# print(gnomad[ig_use,c("Protein.Consequence","var","aa1_wt","resi","aa1_mut","Allele.Count","Allele.Frequency","Filters...joint","ClinVar.Germline.Classification")])

used_transcripts = c()
tt = table(gnomad$Transcript)
print(sprintf("gnomad data has %d transcripts", length(tt)))
print(table(gnomad[ig_use,"Transcript"]))

for (transcript in names(tt)) {
    igt = ig_use[which(gnomad[ig_use,"Transcript"] == transcript)]
    wt_is_match = gnomad[igt,"aa1_wt"] == pfr$wt[gnomad[igt,"resi"]]
    if (all(wt_is_match)) {
        print(sprintf("    all %d varinats match for transcript %s", length(igt),transcript))
	used_transcripts = c(transcript,used_transcripts)
    } else {
        print(sprintf("    of %d varinats, only %d match transcript %s", length(igt), sum(wt_is_match), transcript))        
    }
}

iig = which(gnomad[ig_use,"Transcript"] %in% used_transcripts)
ts_str = paste0(used_transcripts,collapse=",")
print(sprintf("Use %d of %d (%.1f%%) missense variants from transcript(s): %s", length(iig), length(ig_use), length(iig)/length(ig_use)*100, ts_str))
ig_use = ig_use[iig]

# The filter value AC0 could mean "allele count zero" but all 3 has genomic filter PASS and one allele count? Not sure if I should remove these
print(table(gnomad[ig_use,c("Source","Filters...joint")]))

stopifnot( length(gnomad[ig_use,"Transcript.Consequence"]) == length(unique(gnomad[ig_use,"Transcript.Consequence"])) )
vt = table(gnomad[ig_use,"Protein.Consequence"])
ii = which(gnomad[ig_use,"Protein.Consequence"] %in% names(vt[vt==1]))
ndf = gnomad[ig_use[ii],c("var","Protein.Consequence","Transcript.Consequence","Allele.Count","Allele.Number","Allele.Frequency")]

for (redun_aa_var in names(vt[vt>1])) {
    ii = which(gnomad[ig_use,"Protein.Consequence"] == redun_aa_var)
    al_count = sum( gnomad[ig_use[ii],"Allele.Count"] )
    # al_num = mean( gnomad[ig_use[ii],"Allele.Number"] )
    ii_max = ii[ which.max(gnomad[ig_use[ii],"Allele.Count"]) ]
    al_num = gnomad[ig_use[ii_max],"Allele.Number"]
    new_af = al_count/al_num
    cut = 1e-3
    # if (any( abs(al_num-gnomad[ig_use[ii],"Allele.Number"])/al_num > cut )) {
    #     print(sprintf("WARNING: Allele Numbers varries >%.3f%% for %s using %.1f and new AF %.5g", cut*100, redun_aa_var, al_num, new_af))
    # 	print(gnomad[ig_use[ii],c("Protein.Consequence","Transcript.Consequence","Allele.Count","Allele.Number","Allele.Frequency")])
    # }
    print(sprintf("Merged %d alleles for %s using AN %.1f and new AF %.5g", length(ii), redun_aa_var, al_num, new_af))
    print(gnomad[ig_use[ii],c("Protein.Consequence","Transcript.Consequence","Allele.Count","Homozygote.Count","Allele.Number","Allele.Frequency")])
    ts_str = paste0(gnomad[ig_use[ii],"Transcript.Consequence"],collapse=":")
    var_aa1 = gnomad[ig_use[ii[1]],"var"]
    df = data.frame(var=var_aa1, Protein.Consequence=redun_aa_var, Transcript.Consequence=ts_str, Allele.Count=al_count, Allele.Number=al_num,
                    Allele.Frequency=new_af)
    ndf = rbind(ndf, df)
}
print(sprintf("After merging we have %d GnomAD 4.1 variants, Clausen24 had %d", nrow(ndf), sum(! is.na(prkn$gnomad21))))

prkn$gnomad41 = NA
# prkn[match(gnomad[ig_use,"var"],prkn$var),"gnomad41"] = gnomad[ig_use,"Allele.Frequency"]
stopifnot( nrow(ndf) == length(unique(ndf$var)) )
prkn[match(ndf$var,prkn$var),"gnomad41"] = ndf$Allele.Frequency

il = which(! is.na(lene$gnomAD_allele_frequency))
new_var = setdiff(gnomad[ig_use,"var"], lene[il,"variant"])
print(sprintf("New gnomad variants %d: %s", length(new_var), paste0(new_var,collapse=",")))
rm_var = setdiff(lene[il,"variant"], gnomad[ig_use,"var"])
print(sprintf("GnomAD variants from Clausen24 not in GnomAD 4.1 %d: %s", length(rm_var), paste0(rm_var,collapse=",")))

gi_path = ig_use[which(grepl("pathogenic",tolower(gnomad[ig_use,"ClinVar.Germline.Classification"])))]
print(sprintf("Distribution of homozygous counts for %d pathogenic variants", length(gi_path)))
print(summary(gnomad[gi_path,"Homozygote.Count"]/gnomad[gi_path,"Allele.Count"]*100))

print("Pathogenic variants according to gnomad41:")
print(gnomad[gi_path,c("Protein.Consequence","var","Allele.Count","Homozygote.Count","Allele.Frequency","Filters...joint","ClinVar.Germline.Classification")])


##
## Low-throughput validation
##
# This is non-log mean fluorescence of 3 replicates from 2 experiments where each has a WT reference for normalization
# I made a new sheet in the xlsx with mean fluorescence in a csv like setup without normalized values so I check the calc below manually
ltp = read.csv("low_throughput_validation_final_clean.csv")

# Mark each measurement with the WT for normalization
ltp$iwt = NA
ltp[which(ltp$Experiment == 1),"iwt"] = which(ltp$Variant == "WT")
ltp[which(ltp$Experiment == 2),"iwt"] = which(ltp$Variant == "WT2")

# For each replicate, calc difference and normalize to WT
for (rep in seq(3)) {
    cn_diff = sprintf("AO-DMSO.%d",rep)
    ltp[,cn_diff] = ltp[,sprintf("Mean.AO.%d",rep)] - ltp[,sprintf("Mean.DMSO.%d",rep)]
    ltp[,sprintf("AO-DMSO.norm.%d",rep)] = ltp[,cn_diff] / ltp[ltp$iwt,cn_diff]
}
cns = sprintf("AO-DMSO.norm.%d",seq(3))
ltp$mean = apply(ltp[,cns], MARGIN=1, mean)
ltp$sd = apply(ltp[,cns], MARGIN=1, sd)

il2h = match(prkn$var,ltp$Variant)
prkn$ltp = ltp[il2h,"mean"]
prkn$ltp_sd = ltp[il2h,"sd"]


##
## Yi et al 2019
##
yi19 = read.csv("Yi19.csv")

nci = nchar(yi19$X)
yi19$wt = substr(yi19$X, 3, 3)
yi19$mut = substr(yi19$X, nci, nci)
yi19$resi = as.numeric( substr(yi19$X, 4, nci-1) )
yi19[1,"wt"] = ""
yi19[1,"mut"] = ""
yi19[1,"resi"] = 0
yi19$var = sprintf("%s%d%s", yi19$wt, yi19$resi, yi19$mut)
yi19[1,"var"] = "WT"

iy = match(prkn$var, yi19$var)
prkn$yi19 = yi19[iy,"Normalized.mitophagy..."]


##
## AlphaMissense
##
# Cheng23
# https://doi.org/10.5281/zenodo.8208688
# https://zenodo.org/records/8208688/files/AlphaMissense_aa_substitutions.tsv.gz
# zgrep O60260 AlphaMissense_aa_substitutions.tsv.gz > AlphaMissense_O60260.tsv
am = read.csv("AlphaMissense_O60260.tsv", sep="\t", header=F)
colnames(am) = c("uniprot_id","protein_variant","am_pathogenicity","am_class")

i_am2prkn = match(prkn$var, am$protein_variant)
prkn$am_score = am[i_am2prkn,"am_pathogenicity"]
prkn$am_class = am[i_am2prkn,"am_class"]


##
## popEVE
##
# Orenbuch23_popEVE_genetic_disorders_preprint
# https://pop.evemodel.org
pe = read.csv("NP_004553.2.csv")

i_pe2prkn = match(prkn$var, pe$mutant)
prkn$peve = pe[i_pe2prkn,"popEVE"]
prkn$eve = pe[i_pe2prkn,"EVE"]


##
## Save
##
cns = c("var","ao","ao_sd","dmso","dmso_sd","activity","activity_sd","act_keima","act_keima_sd","ltp","ltp_sd","yi19","abun","abun_sd","plddt","rasa","gemme","rosetta","esm1b","esmif","gnomad21","clinvar","gnomad41","am_score","am_class","peve","eve")
write.csv(prkn[,cns], "prkn.csv", row.names=F)

