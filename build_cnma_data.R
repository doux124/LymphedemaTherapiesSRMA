#!/usr/bin/env Rscript
# =====================================================================
# Build the arm-level CNMA dataset + imputation log for the additive
# component network meta-analysis of MLD and IPC on a compression
# backbone (estimand: Hedges' g SMD of a limb-volume outcome at end of
# intervention; higher SMD = more volume reduction = better).
#
# Every number is page-anchored to the source PDF/full text. SD
# imputation follows the prespecified hierarchy (REFINED_ANALYSIS_PLAN
# section 4.3): reported SD -> SE/CI/p -> median/IQR or range (Wan 2014,
# Luo 2018) -> change-score SD via borrowed pre/post r -> borrow SD.
#
# Nodes (compression backbone = reference, NOT an estimable component):
#   CB           compression backbone only
#   CB+MLD       + manual lymphatic drainage
#   CB+IPC       + intermittent pneumatic compression
#   CB+MLD+IPC   + both (combination)
# =====================================================================
suppressPackageStartupMessages({library(metafor); library(dplyr)})
root <- normalizePath(getwd(), mustWork=TRUE)
outdir <- file.path(root,"final_analysis","outputs"); dir.create(outdir,showWarnings=FALSE,recursive=TRUE)

## ---- Wan 2014 / Luo 2018 median-based conversions -------------------
# Scenario S1: min (a), median (m), max (b), n  -> mean, SD
wan_S1 <- function(a,m,b,n){
  mean <- (a+2*m+b)/4 + (a-2*m+b)/(4*n)            # Luo 2018 mean
  xi   <- 2*qnorm((n-0.375)/(n+0.25))              # Wan 2014 SD denominator
  sd   <- (b-a)/xi
  c(mean=mean, sd=sd)
}
# Scenario S2: q1, median (m), q3, n  -> mean, SD
wan_S2 <- function(q1,m,q3,n){
  mean <- (q1+m+q3)/3                              # Luo 2018 mean
  eta  <- 2*qnorm((0.75*n-0.125)/(n+0.25))         # Wan 2014 SD denominator
  sd   <- (q3-q1)/eta
  c(mean=mean, sd=sd)
}
# change-score SD from pre & post SD and assumed correlation r
sd_change <- function(sd_pre, sd_post, r) sqrt(sd_pre^2 + sd_post^2 - 2*r*sd_pre*sd_post)

BASE_R <- 0.5     # base-case pre/post correlation (sensitivity: 0.3, 0.7)

imp <- list()  # imputation log rows
logimp <- function(study, arm, method, detail, page) {
  imp[[length(imp)+1]] <<- data.frame(study_id=study, arm=arm, hierarchy_step=method,
                                      detail=detail, source_page=page, stringsAsFactors=FALSE)
}

# container for arm-level rows
A <- list()
add <- function(study_id,citation,edge,index_node,ref_node,
                n1,m1,sd1,n2,m2,sd2,outcome,direction,timepoint,phase,
                sd_source,page,flags="") {
  A[[length(A)+1]] <<- data.frame(
    study_id,citation,edge,index_node,ref_node,
    n_index=n1, mean_index=m1, sd_index=sd1,
    n_ref=n2,   mean_ref=m2,   sd_ref=sd2,
    outcome_construct=outcome, scale_direction=direction,
    timepoint, phase, sd_source, source_page=page, flags, stringsAsFactors=FALSE)
}

# ============== MLD add-on edges (CB+MLD vs CB) ======================
# McNeely 2004 (SRMA-0457): PREV; arm means+SD reported. p.8 Table 3.
add("SRMA-0457","McNeely 2004","MLD_vs_CB","CB+MLD","CB",
    23,46.1,22.6, 21,38.6,16.1,
    "PREV (% reduction excess volume)","higher_better","end intensive tx","intensive",
    "reported arm SD","p.8 Table 3")

# Tambour 2018 (SRMA-0243): PREV; only between-group MD + 95% CI. p.8 primary table.
# MD = +1.9 (benefit dir), 95% CI -19 to 22.9 -> SE 10.689; pooled SD from SE & n.
se_tam <- (22.9-(-19))/(2*1.96)
sdpool_tam <- se_tam/sqrt(1/39+1/38)
logimp("SRMA-0243","both","2: SE/CI -> SD",
       sprintf("Between-group MD +1.9, 95%% CI -19..22.9 -> SE=%.2f; pooled SD=%.1f from SE and n=39/38",se_tam,sdpool_tam),
       "p.8 primary-outcome table")
add("SRMA-0243","Tambour 2018","MLD_vs_CB","CB+MLD","CB",
    39,1.9,sdpool_tam, 38,0.0,sdpool_tam,
    "PREV (between-group MD representation)","higher_better","1 month","intensive",
    "derived from reported 95% CI of between-group MD (pooled SD)","p.8 primary table")

# Andersen 2000 (SRMA-0901): % reduction absolute edema volume; group means + 95% CI. p.6 (means/CI), p.5 (n).
sd_and_mld <- (65-32)/(2*1.96)*sqrt(20)
sd_and_cb  <- (78-43)/(2*1.96)*sqrt(22)
logimp("SRMA-0901","CB+MLD","2: SE/CI -> SD",
       sprintf("Mean 48%%, 95%% CI 32-65, n=20 -> group SD=%.1f",sd_and_mld),"p.6 Results")
logimp("SRMA-0901","CB","2: SE/CI -> SD",
       sprintf("Mean 60%%, 95%% CI 43-78, n=22 -> group SD=%.1f",sd_and_cb),"p.6 Results")
add("SRMA-0901","Andersen 2000","MLD_vs_CB","CB+MLD","CB",
    20,48.0,sd_and_mld, 22,60.0,sd_and_cb,
    "% reduction in absolute edema volume","higher_better","3 months","intensive (garment)",
    "group SDs derived from reported 95% CIs of group mean % reduction","p.6 Results; p.5 Table 2 (n)",
    "MLD negative (did not contribute)")

# De Vrieze 2022 (SRMA-0144): traditional MLD vs PLACEBO (sham) MLD on a DLT backbone;
# primary endpoint = reduction in % excess arm/hand volume at end of intensive phase.
# Between-group MD = -0.2 pp, 95% CI -2.1 to 1.8 -> SE/CI conversion. p.1 Abstract.
se_dv <- (1.8-(-2.1))/(2*1.96); sdpool_dv <- se_dv/sqrt(1/64+1/65)
logimp("SRMA-0144","both","2: SE/CI -> SD",
       sprintf("Traditional MLD vs placebo MLD MD -0.2 pp, 95%% CI -2.1..1.8 -> SE=%.2f; pooled SD=%.2f from SE and n=64/65",se_dv,sdpool_dv),
       "p.1 Abstract (primary outcome, end of intensive phase)")
add("SRMA-0144","De Vrieze 2022","MLD_vs_CB","CB+MLD","CB",
    64,-0.2,sdpool_dv, 65,0.0,sdpool_dv,
    "reduction in % excess arm/hand volume","higher_better","end intensive phase","intensive (DLT)",
    "derived from reported 95% CI of between-group MD (pooled SD); comparator = placebo/sham MLD","p.1 Abstract; p.4 (primary endpoint)",
    "sham-controlled MLD contrast (largest, cleanest MLD-vs-no-MLD edge, n=129)")

# ============== IPC add-on edges (CB+IPC vs CB) ======================
# Szuba 2002 (SRMA-0469): PREV; arm means+SD reported. p.4 Fig 1.
add("SRMA-0469","Szuba 2002","IPC_vs_CB","CB+IPC","CB",
    12,45.3,18.2, 11,26.0,22.1,
    "PREV (% reduction excess volume)","higher_better","end 2-wk DLT","intensive",
    "reported arm SD","p.4 Figure 1")

# Uzkeser 2015 (SRMA-0904): limb volume (water displacement) median(range); change from baseline. p.4.
uz_ipc_base <- wan_S1(220,840,3460,15); uz_ipc_end <- wan_S1(60,500,2160,15)
uz_cb_base  <- wan_S1(180,630,1820,15); uz_cb_end  <- wan_S1(0,480,1410,15)
red_ipc <- uz_ipc_base["mean"]-uz_ipc_end["mean"]      # reduction (benefit)
red_cb  <- uz_cb_base["mean"] -uz_cb_end["mean"]
sdr_ipc <- sd_change(uz_ipc_base["sd"], uz_ipc_end["sd"], BASE_R)
sdr_cb  <- sd_change(uz_cb_base["sd"],  uz_cb_end["sd"],  BASE_R)
logimp("SRMA-0904","CB+IPC","3: median/range -> Wan2014; 4: change SD via r",
       sprintf("baseline 840(220-3460)->mean %.0f sd %.0f; endpoint 500(60-2160)->mean %.0f sd %.0f; reduction %.0f, change-SD(r=0.5) %.0f",
               uz_ipc_base["mean"],uz_ipc_base["sd"],uz_ipc_end["mean"],uz_ipc_end["sd"],red_ipc,sdr_ipc),"p.4 Results")
logimp("SRMA-0904","CB","3: median/range -> Wan2014; 4: change SD via r",
       sprintf("baseline 630(180-1820)->mean %.0f sd %.0f; endpoint 480(0-1410)->mean %.0f sd %.0f; reduction %.0f, change-SD(r=0.5) %.0f",
               uz_cb_base["mean"],uz_cb_base["sd"],uz_cb_end["mean"],uz_cb_end["sd"],red_cb,sdr_cb),"p.4 Results")
add("SRMA-0904","Uzkeser 2015","IPC_vs_CB","CB+IPC","CB",
    15,as.numeric(red_ipc),as.numeric(sdr_ipc), 15,as.numeric(red_cb),as.numeric(sdr_cb),
    "limb-volume reduction (water displacement, mL)","higher_better","3 weeks","intensive (CDT)",
    "median/range -> Wan2014 mean/SD (baseline & endpoint); change-score SD via r=0.5","p.4 Results",
    "IPC negative (pump did not contribute); baseline imbalance handled via change scores")

# Tastaban 2020 (SRMA-0210): PREV median(IQR). p.5 Table 2.
ta_ipc <- wan_S2(48.3,54.6,58.2,38); ta_cb <- wan_S2(48.3,49.6,54.9,38)
logimp("SRMA-0210","CB+IPC","3: median/IQR -> Wan2014/Luo2018",
       sprintf("PREV 54.6 (IQR 48.3-58.2), n=38 -> mean %.1f sd %.1f",ta_ipc["mean"],ta_ipc["sd"]),"p.5 Table 2")
logimp("SRMA-0210","CB","3: median/IQR -> Wan2014/Luo2018",
       sprintf("PREV 49.6 (IQR 48.3-54.9), n=38 -> mean %.1f sd %.1f",ta_cb["mean"],ta_cb["sd"]),"p.5 Table 2")
add("SRMA-0210","Tastaban 2020","IPC_vs_CB","CB+IPC","CB",
    38,as.numeric(ta_ipc["mean"]),as.numeric(ta_ipc["sd"]), 38,as.numeric(ta_cb["mean"]),as.numeric(ta_cb["sd"]),
    "PREV (% reduction excess volume)","higher_better","4 weeks (20 sessions)","intensive (CDT)",
    "median/IQR -> Wan2014/Luo2018 mean/SD","p.5 Table 2")

# ============== Direct MLD vs IPC edges =============================
# Johansson 1998 (SRMA-0488): lymphedema volume mL; Part II change; pre/post SD reported. p.5 Table 3.
joh_mld_red <- 579-504; joh_spc_red <- 411-382
joh_mld_sd  <- sd_change(258,252,BASE_R); joh_spc_sd <- sd_change(203,193,BASE_R)
logimp("SRMA-0488","CB+MLD","4: change SD via r",
       sprintf("test2 579+/-258 -> test3 504+/-252 (n=12); reduction %d, change-SD(r=0.5) %.0f",joh_mld_red,joh_mld_sd),"p.5 Table 3")
logimp("SRMA-0488","CB+IPC","4: change SD via r",
       sprintf("test2 411+/-203 -> test3 382+/-193 (n=12); reduction %d, change-SD(r=0.5) %.0f",joh_spc_red,joh_spc_sd),"p.5 Table 3")
add("SRMA-0488","Johansson 1998","MLD_vs_IPC","CB+MLD","CB+IPC",
    12,joh_mld_red,joh_mld_sd, 12,joh_spc_red,joh_spc_sd,
    "lymphedema-volume reduction (mL)","higher_better","Part II (after MLD/SPC)","intensive",
    "change-score SD from reported test2/test3 SDs via r=0.5","p.5 Table 3",
    "SPC = sequential pneumatic compression (IPC)")

# Haghighat 2010 (SRMA-0405): PREV; arm means+SD reported. p.5 Table 3. IPC arm replaced most MLD.
add("SRMA-0405","Haghighat 2010","IPC_vs_MLD","CB+IPC","CB+MLD",
    56,37.5,14.4, 56,43.1,13.7,
    "PREV (% reduction excess volume)","higher_better","end intensive tx","intensive",
    "reported arm SD","p.5 Table 3",
    "confounded: IPC replaced most arm MLD dose -> excluded in clean-add-on sensitivity")

# Gurdal 2012 (SRMA-0905): total arm volume mL; reduction; baseline SD reported, post SD imputed=baseline. p.3-4.
gur_mld_red <- 3533-3004; gur_ipc_red <- 3581-3142
gur_mld_sd  <- sd_change(739,739,BASE_R); gur_ipc_sd <- sd_change(783,783,BASE_R)  # post SD imputed = baseline SD
logimp("SRMA-0905","CB+MLD","4/5: post SD imputed=baseline; change SD via r",
       sprintf("baseline 3533+/-739 -> endpoint 3004 (post SD imputed=739, n=15); reduction %d, change-SD(r=0.5) %.0f",gur_mld_red,gur_mld_sd),"p.3-4 Table 3")
logimp("SRMA-0905","CB+IPC","4/5: post SD imputed=baseline; change SD via r",
       sprintf("baseline 3581+/-783 -> endpoint 3142 (post SD imputed=783, n=15); reduction %d, change-SD(r=0.5) %.0f",gur_ipc_red,gur_ipc_sd),"p.3-4 Table 3")
add("SRMA-0905","Gurdal 2012","MLD_vs_IPC","CB+MLD","CB+IPC",
    15,gur_mld_red,gur_mld_sd, 15,gur_ipc_red,gur_ipc_sd,
    "total-arm-volume reduction (mL)","higher_better","6 weeks","intensive",
    "post-treatment SD imputed = reported baseline SD; change-score SD via r=0.5","p.3 (baseline SD); p.4 Table 3 (means)",
    "head-to-head + randomization contradiction + SLD confound -> excluded in clean-add-on sensitivity")

# ============== Combination edge (CB+MLD+IPC vs CB+MLD) =============
# Szolnoky 2009 (SRMA-0407): % total arm volume reduction; NO dispersion reported.
# DOSE-SUBSTITUTION CONTRAST: the arms are 60 min MLD vs 30 min MLD + 30 min IPC
# (source p.3) -- IPC was substituted for half the MLD session, not added to an
# identical regimen. It is coded here as adding IPC to a CB+MLD backbone, which is
# valid only under a dose-equivalence assumption (30 == 60 min MLD). It is RETAINED
# in the primary network because it is the SOLE within-study combination contrast,
# i.e. the only anchor for the CB+MLD+IPC node and the only direct test of
# additivity; the assumption is justified by the precise-null MLD component
# (bounding the 30-min reduction near zero) and examined in the full-network
# sensitivity (which re-adds Haghighat/Gurdal). Audit & author-approved 2026-06-17.
# Impute SD from the reported between-group p<0.05 evaluated at the p=0.05 boundary
# (largest SD consistent with the reported significance -> most conservative / widest CI).
sz_diff <- 7.93-3.06
df_sz <- 13+14-2
t_crit <- qt(0.975, df_sz)
se_diff_sz <- sz_diff/t_crit
sdpool_sz <- se_diff_sz/sqrt(1/14+1/13)
logimp("SRMA-0407","both","5: borrow/derive SD from reported p<0.05 boundary",
       sprintf("means 7.93%% (n=14) vs 3.06%% (n=13), between-group p<0.05; at p=0.05 boundary t=%.2f -> SE_diff=%.2f -> pooled SD=%.2f (arm SDs set equal)",t_crit,se_diff_sz,sdpool_sz),
       "p.4 Results (means + p<0.05)")
add("SRMA-0407","Szolnoky 2009","MLDIPC_vs_MLD","CB+MLD+IPC","CB+MLD",
    14,7.93,sdpool_sz, 13,3.06,sdpool_sz,
    "% reduction total arm volume","higher_better","end of therapy","intensive (CDP)",
    "no dispersion reported; SD derived from reported p<0.05 at the boundary (most conservative)","p.4 Results",
    "dose-substitution combination edge (60min MLD vs 30min MLD+30min IPC); sole combination-node anchor, retained under dose-equivalence assumption; weakest imputation; excluded in 'exclude imputed-variance' sensitivity")

arm <- bind_rows(A)
implog <- bind_rows(imp)

# network role: the two confounded/non-randomised head-to-heads (Haghighat =
# IPC-replaces-MLD; Gurdal = randomization contradiction + SLD confound) are
# DEMOTED to a sensitivity role because the clean MLD-vs-IPC edge they would add is
# already supplied without dose substitution by Johansson. Szolnoky is ALSO a
# dose-substitution contrast but is RETAINED in primary because it is the only
# evidence at the combination node (no clean alternative exists anywhere) -- see the
# block above. The primary network is therefore eight clean add-on/component
# contrasts plus the one combination contrast (Szolnoky). Audit-reviewed 2026-06-17.
arm$network_role <- ifelse(arm$study_id %in% c("SRMA-0405","SRMA-0905"),
                           "sensitivity_confounded", "primary")

## ---- compute per-contrast SMD (Hedges g) + SE -----------------------
es <- escalc(measure="SMD",
             m1i=mean_index, sd1i=sd_index, n1i=n_index,
             m2i=mean_ref,   sd2i=sd_ref,   n2i=n_ref, data=arm)
arm$smd  <- as.numeric(es$yi)      # Hedges g, index vs ref (positive = index better)
arm$se   <- sqrt(as.numeric(es$vi))
arm$smd_ci_low  <- arm$smd - 1.96*arm$se
arm$smd_ci_high <- arm$smd + 1.96*arm$se

## ---- node component design map -------------------------------------
nodes <- data.frame(
  node = c("CB","CB+MLD","CB+IPC","CB+MLD+IPC"),
  MLD  = c(0,1,0,1), IPC = c(0,0,1,1), stringsAsFactors=FALSE)

write.csv(arm,    file.path(outdir,"cnma_arm_level.csv"), row.names=FALSE)
write.csv(implog, file.path(outdir,"imputation_log.csv"), row.names=FALSE)
write.csv(nodes,  file.path(outdir,"cnma_node_components.csv"), row.names=FALSE)

cat("=== cnma_arm_level.csv (",nrow(arm)," contrasts) ===\n")
print(arm[,c("study_id","edge","index_node","ref_node","smd","se","smd_ci_low","smd_ci_high")], digits=3)
cat("\n=== imputation_log.csv (",nrow(implog)," entries) ===\n")
cat("Edges per type:\n"); print(table(arm$edge))
cat("\nTrials with imputed variance:",
    paste(unique(implog$study_id), collapse=", "), "\n")
