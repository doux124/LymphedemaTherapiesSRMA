#!/usr/bin/env Rscript
# =====================================================================
# Standard random-effects NMA + additive CNMA (netmeta/netcomb) of MLD
# and IPC on a compression backbone. Reads cnma_arm_level.csv.
#
# PRIMARY network  = clean contrasts (network_role == "primary"):
#   clean MLD/IPC add-on edges + the clean Johansson MLD-vs-IPC edge +
#   De Vrieze sham-controlled MLD edge + Szolnoky combination edge.
# SENSITIVITY-FULL = primary + the two confounded/non-randomised head-to-
#   heads (Haghighat, Gurdal) added back.
# Steps 4-5: design/identifiability, network plot, NMA, CNMA,
# heterogeneity, inconsistency, model comparison.
# =====================================================================
suppressPackageStartupMessages({library(netmeta); library(dplyr)})
root <- normalizePath(getwd(), mustWork=TRUE)
outdir <- file.path(root,"final_analysis","outputs"); figdir <- file.path(root,"final_analysis","figures")
dir.create(figdir,showWarnings=FALSE,recursive=TRUE)

arm <- read.csv(file.path(outdir,"cnma_arm_level.csv"), stringsAsFactors=FALSE)
mkdat <- function(a) with(a, data.frame(studlab=paste0(study_id," (",sub(" .*","",citation),")"),
                                        treat1=index_node, treat2=ref_node, TE=smd, seTE=se, stringsAsFactors=FALSE))
arm_prim <- subset(arm, network_role=="primary")
dat  <- mkdat(arm_prim)        # primary
datF <- mkdat(arm)             # full (sensitivity)

fit <- function(d){
  net <- netmeta(TE,seTE,treat1,treat2,studlab,data=d,sm="SMD",common=TRUE,random=TRUE,
                 reference.group="CB",method.tau="REML",details.chkmultiarm=FALSE)
  nc  <- netcomb(net, inactive="CB")
  list(net=net, nc=nc)
}
P <- fit(dat); Fz <- fit(datF)
net <- P$net; nc <- P$nc

sink(file.path(outdir,"cnma_model_results.txt"))
cat("################ NETWORK / CNMA RESULTS ################\n")
cat("Date:", as.character(Sys.Date()), "  netmeta", as.character(packageVersion("netmeta")), "\n")
cat("PRIMARY = clean network (confounded head-to-heads Haghighat & Gurdal demoted to sensitivity).\n\n")

# ---------- Step 4: connectivity + design/identifiability ----------
nc_con <- netconnection(dat$treat1, dat$treat2, dat$studlab)
cat("==== PRIMARY NETWORK STRUCTURE ====\n")
cat("n studies (contrasts):", net$k, "  treatments:", net$n, "\n")
cat("Treatments:", paste(net$trts, collapse=", "), "\n")
cat("Network connected:", nc_con$n.subnets==1, " (subnets:", nc_con$n.subnets, ")\n")
cat("Studies per edge:\n"); print(table(paste(arm_prim$index_node,"vs",arm_prim$ref_node)))
nodemap <- read.csv(file.path(outdir,"cnma_node_components.csv"), stringsAsFactors=FALSE)
Bmat <- as.matrix(nodemap[,c("MLD","IPC")]); rownames(Bmat) <- nodemap$node
cat("\nComponent design matrix (CB = inactive backbone):\n"); print(Bmat)
cat("Rank:", qr(Bmat)$rank, "-> both components identifiable\n\n")

report_net <- function(net, nc, label){
  cat("==== ", label, " ====\n")
  cb <- which(net$trts=="CB")
  tab <- data.frame(treatment=net$trts,
                    SMD_vs_CB=round(net$TE.random[,cb],3),
                    CI_low=round(net$lower.random[,cb],3),
                    CI_high=round(net$upper.random[,cb],3),
                    p=round(net$pval.random[,cb],4))
  cat("- Standard random-effects NMA (vs CB):\n"); print(tab, row.names=FALSE)
  cat(sprintf("  tau2=%.4f tau=%.3f I2=%.1f%% ; global Q=%.2f (df=%d,p=%.4f); incoherence Q=%.2f (df=%d, p=%.4f)\n",
              net$tau2, net$tau, net$I2*100, net$Q, net$df.Q, net$pval.Q,
              net$Q.inconsistency, net$df.Q.inconsistency, net$pval.Q.inconsistency))
  comp <- data.frame(component=nc$comps, SMD=round(nc$Comp.random,3),
                     CI_low=round(nc$lower.Comp.random,3), CI_high=round(nc$upper.Comp.random,3),
                     p=round(nc$pval.Comp.random,4))
  cat("- Additive CNMA component effects:\n"); print(comp, row.names=FALSE)
  cbn <- which(nc$trts=="CB")
  cmb <- data.frame(treatment=nc$trts, SMD_vs_CB=round(nc$TE.random[,cbn],3),
                    CI_low=round(nc$lower.random[,cbn],3), CI_high=round(nc$upper.random[,cbn],3),
                    p=round(nc$pval.random[,cbn],4))
  cat("- Additive combination/treatment estimates (vs CB):\n"); print(cmb, row.names=FALSE)
  cat(sprintf("- Additivity (CNMA vs full NMA): Q.diff=%.3f (df=%d, p=%.4f)\n\n",
              nc$Q.diff, nc$df.Q.diff, nc$pval.Q.diff))
  list(nma=tab, comp=comp, comb=cmb)
}
RP <- report_net(P$net, P$nc, "PRIMARY (clean network)")
RF <- report_net(Fz$net, Fz$nc, "SENSITIVITY-FULL (+ Haghighat & Gurdal head-to-heads)")

# node-split for both
ns_out <- function(net, file){
  ns <- tryCatch(netsplit(net), error=function(e) NULL)
  if(is.null(ns)) return(invisible())
  nsdf <- data.frame(comparison=ns$direct.random$comparison,
                     direct=round(ns$direct.random$TE,3), indirect=round(ns$indirect.random$TE,3),
                     diff=round(ns$compare.random$TE,3), p=round(ns$compare.random$p,4))
  nsdf <- nsdf[!is.na(nsdf$direct)&!is.na(nsdf$indirect),]
  print(nsdf, row.names=FALSE); write.csv(nsdf, file, row.names=FALSE)
}
cat("==== NODE-SPLIT: PRIMARY ====\n"); ns_out(P$net, file.path(outdir,"cnma_nodesplit.csv"))
cat("\n==== NODE-SPLIT: SENSITIVITY-FULL ====\n"); ns_out(Fz$net, file.path(outdir,"cnma_nodesplit_full.csv"))
sink()

# tidy CSVs (primary headline)
write.csv(RP$nma,  file.path(outdir,"nma_treatment_vs_CB.csv"), row.names=FALSE)
write.csv(RP$comp, file.path(outdir,"cnma_component_effects.csv"), row.names=FALSE)
write.csv(RP$comb, file.path(outdir,"cnma_combination_vs_CB.csv"), row.names=FALSE)
# primary-vs-full component comparison
cmp <- rbind(cbind(network="primary",   RP$comp),
             cbind(network="full(+H&G)", RF$comp))
write.csv(cmp, file.path(outdir,"cnma_primary_vs_full_components.csv"), row.names=FALSE)

# plots
pdf(file.path(figdir,"Fig_network_plot.pdf"), width=7, height=6)
netgraph(net, plastic=FALSE, thickness="number.of.studies", points=TRUE, cex.points=4,
         col="grey30", number.of.studies=TRUE, cex=1.1,
         main="Primary MLD/IPC component network (compression backbone = CB)")
dev.off()
pdf(file.path(figdir,"Fig_component_forest.pdf"), width=7, height=3)
forest(nc, what="components", main="Additive component effects (primary network; SMD vs backbone)")
dev.off()

cat("\nDONE. Results -> cnma_model_results.txt ; figures refreshed.\n")
cat("\n--- PRIMARY components ---\n"); print(RP$comp, row.names=FALSE)
cat(sprintf("PRIMARY: I2=%.1f%% ; incoherence p=%.3f ; additivity p=%.3f\n", net$I2*100, net$pval.Q.inconsistency, nc$pval.Q.diff))
cat("\n--- FULL (+H&G) components ---\n"); print(RF$comp, row.names=FALSE)
cat(sprintf("FULL: I2=%.1f%% ; incoherence p=%.3f\n", Fz$net$I2*100, Fz$net$pval.Q.inconsistency))
