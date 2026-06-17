#!/usr/bin/env Rscript
# Merge author-approved RoB 2 for the original 80 studies with the de novo
# assessments of the newly ingested gap trials, then recompute the RoB
# distribution and inter-rater agreement (% + Cohen's kappa) over the
# expanded set. The original 80 consensus judgments are NOT altered.
# Moattari (SRMA-0903) is documented but EXCLUDED (non-randomized single arm).

suppressMessages({library(dplyr)})

root <- normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(FALSE),value=TRUE)[1])), "..", ".."), mustWork=FALSE)
if (!dir.exists(file.path(root,"final_analysis","data"))) root <- getwd()

canon <- read.csv(file.path(root,"final_analysis","data","rob2_assessment.csv"), stringsAsFactors=FALSE)
new   <- read.csv(file.path(root,"final_analysis","outputs","rob2_new_trials.csv"), stringsAsFactors=FALSE)

doms <- c("D1_randomization","D2_deviations","D3_missing_data","D4_outcome_measure","D5_selective_report")
keep_cols <- c("study_id","citation","report_id","page_count",
  as.vector(t(outer(doms, c("__reviewer1","__reviewer2","__consensus","__evidence"), paste0))),
  "overall_consensus","verification_status","note")

canon <- canon[, keep_cols]
new   <- new[,   keep_cols]

# New trials added to the broad SR (Moattari excluded: not an RCT)
new_included <- new %>% filter(study_id != "SRMA-0903")
expanded <- bind_rows(canon, new_included)

agree_kappa <- function(df){
  r1 <- unlist(df[paste0(doms,"__reviewer1")]); r2 <- unlist(df[paste0(doms,"__reviewer2")])
  ok <- !is.na(r1)&!is.na(r2)&r1!=""&r2!=""
  r1 <- r1[ok]; r2 <- r2[ok]
  tab <- table(factor(r1, levels=c("Low","Some concerns","High")),
               factor(r2, levels=c("Low","Some concerns","High")))
  po <- sum(diag(tab))/sum(tab); pe <- sum(rowSums(tab)*colSums(tab))/sum(tab)^2
  list(n_studies=nrow(df), n_cells=sum(ok), pct_agree=round(po*100,2), kappa=round((po-pe)/(1-pe),3))
}

dist_tab <- function(df) {
  t <- table(factor(df$overall_consensus, levels=c("Low","Some concerns","High")))
  as.data.frame(t, stringsAsFactors=FALSE) |> setNames(c("overall","n"))
}

orig <- agree_kappa(canon); exp <- agree_kappa(expanded)

cat("================ RoB 2 — ORIGINAL 80 ================\n")
print(dist_tab(canon)); cat(sprintf("studies=%d  per-domain cells=%d  agreement=%.2f%%  kappa=%.3f\n\n",
  orig$n_studies, orig$n_cells, orig$pct_agree, orig$kappa))
cat("================ RoB 2 — EXPANDED (",nrow(expanded),") ================\n")
print(dist_tab(expanded)); cat(sprintf("studies=%d  per-domain cells=%d  agreement=%.2f%%  kappa=%.3f\n\n",
  exp$n_studies, exp$n_cells, exp$pct_agree, exp$kappa))
cat("New trials added:", paste(new_included$study_id, collapse=", "), "\n")
cat("Documented but EXCLUDED (non-RCT):", paste(new$study_id[new$study_id=="SRMA-0903"], collapse=", "), "\n")

# Write merged set + a robvis-ready summary (consensus per domain) for figure rebuild
out_dir <- file.path(root,"final_analysis","outputs"); dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
write.csv(expanded, file.path(out_dir,"rob2_assessment_expanded.csv"), row.names=FALSE)

summ <- data.frame(
  analysis = c("original_80","expanded"),
  n_studies = c(orig$n_studies, exp$n_studies),
  n_low = c(sum(canon$overall_consensus=="Low"), sum(expanded$overall_consensus=="Low")),
  n_some = c(sum(canon$overall_consensus=="Some concerns"), sum(expanded$overall_consensus=="Some concerns")),
  n_high = c(sum(canon$overall_consensus=="High"), sum(expanded$overall_consensus=="High")),
  pct_agreement = c(orig$pct_agree, exp$pct_agree),
  cohen_kappa = c(orig$kappa, exp$kappa)
)
write.csv(summ, file.path(out_dir,"rob2_kappa_summary.csv"), row.names=FALSE)
cat("\nWrote rob2_assessment_expanded.csv and rob2_kappa_summary.csv\n")
