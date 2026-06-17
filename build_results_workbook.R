#!/usr/bin/env Rscript
# Assemble the CNMA results workbook (one sheet per analysis) from the
# reproducible CSV outputs in final_analysis/outputs/.
suppressPackageStartupMessages({library(openxlsx)})
root <- normalizePath(getwd(), mustWork=TRUE)
outdir <- file.path(root,"final_analysis","outputs")
rd <- function(f) tryCatch(read.csv(file.path(outdir,f), stringsAsFactors=FALSE, check.names=FALSE), error=function(e) data.frame(note=paste("missing:",f)))

overview <- data.frame(
  item=c("Title","Estimand","Network nodes","Components (estimable)","Reference (inactive)",
         "PRIMARY network","Demoted to sensitivity","New gap trials added","Expanded RoB set","Kappa (expanded)",
         "PRIMARY component IPC","PRIMARY component MLD","PRIMARY incoherence / I2","Full network (sensitivity)","Additivity test","Status"),
  value=c("Additive CNMA of MLD and IPC on a compression backbone (BCRL, JCM)",
          "Hedges' g SMD of a limb-volume outcome at end of intervention (higher = more reduction)",
          "CB, CB+MLD, CB+IPC, CB+MLD+IPC",
          "MLD, IPC","CB (compression backbone)",
          "9 clean contrasts (incl. De Vrieze sham-controlled MLD; clean Johansson MLD-vs-IPC edge)",
          "Haghighat & Gurdal (confounded/non-randomised head-to-heads)",
          "Andersen 2000, Uzkeser 2015, Gurdal 2012 (network); Dini 1998 (SWiM/benchmark); Moattari 2012 EXCLUDED (single-arm)",
          "84 studies (80 author-approved + 4 new; Moattari excluded)","98.57% agreement, Cohen kappa 0.968",
          "SMD 0.46 (95% CI 0.17-0.76, p=0.002) -- significant, robust",
          "SMD 0.06 (95% CI -0.17-0.28, p=0.61) -- precise null",
          "incoherence p=0.13 (resolved); I2 = 0%",
          "adding Haghighat+Gurdal -> incoherence p=0.002, I2=42% (confirms they are the distortion)",
          "Q.diff p=0.41 (additivity acceptable)",
          "RESULTS ONLY - no manuscript prose written"))

sheets <- list(
  Overview              = overview,
  New_trials_provenance = rd("new_trials_provenance.csv"),
  RoB2_new_trials       = rd("rob2_new_trials.csv"),
  RoB2_kappa_summary    = rd("rob2_kappa_summary.csv"),
  RoB2_expanded_84      = rd("rob2_assessment_expanded.csv"),
  CNMA_arm_level        = rd("cnma_arm_level.csv"),
  Imputation_log        = rd("imputation_log.csv"),
  Node_components       = rd("cnma_node_components.csv"),
  NMA_treatment_vs_CB   = rd("nma_treatment_vs_CB.csv"),
  CNMA_component_effects= rd("cnma_component_effects.csv"),
  CNMA_combination      = rd("cnma_combination_vs_CB.csv"),
  Primary_vs_full_comp  = rd("cnma_primary_vs_full_components.csv"),
  Node_split_primary    = rd("cnma_nodesplit.csv"),
  Node_split_full       = rd("cnma_nodesplit_full.csv"),
  Benchmark             = rd("benchmark_table.csv"),
  Transitivity          = rd("transitivity_table.csv"),
  Sensitivity           = rd("cnma_sensitivity.csv"),
  Certainty_CINeMA      = rd("certainty_cinema_grade.csv")
)

wb <- createWorkbook()
hs <- createStyle(textDecoration="bold", fgFill="#1F4E79", fontColour="white", halign="left")
for(nm in names(sheets)){
  addWorksheet(wb, nm)
  writeData(wb, nm, sheets[[nm]], headerStyle=hs)
  freezePane(wb, nm, firstRow=TRUE)
  setColWidths(wb, nm, cols=1:max(1,ncol(sheets[[nm]])), widths="auto")
}
saveWorkbook(wb, file.path(outdir,"CNMA_results_workbook.xlsx"), overwrite=TRUE)
cat("Wrote CNMA_results_workbook.xlsx with", length(sheets), "sheets:\n  ", paste(names(sheets),collapse=", "), "\n")
