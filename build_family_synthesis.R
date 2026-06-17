#!/usr/bin/env Rscript

# Build the structured, full-corpus, non-pooled synthesis dataset for all 80
# included studies, organised by intervention family. Combines:
#   - canonical all-study RoB 2 (data/rob2_assessment.csv)
#   - full-text quantitative reassessment family/role (data/full_text_quantitative_reassessment.csv)
#   - source-verified component-contrast effects (data/component_contrast_audit.csv)
#   - source-verified priority non-pooled extractions (data/priority_nonpooled_extractions.csv)
#   - source-verified pooled effects (verified_primary_effects.csv, verified_exercise_effects.csv)
# Outputs:
#   - data/structured_family_synthesis.csv  (one row per study)
#   - data/family_synthesis_summary.csv     (one row per intervention family)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
})

root <- normalizePath(getwd(), mustWork = TRUE)

rob   <- read_csv(file.path(root, "final_analysis/data/rob2_assessment.csv"), show_col_types = FALSE)
reass <- read_csv(file.path(root, "final_analysis/data/full_text_quantitative_reassessment.csv"), show_col_types = FALSE)
audit <- read_csv(file.path(root, "final_analysis/data/component_contrast_audit.csv"), show_col_types = FALSE)
prio  <- read_csv(file.path(root, "final_analysis/data/priority_nonpooled_extractions.csv"), show_col_types = FALSE)
emap  <- read_csv(file.path(root, "final_analysis/data/evidence_map_extractions.csv"), show_col_types = FALSE)
vpe   <- read_csv(file.path(root, "final_analysis/data/verified_primary_effects.csv"), show_col_types = FALSE)
vee   <- read_csv(file.path(root, "final_analysis/data/verified_exercise_effects.csv"), show_col_types = FALSE)
reg   <- read_csv(file.path(root, "final_analysis/data/study_register.csv"), show_col_types = FALSE)

# ---- Map the 8 reassessment families to the 7 manuscript families ----------
family_map <- c(
  "MLD or drainage comparison"                                              = "MLD / drainage",
  "IPC or pneumatic-device comparison"                                      = "IPC / pneumatic devices",
  "Exercise versus usual care: compatible relative interlimb-volume scale"  = "Exercise / movement",
  "Exercise or movement: incompatible comparator, estimand, or reporting"   = "Exercise / movement",
  "Compression modality, dose, or taping comparison"                        = "Compression modality / dose / taping",
  "Nighttime compression add-on"                                            = "Nighttime compression",
  "Laser or photobiomodulation comparison"                                  = "Laser / photobiomodulation",
  "Multicomponent, self-care, or structural comparison"                     = "Multicomponent / self-care / structural"
)
family_order <- c(
  "MLD / drainage", "IPC / pneumatic devices", "Exercise / movement",
  "Compression modality / dose / taping", "Nighttime compression",
  "Laser / photobiomodulation", "Multicomponent / self-care / structural"
)

base <- reass %>%
  transmute(
    study_id, year, title,
    family_detail = synthesis_family,
    family = recode(synthesis_family, !!!family_map),
    recommended_role = recommended_role,
    reassessment_reason = reason
  ) %>%
  left_join(rob %>% select(study_id, citation, rob_overall = overall_consensus), by = "study_id")

# ---- Per-study effect / comparison detail ----------------------------------
# Priority extractions (source-verified)
prio_detail <- prio %>%
  transmute(study_id,
            design,
            comparison,
            outcome_scale, timepoint,
            effect_summary,
            uncertainty,
            source_anchor,
            pooling_status = pooling_status,
            pooling_reason)

# Component-contrast audit (source-verified). Uncertainty is parsed from the
# audited effect summary (CI/SE/median descriptors); the timepoint is the
# review's central audited volume contrast (end of randomized intensive treatment).
extract_uncertainty <- function(x) {
  ci  <- regmatches(x, regexpr("95% ?CI[^;]*", x))
  se  <- regmatches(x, regexpr("SE ?[0-9.]+", x))
  med <- ifelse(grepl("median|IQR", x, ignore.case = TRUE), "medians/IQR only", "")
  out <- paste(c(ci, se, med)[nzchar(c(ci, se, med))], collapse = "; ")
  ifelse(nzchar(out), out, "Not separately reported (see effect summary)")
}
audit_detail <- audit %>%
  transmute(study_id,
            design = "Randomized controlled trial",
            comparison = component_identification,
            outcome_scale = outcome_estimand,
            timepoint = "End of randomized intensive treatment (central audited volume contrast)",
            effect_summary,
            uncertainty = vapply(effect_summary, extract_uncertainty, character(1)),
            source_anchor = source_location,
            pooling_status = ifelse(primary_pool_status == "Included",
                                    "Pooled (primary MLD)", "Not pooled"),
            pooling_reason = ifelse(is.na(reason_not_primary_pool) | reason_not_primary_pool == "",
                                    "Entered the primary MLD add-on PREV pool", reason_not_primary_pool))

# Evidence-map full-text extractions (source-verified, page-anchored)
em_detail <- emap %>%
  transmute(study_id, design, comparison, outcome_scale, timepoint,
            effect_summary, uncertainty, source_anchor,
            pooling_status, pooling_reason)

# Pooled exercise effects (source-verified)
ex_detail <- vee %>%
  transmute(study_id,
            design = NA_character_,
            comparison,
            outcome_scale = outcome, timepoint = analysis_timepoint,
            effect_summary = sprintf("%+0.2f pp (benefit direction)", effect_pp),
            uncertainty = sprintf("95%% CI %0.2f to %0.2f; SE %0.3f", ci_low, ci_high, se),
            source_anchor = source_page_table_figure,
            pooling_status = "Pooled (secondary exercise)",
            pooling_reason = "Entered the secondary exercise vs usual-care relative interlimb-volume pool")

# Precedence: priority > exercise pool > evidence-map > component audit (disjoint)
detail <- bind_rows(prio_detail, ex_detail, em_detail, audit_detail) %>%
  group_by(study_id) %>% slice(1) %>% ungroup()

structured <- base %>%
  left_join(detail, by = "study_id") %>%
  mutate(
    design = ifelse(is.na(design) | design == "", "Randomized controlled trial", design),
    pooling_status = ifelse(is.na(pooling_status), "Not pooled", pooling_status),
    pooling_reason = ifelse(is.na(pooling_reason) | pooling_reason == "",
                            reassessment_reason, pooling_reason),
    family = factor(family, levels = family_order)
  ) %>%
  arrange(family, desc(pooling_status == "Pooled (primary MLD)"),
          desc(pooling_status == "Pooled (secondary exercise)"), citation) %>%
  select(study_id, citation, year, family, family_detail, design, comparison,
         outcome_scale, timepoint, effect_summary, uncertainty,
         rob_overall, recommended_role, pooling_status, pooling_reason, source_anchor)

stopifnot(nrow(structured) == 80, length(unique(structured$study_id)) == 80)
stopifnot(all(!is.na(structured$rob_overall)))

# completeness: every study source-verified across the core fields (no NA/placeholders)
core <- c("comparison", "outcome_scale", "timepoint", "effect_summary", "source_anchor")
for (col in core) {
  bad <- structured$study_id[is.na(structured[[col]]) | structured[[col]] %in%
           c("", "NA", "See Supplementary Table 2 (characteristics)")]
  if (length(bad)) stop("Incomplete ", col, " for: ", paste(bad, collapse = ", "))
}
# every source anchor cites a page number
nopage <- structured$study_id[!grepl("p\\.\\s?[0-9]", structured$source_anchor)]
if (length(nopage)) stop("Source anchor without a page number for: ", paste(nopage, collapse = ", "))

write_csv(structured, file.path(root, "final_analysis/data/structured_family_synthesis.csv"))

# ---- Family-level summary ---------------------------------------------------
fam_summary <- structured %>%
  group_by(family) %>%
  summarise(
    n_studies = n(),
    rob_high = sum(rob_overall == "High"),
    rob_some = sum(rob_overall == "Some concerns"),
    rob_low  = sum(rob_overall == "Low"),
    n_quantitative = sum(str_detect(pooling_status, "Pooled")),
    n_source_verified_nonpooled = sum(recommended_role %in%
        c("Source-verified non-pooled evidence", "Priority for source-verified non-pooled extraction")),
    n_evidence_map = sum(recommended_role == "Evidence-map row only"),
    .groups = "drop"
  ) %>%
  mutate(rob_pattern = sprintf("%d high / %d some concerns / %d low", rob_high, rob_some, rob_low))

pooling_rationale <- c(
  "MLD / drainage" = "Two clean MLD add-on PREV trials with a source-verified mean and variance were pooled (primary analysis); the remaining drainage trials used incompatible comparators, estimands, or lacked an extractable PREV variance.",
  "IPC / pneumatic devices" = "Not pooled: only one clean IPC add-on PREV trial reported a usable mean and variance (Szuba 2002); other IPC contrasts changed MLD dose together, reported medians/IQRs, or used a different estimand.",
  "Exercise / movement" = "Three exercise-versus-usual-care trials sharing a relative interlimb-volume percentage-point scale were pooled as a secondary/exploratory analysis; the remaining exercise trials used incompatible comparators, doses, timepoints, or summary statistics.",
  "Compression modality / dose / taping" = "Not pooled: trials compared compression types, doses, or taping rather than a common added component, across heterogeneous phases and outcome scales (raw mL, RVC, medians).",
  "Nighttime compression" = "Not pooled: clean nighttime-compression contrasts exist but do not all report a common effect and variance on the same scale (relative-change, rebound, and adjusted between-group estimands).",
  "Laser / photobiomodulation" = "Not pooled: laser trials use heterogeneous co-interventions, timepoints, and raw-volume/circumference outcomes, frequently without an extractable between-group variance.",
  "Multicomponent / self-care / structural" = "Not pooled: multicomponent and structural contrasts differ in several components at once and report incompatible outcome scales."
)
fam_summary <- fam_summary %>%
  mutate(pooling_rationale = pooling_rationale[as.character(family)]) %>%
  arrange(family)

write_csv(fam_summary, file.path(root, "final_analysis/data/family_synthesis_summary.csv"))

cat("Structured family synthesis written:\n")
cat(sprintf("  data/structured_family_synthesis.csv (%d studies)\n", nrow(structured)))
cat(sprintf("  data/family_synthesis_summary.csv (%d families)\n", nrow(fam_summary)))
cat("\nFamily-level RoB and roles:\n")
print(as.data.frame(fam_summary %>% select(family, n_studies, rob_pattern,
       n_quantitative, n_source_verified_nonpooled, n_evidence_map)))
cat(sprintf("\nTotals: %d studies; RoB %d high / %d some / %d low; %d quantitative\n",
    sum(fam_summary$n_studies), sum(fam_summary$rob_high), sum(fam_summary$rob_some),
    sum(fam_summary$rob_low), sum(fam_summary$n_quantitative)))
