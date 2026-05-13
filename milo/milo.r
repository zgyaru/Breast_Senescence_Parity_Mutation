library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(ggbeeswarm)

source("../Figure2/milo_func/helper.r")

match_dict <- c(
  "old_Parous_WT_Aging_atlas_None" = "Old_WT_Parous",
  "old_NP_WT_Aging_atlas_None" = "Old_WT_NP",
  "young_NP_WT_Aging_atlas_None" = "Young_WT_NP",
  "mid_NP_BRCA1_Tumour_Aging_atlas_None" = "BRCA1_Tumour",
  "mid_NP_BRCA1_Aging_atlas_None" = "BRCA1",
  "mid_NP_WT_Aging_atlas_None" = "Mid_WT_NP",
  "young_9.5dG_WT_Aging_atlas_None" = "Young_WT_9.5dG",
  "young_4.5dG_WT_Aging_atlas_None" = "Young_WT_4.5dG",
  "young_14.5dG_WT_Aging_atlas_None" = "Young_WT_14.5dG",
  "mid_NP_PALB2_BRCA2_PALB2_None" = "PALB2",
  "mid_NP_BRCA2_BRCA2_PALB2_None" = "BRCA2",
  "mid_NP_BRCA2_Tumour_BRCA2_PALB2_None" = "BRCA2_Tumour",
  "mid_NP_WT_BRCA2_PALB2_None" = "BRCA2_WT",
  "mid_NP_BRCA1_Senolytics_Seno_ABT737_control" = "ABT737_Control",
  "mid_NP_BRCA1_Senolytics_Seno_ABT737" = "ABT737_Treated",
  "mid_NP_BRCA1_WKAA_First_None" = "BRCA1"
)

process_adata <- function(adata, match_dict) {
  meta <- as.data.frame(colData(adata))

  meta$Age <- suppressWarnings(as.numeric(meta$Age))
  meta$age_group <- ifelse(meta$Age < 20, "young",
                           ifelse(meta$Age < 80, "mid", "old"))

  meta$Genotype_old <- meta$Genotype
  meta$Genotype <- as.character(meta$Genotype)
  meta$Genotype[meta$Genotype %in% c("WKBR", "WKAA")] <- "BRCA1"

  meta$Group_long <- paste(
    meta$age_group,
    meta$Parity,
    meta$Genotype,
    meta$Dataset,
    meta$Treatment,
    sep = "_"
  )
  meta$Group <- unname(match_dict[meta$Group_long])

  unmatched <- unique(meta$Group_long[is.na(meta$Group)])
  if (length(unmatched) > 0) {
    message("Unmatched Group_long values:")
    print(unmatched)
  }

  colData(adata) <- S4Vectors::DataFrame(meta)
  adata
}

prepare_milo <- function(milo) {
  milo <- process_adata(milo, match_dict)

  colData(milo)$ct_level3 <- as.character(colData(milo)$ct_level3)
  colData(milo)$ct_level3[colData(milo)$ct_level3 == "Aged_Il33+"] <- "Il33+"

  keep <- milo$Parity %in% c("NP", "Parous") &
    milo$Genotype %in% c("WT", "BRCA1", "BRCA2", "PALB2") &
    !(milo$ct_level3 %in% c("Doublet", "Tumour", "DDC", "Low_Quality")) &
    !is.na(milo$Age) &
    !is.na(milo$ct_level3) &
    !is.na(milo$Genotype) &
    !is.na(milo$Parity)

  milo <- milo[, keep]
  colData(milo) <- droplevels(colData(milo))

  milo$ct_level3 <- factor(milo$ct_level3)
  milo$Genotype <- factor(milo$Genotype, levels = c("WT", "BRCA1", "BRCA2", "PALB2"))
  milo$Parity <- factor(milo$Parity, levels = c("NP", "Parous"))

  milo$NewParity <- factor(milo$Parity, levels = c("NP", "Parous"))
  milo$Age_group <- factor(ifelse(milo$Age < 70, "Young", "Old"), levels = c("Young", "Old"))

  age_min <- min(milo$Age, na.rm = TRUE)
  age_range <- max(milo$Age, na.rm = TRUE) - age_min
  milo$Age_scaled <- if (age_range == 0) 0 else (milo$Age - age_min) / age_range

  milo$NewGenotype <- factor(milo$Genotype != "WT", levels = c(FALSE, TRUE), labels = c("WT", "Mutant"))
  milo$BRCA1 <- factor(milo$Genotype == "BRCA1", levels = c(FALSE, TRUE), labels = c("Other", "BRCA1"))
  milo$BRCA2 <- factor(milo$Genotype == "BRCA2", levels = c(FALSE, TRUE), labels = c("Other", "BRCA2"))
  milo$PALB2 <- factor(milo$Genotype == "PALB2", levels = c(FALSE, TRUE), labels = c("Other", "PALB2"))
  milo$Treat <- factor(milo$Treatment == "Seno_ABT737", levels = c(FALSE, TRUE), labels = c("Control", "Treated"))

  milo
}

run_milo_pipeline <- function(compartment, input_rds, figure_dir, result_dir, k = 40) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(figure_dir, compartment), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(figure_dir, compartment, "shared"), recursive = TRUE, showWarnings = FALSE)

  milo <- readRDS(input_rds)
  milo <- prepare_milo(milo)

  saveRDS(milo, file.path(result_dir, paste0("processed_", compartment, ".rds")))

  p <- plotReducedDim(
    milo,
    dimred = "UMAP",
    text_by = "ct_level3",
    colour_by = "ct_level3",
    point_size = 1
  ) + guides(fill = "none")

  ggsave(
    filename = file.path(figure_dir, compartment, paste0("ct_level3_", compartment, "_umap_new.png")),
    plot = p,
    height = 6,
    width = 5.15
  )

  result_list_ind <- list()
  result_list_ind[["Parity"]] <- milo_func(milo, compartment, k, "Parity", figure_dir)
  result_list_ind[["Age"]] <- milo_func(milo, compartment, k, "Age", figure_dir)
  result_list_ind[["BRCA1"]] <- milo_func(milo, compartment, k, "BRCA1", figure_dir)
  result_list_ind[["BRCA2"]] <- milo_func(milo, compartment, k, "BRCA2", figure_dir)
  result_list_ind[["PALB2"]] <- milo_func(milo, compartment, k, "PALB2", figure_dir)

  saveRDS(
    result_list_ind,
    file.path(result_dir, paste0("milo_res_", compartment, "_ind.rds"))
  )

  shared_dir <- file.path(figure_dir, compartment, "shared")
  milo_built <- construct_milo(milo, compartment, k, shared_dir)

  result_list_shared <- list()
  result_list_shared[["Parity"]] <- milo_func_shared(milo_built, compartment, k, "Parity", shared_dir)
  result_list_shared[["Age"]] <- milo_func_shared(milo_built, compartment, k, "Age", shared_dir)
  result_list_shared[["BRCA1"]] <- milo_func_shared(milo_built, compartment, k, "BRCA1", shared_dir)
  result_list_shared[["BRCA2"]] <- milo_func_shared(milo_built, compartment, k, "BRCA2", shared_dir)
  result_list_shared[["PALB2"]] <- milo_func_shared(milo_built, compartment, k, "PALB2", shared_dir)

  saveRDS(
    result_list_shared,
    file.path(result_dir, paste0("milo_res_", compartment, "_shared.rds"))
  )

  invisible(list(
    milo = milo,
    milo_built = milo_built,
    individual = result_list_ind,
    shared = result_list_shared
  ))
}

run_milo_pipeline(
  compartment = "epithelial", ##change this for different compartments
  input_rds = "../Figure2/data/epi.rds",
  figure_dir = "../Figure2/figure",
  result_dir = "../Figure2/milo_res",
  k = 40
)