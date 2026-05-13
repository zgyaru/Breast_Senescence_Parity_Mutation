library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(ggbeeswarm)


lym_res <- readRDS('../Figure2/milo_res/milo_res_lym_ind.rds')
epi_res <- readRDS('../Figure2/milo_res/milo_res_epi_ind.rds')
mye_res <- readRDS('../Figure2/milo_res/milo_res_mye_ind.rds')
stro_res <- readRDS('../Figure2/milo_res/milo_res_stro_ind.rds')


comparisons <- names(epi_res)

# Create a named list to store merged results
merged_results <- list()

for (cmp in comparisons) {
  e <- epi_res[[cmp]];  e$compartment <- "epi"
  l <- lym_res[[cmp]];  l$compartment <- "lym"
  m <- mye_res[[cmp]];  m$compartment <- "mye"
  s <- stro_res[[cmp]]; s$compartment <- "stro"
  
  merged_results[[cmp]] <- do.call(rbind, list(e, l, m, s))
}

total_cells <- c("PV3", "PV2", "PV1", "LEC_2", "LEC_1", "VEV", "VEAT","VEA_2", "VEA_1", "VEC", "Fb8", "Fb7", "Fb6", "Fb5", "Fb4", "Fb3", "Fb2", "Fb1", "Tam_1", "Tam_2", "Tam_3",
                "Mast Cell","Neutrophil","pDC","mDC","DC2","DC1","Non_Classical_Monocyte","Classical_Monocyte","Mo3_3","Mo3_2","Mo3_1","Mo2","Mo1", "Proliferating_T","CD8_Trm/NKT",
                "NKT17","NK","DN_NKT","B_cells","ILC2","CD4_Treg","CD4_Th2","CD4_Th1","CD8_Exhausted" ,"CD8_NKT","CD8_CTL","CD8_Trm", "CD8_Tem", "CD8_Tcm","T_naive", 'BMYO4', 'BMYO3', 
                'BMYO2', 'BMYO1', 'LHS4', 'LHS3', 'LHS2', 'LHS1', 'Il33+', 'LASP7', 'LASP6', 'LASP5', 'LASP4', 'LASP3', 'LASP2', 'LASP1')

for (cmp in names(merged_results)) {
  df <- merged_results[[cmp]]
  df$CellType <- factor(df$CellType, levels = total_cells)
  merged_results[[cmp]] <- df[order(df$CellType), ]
}

plotDAbeeswarm_1 <- function(da.res, group.by = NULL, alpha = 0.1, subset.nhoods = NULL, group.levels = NULL) {
  
  # Validate and prepare grouping column
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgotten to run annotateNhoods(x, da.res, ", group.by, ")?")
    }
    if (is.numeric(da.res[[group.by]])) {
      warning(group.by, " is numeric. Consider binning it before plotting.")
    }
    da.res <- mutate(da.res, group_by = da.res[[group.by]])
  } else {
    da.res <- mutate(da.res, group_by = "g1")
  }

  # Set factor levels explicitly using group.levels (e.g. total_cells)
  if (!is.null(group.levels)) {
    da.res <- mutate(da.res, group_by = factor(group_by, levels = group.levels))
  } else {
    da.res <- mutate(da.res, group_by = factor(group_by))  # fallback
  }

  # Subset neighborhoods if specified
  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods, ]
  }

  # Prepare plot
  p <- da.res %>%
    mutate(
      is_signif = ifelse(SpatialFDR < alpha, 1, 0),
      logFC_color = ifelse(is_signif == 1, logFC, NA)
    ) %>%
    ggplot(aes(x = logFC, y = group_by, color = logFC_color)) +
    geom_quasirandom(groupOnY = TRUE) +
    scale_y_discrete(drop = FALSE) +
    scale_color_gradient2() +
    guides(color = "none") +
    xlab("Log Fold Change") +
    ylab(group.by) +
    theme_bw(base_size = 22) +
    theme(strip.text.y = element_text(angle = 0))

  return(p)
}


for (cmp in names(merged_results)){
    p <- plotDAbeeswarm_1(merged_results[[cmp]], group.by = "CellType", group.levels = total_cells)
    height <- 0.5 * length(total_cells)
    ggsave(paste0('../summary_figures/dot/', cmp, '_dotplot.pdf'), plot = p, width = 10, height = height)
}



