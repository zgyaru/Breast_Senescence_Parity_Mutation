library(dplyr)
library(ggplot2)

out_dir <- "/home/rx238/rds/hpc-work/final_code/Figure2/summary_figures/signature"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helpers
# ----------------------------
validate_res_df <- function(df, name = "result") {
  required_cols <- c("CellType", "SpatialFDR", "logFC")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "%s is missing required columns: %s",
      name,
      paste(missing_cols, collapse = ", ")
    ))
  }
}

agg_by_celltype <- function(df, label_col, label_value, use_sig = FALSE, fdr_cutoff = 0.1) {
  validate_res_df(df, label_value)

  df2 <- df %>%
    dplyr::select(CellType, SpatialFDR, logFC)

  if (use_sig) {
    df2 <- df2 %>% dplyr::filter(SpatialFDR < fdr_cutoff)
  }

  df2 %>%
    dplyr::group_by(CellType) %>%
    dplyr::summarise(
      mean = if (all(is.na(logFC))) 0 else mean(logFC, na.rm = TRUE),
      sd   = if (sum(!is.na(logFC)) > 1) sd(logFC, na.rm = TRUE) else 0,
      n    = sum(!is.na(logFC)),
      se   = ifelse(n > 1 & is.finite(sd), sd / sqrt(n), 0),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      !!label_col := label_value,
      In_data = TRUE
    )
}

add_missing_celltypes <- function(df, all_celltypes, label_col, label_value) {
  missing_celltypes <- setdiff(all_celltypes, df$CellType)

  if (length(missing_celltypes) == 0) {
    return(df)
  }

  missing_df <- data.frame(
    CellType = missing_celltypes,
    mean = 0,
    sd = 0,
    n = 0,
    se = 0,
    In_data = FALSE
  )
  missing_df[[label_col]] <- label_value

  dplyr::bind_rows(df, missing_df)
}

check_celltype_order <- function(df_list, celltype_order, label = "") {
  observed <- sort(unique(unlist(lapply(df_list, function(x) x$CellType))))
  missing_in_order <- setdiff(observed, celltype_order)
  missing_in_data <- setdiff(celltype_order, observed)

  if (length(missing_in_order) > 0) {
    message("Cell types in data but not in celltype_order for ", label, ":")
    print(missing_in_order)
  }

  if (length(missing_in_data) > 0) {
    message("Cell types in celltype_order but not in data for ", label, ":")
    print(missing_in_data)
  }
}

build_agg_joint <- function(res_list, label_col, celltype_order, use_sig = FALSE, fdr_cutoff = 0.1) {
  check_celltype_order(res_list, celltype_order, label = label_col)

  agg_list <- lapply(names(res_list), function(lbl) {
    agg_by_celltype(res_list[[lbl]], label_col, lbl, use_sig = use_sig, fdr_cutoff = fdr_cutoff)
  })
  names(agg_list) <- names(res_list)

  agg_list <- lapply(names(agg_list), function(lbl) {
    add_missing_celltypes(agg_list[[lbl]], unique(celltype_order), label_col, lbl)
  })
  names(agg_list) <- names(res_list)

  joint <- dplyr::bind_rows(agg_list)
  joint$CellType <- factor(joint$CellType, levels = celltype_order)
  joint <- joint %>% dplyr::filter(!is.na(CellType))

  joint
}

plot_graph_generic <- function(
    df,
    group_col,
    color_values,
    xlim_vec = c(-8, 8),
    show_y_labels = TRUE,
    show_legend = TRUE,
    legend_title = "Significant Nhood"
) {
  df <- df %>%
    dplyr::mutate(
      sd = ifelse(is.na(sd), 0, sd),
      se = ifelse(is.na(se), 0, se)
    ) %>%
    dplyr::arrange(.data[[group_col]], CellType)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      y = CellType,
      x = mean,
      group = .data[[group_col]],
      color = .data[[group_col]]
    )
  ) +
    ggplot2::geom_point(ggplot2::aes(shape = In_data), size = 2.5) +
    ggplot2::geom_path(alpha = 0.5, linewidth = 0.8) +
    ggplot2::geom_errorbarh(
      data = df %>% dplyr::filter(In_data),
      ggplot2::aes(xmin = mean - se, xmax = mean + se),
      height = 0.2,
      alpha = 0.7,
      linewidth = 0.3
    ) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 4, `TRUE` = 16)) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    ggplot2::coord_cartesian(xlim = xlim_vec, expand = FALSE, clip = "off") +
    ggplot2::scale_x_continuous(breaks = seq(xlim_vec[1], xlim_vec[2], by = 4)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.2, color = "grey80"),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 6, 6, 6),
      axis.title.y = if (show_y_labels) ggplot2::element_text() else ggplot2::element_blank(),
      axis.text.y = if (show_y_labels) ggplot2::element_text() else ggplot2::element_blank(),
      legend.position = if (show_legend) "right" else "none"
    ) +
    ggplot2::labs(x = "", y = "", shape = legend_title, title = "")

  p
}

save_dual_plots <- function(df, group_col, color_values, prefix, xlim_vec = c(-8, 8)) {
  n_celltypes <- length(unique(df$CellType))
  height_plot <- max(4, n_celltypes * 0.3)

  grDevices::cairo_pdf(paste0(prefix, "_label.pdf"), width = 6, height = height_plot)
  print(plot_graph_generic(
    df = df,
    group_col = group_col,
    color_values = color_values,
    xlim_vec = xlim_vec,
    show_y_labels = TRUE,
    show_legend = TRUE
  ))
  dev.off()

  grDevices::cairo_pdf(paste0(prefix, "_nolabel.pdf"), width = 4, height = height_plot)
  print(plot_graph_generic(
    df = df,
    group_col = group_col,
    color_values = color_values,
    xlim_vec = xlim_vec,
    show_y_labels = FALSE,
    show_legend = FALSE,
    legend_title = ""
  ))
  dev.off()
}

graph_geno_all <- function(res, celltype_order, celltype, out_dir,
                           use_sig = FALSE, fdr_cutoff = 0.1,
                           xlim_vec = c(-8, 8)) {
  required <- c("BRCA1", "BRCA2", "PALB2")
  missing <- setdiff(required, names(res))
  if (length(missing) > 0) {
    stop("Missing genotype result(s): ", paste(missing, collapse = ", "))
  }

  da_list <- list(
    "Brca1 LOF" = res[["BRCA1"]],
    "Brca2 LOF" = res[["BRCA2"]],
    "Palb2 LOF" = res[["PALB2"]]
  )

  da_agg_joint <- build_agg_joint(
    res_list = da_list,
    label_col = "Genotype",
    celltype_order = celltype_order,
    use_sig = use_sig,
    fdr_cutoff = fdr_cutoff
  ) %>% dplyr::arrange(Genotype, CellType)

  print(da_agg_joint, n = Inf)

  prefix <- file.path(out_dir, paste0("geno_", celltype))
  save_dual_plots(
    df = da_agg_joint,
    group_col = "Genotype",
    color_values = c(
      "Brca1 LOF" = "#ffb4a2",
      "Brca2 LOF" = "#B5838D",
      "Palb2 LOF" = "#6D6875"
    ),
    prefix = prefix,
    xlim_vec = xlim_vec
  )

  utils::write.csv(da_agg_joint, paste0(prefix, ".csv"), row.names = FALSE)
}

graph_age_all <- function(res, celltype_order, celltype, out_dir,
                          use_sig = FALSE, fdr_cutoff = 0.1,
                          xlim_vec = c(-8, 8)) {
  required <- c("Age", "Parity")
  missing <- setdiff(required, names(res))
  if (length(missing) > 0) {
    stop("Missing age/parity result(s): ", paste(missing, collapse = ", "))
  }

  da_list <- list(
    "Age" = res[["Age"]],
    "Parity" = res[["Parity"]]
  )

  da_agg_joint <- build_agg_joint(
    res_list = da_list,
    label_col = "Group",
    celltype_order = celltype_order,
    use_sig = use_sig,
    fdr_cutoff = fdr_cutoff
  ) %>% dplyr::arrange(Group, CellType)

  print(da_agg_joint, n = Inf)

  prefix <- file.path(out_dir, paste0("age_", celltype))
  save_dual_plots(
    df = da_agg_joint,
    group_col = "Group",
    color_values = c(
      "Age" = "#9e2a2b",
      "Parity" = "#e09f3e"
    ),
    prefix = prefix,
    xlim_vec = xlim_vec
  )

  utils::write.csv(da_agg_joint, paste0(prefix, ".csv"), row.names = FALSE)
}

# ----------------------------
# Inputs
# ----------------------------
res_paths <- list(
  lym  = "../milo_res/milo_res_lym_ind.rds",
  epi  = "../milo_res/milo_res_epi_ind.rds",
  mye  = "../milo_res/milo_res_mye_ind.rds",
  stro = "../milo_res/milo_res_stro_ind.rds"
)

celltype_orders <- list(
  lym = c(
    "Proliferating_T", "CD8_Trm/NKT", "NKT17", "NK", "DN_NKT", "B_cells",
    "ILC2", "CD4_Treg", "CD4_Th2", "CD4_Th1", "CD8_Exhausted", "CD8_NKT",
    "CD8_CTL", "CD8_Trm", "CD8_Tem", "CD8_Tcm", "T_naive"
  ),
  mye = c(
    "Tam_1", "Tam_2", "Tam_3", "Mast Cell", "Neutrophil", "pDC", "mDC", "DC2", "DC1",
    "Non_Classical_Monocyte", "Classical_Monocyte", "Mo3_3", "Mo3_2", "Mo3_1", "Mo2", "Mo1"
  ),
  epi = c(
    "BMYO4", "BMYO3", "BMYO2", "BMYO1", "LHS4", "LHS3", "LHS2", "LHS1", "Il33+",
    "LASP7", "LASP6", "LASP5", "LASP4", "LASP3", "LASP2", "LASP1"
  ),
  stro = c(
    "PV3", "PV2", "PV1", "LEC_2", "LEC_1", "VEV", "VEAT", "VEA_2", "VEA_1", "VEC",
    "Fb8", "Fb7", "Fb6", "Fb5", "Fb4", "Fb3", "Fb2", "Fb1"
  )
)

res_list <- lapply(res_paths, readRDS)

for (nm in names(res_list)) {
  graph_geno_all(
    res = res_list[[nm]],
    celltype_order = celltype_orders[[nm]],
    celltype = paste0(nm, "_all"),
    out_dir = out_dir
  )

  graph_age_all(
    res = res_list[[nm]],
    celltype_order = celltype_orders[[nm]],
    celltype = paste0(nm, "_all"),
    out_dir = out_dir
  )
}
