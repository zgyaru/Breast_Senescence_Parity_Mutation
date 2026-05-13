
construct_milo <- function(adata, celltype, k_value, dir){
    set.seed(123)
      ##Create dir
      if (!dir.exists(paste0(dir, celltype))){
        dir.create(paste0(dir, celltype), showWarnings = FALSE)}
      print('count')
      milo <- Milo(adata)
      names(assays(milo))=c('logcounts2', 'logcounts','raw')
      milo <- buildGraph(milo, k = k_value, d = 20, reduced.dim = "X_scVI", transposed = TRUE)
      milo <- makeNhoods(milo, prop = 0.1, k = k_value, d=20, refined = TRUE, reduced_dims = "X_scVI", refinement_scheme = "graph")
      plotNhoodSizeHist(milo)
      k_value_str = as.character(k_value)
      ggsave(paste0(dir,celltype,'/', celltype, '_', k_value_str, ' hist_new.png'),height=6, width=5.15)
      milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="Mouse_ID")
      return(milo)
}


milo_func <- function(milo, celltype, k_value, comparison, dir){
    set.seed(123)
    if (!dir.exists(paste0(dir, celltype))){
            dir.create(paste0(dir, celltype), showWarnings = FALSE)}
     ##Create dMouse_IDir
    n_celltype <- length(unique(colData(milo)$ct_level3))
    height_milo <- n_celltype * 0.4
    print(paste("Number of cell types detected:", n_celltype))
    set.seed(123)
    print('count')
    if (comparison == "Parity"){
        if (!dir.exists(paste0(dir, celltype, '/', comparison))){
        dir.create(paste0(dir, celltype, '/', comparison), showWarnings = FALSE)}
        print('parity')
        design_df <- data.frame(colData(milo)) %>%
        filter(Parity == 'Parous' | (Genotype == 'WT' & Age < 70)) %>%
        select(Mouse_ID, Age, NewParity)
        design_df <- distinct(design_df)
        design_df <- droplevels(design_df)
        rownames(design_df) = design_df$Mouse_ID
        print('count')
        milo <- Milo(milo)
        milo <- milo[, colData(milo)$Mouse_ID %in% design_df$Mouse_ID]
        colData(milo) <- droplevels(colData(milo))
        mouse_ids <- unique(colData(milo)$Mouse_ID)
        print(mouse_ids)
        names(assays(milo))=c('logcounts2', 'logcounts','raw')
        milo <- buildGraph(milo, k = k_value, d = 20, reduced.dim = "X_scVI", transposed = TRUE)
        milo <- makeNhoods(milo, prop = 0.1, k = k_value, d=20, refined = TRUE, reduced_dims = "X_scVI", refinement_scheme = "graph")
        plotNhoodSizeHist(milo)
        k_value_str = as.character(k_value)
        ggsave(paste0(dir,celltype,'/', celltype, '_', comparison, k_value_str, ' hist_new.png'),height=6, width=5.15)
        milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="Mouse_ID")
        milo.res <- testNhoods(milo, design = ~ NewParity, design.df = design_df, reduced.dim='X_scVI', fdr.weighting="graph-overlap")
    }
    else if (comparison == "Age"){
        if (!dir.exists(paste0(dir, celltype, '/', comparison))){
        dir.create(paste0(dir, celltype, '/', comparison), showWarnings = FALSE)}
        print('age')
        design_df <- data.frame(colData(milo)) %>%
        filter((Genotype == 'WT') & (Parity == "NP")) %>%
        select(Mouse_ID, Age_group_bool)
        design_df <- distinct(design_df)
        design_df <- droplevels(design_df)
        rownames(design_df) = design_df$Mouse_ID
        print('count')
        milo <- Milo(milo)
        milo <- milo[, colData(milo)$Mouse_ID %in% design_df$Mouse_ID]
        mouse_ids <- unique(colData(milo)$Mouse_ID)
        print(mouse_ids)
        colData(milo) <- droplevels(colData(milo))
        names(assays(milo))=c('logcounts2', 'logcounts','raw')
        milo <- buildGraph(milo, k = k_value, d = 20, reduced.dim = "X_scVI", transposed = TRUE)
        milo <- makeNhoods(milo, prop = 0.1, k = k_value, d=20, refined = TRUE, reduced_dims = "X_scVI", refinement_scheme = "graph")
        plotNhoodSizeHist(milo)
        k_value_str = as.character(k_value)
        ggsave(paste0(dir,celltype,'/', celltype, '_', comparison, k_value_str, ' hist_new.png'),height=6, width=5.15)
        milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="Mouse_ID")
        milo.res <- testNhoods(milo, design = ~ Age_group_bool, design.df = design_df, reduced.dim='X_scVI', fdr.weighting="graph-overlap")
    }
    else if (comparison == "BRCA1"){
        if (!dir.exists(paste0(dir, celltype, '/', comparison))){
        dir.create(paste0(dir, celltype, '/', comparison), showWarnings = FALSE)}
        print('BRCA1')
        sen_samples <- c('AA21.4d', 'AA21.4e', 'AA21.4f', 'AA23.4h', 'AA24.4c', 'AA24.4d_R')
        design_df <- data.frame(colData(milo)) %>%
        filter((Parity == "NP") & (Age < 70) & (Genotype %in% c('BRCA1', 'WT')) & !(Mouse_ID %in% sen_samples)) %>%
        select(Mouse_ID, NewGenotype)
        design_df <- distinct(design_df)
        design_df <- droplevels(design_df)
        rownames(design_df) = design_df$Mouse_ID
        print('count')
        milo <- Milo(milo)
        milo <- milo[, colData(milo)$Mouse_ID %in% design_df$Mouse_ID]
        colData(milo) <- droplevels(colData(milo))
        mouse_ids <- unique(colData(milo)$Mouse_ID)
        print(mouse_ids)
        names(assays(milo))=c('logcounts2', 'logcounts','raw')
        milo <- buildGraph(milo, k = k_value, d = 20, reduced.dim = "X_scVI", transposed = TRUE)
        milo <- makeNhoods(milo, prop = 0.1, k = k_value, d=20, refined = TRUE, reduced_dims = "X_scVI", refinement_scheme = "graph")
        plotNhoodSizeHist(milo)
        k_value_str = as.character(k_value)
        ggsave(paste0(dir,celltype,'/', celltype, '_', k_value_str, comparison, ' hist_new.png'),height=6, width=5.15)
        milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="Mouse_ID")
        milo.res <- testNhoods(milo, design = ~ NewGenotype, design.df = design_df, reduced.dim='X_scVI', fdr.weighting="graph-overlap")
    }
    else if (comparison == "BRCA2"){
        if (!dir.exists(paste0(dir, celltype, '/', comparison))){
        dir.create(paste0(dir, celltype, '/', comparison), showWarnings = FALSE)}
        print('BRCA2')
        design_df <- data.frame(colData(milo)) %>%
        filter((Parity == "NP") & (Age < 70) & (Genotype %in% c('BRCA2', 'WT'))) %>%
        select(Mouse_ID, NewGenotype)
        design_df <- distinct(design_df)
        design_df <- droplevels(design_df)
        rownames(design_df) = design_df$Mouse_ID
        print('count')
        milo <- Milo(milo)
        milo <- milo[, colData(milo)$Mouse_ID %in% design_df$Mouse_ID]
        colData(milo) <- droplevels(colData(milo))
        mouse_ids <- unique(colData(milo)$Mouse_ID)
        print(mouse_ids)
        names(assays(milo))=c('logcounts2', 'logcounts','raw')
        milo <- buildGraph(milo, k = k_value, d = 20, reduced.dim = "X_scVI", transposed = TRUE)
        milo <- makeNhoods(milo, prop = 0.1, k = k_value, d=20, refined = TRUE, reduced_dims = "X_scVI", refinement_scheme = "graph")
        plotNhoodSizeHist(milo)
        k_value_str = as.character(k_value)
        ggsave(paste0(dir,celltype,'/', celltype, '_', k_value_str, comparison, ' hist_new.png'),height=6, width=5.15)
        milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="Mouse_ID")
        milo.res <- testNhoods(milo, design = ~  NewGenotype, design.df = design_df, reduced.dim='X_scVI', fdr.weighting="graph-overlap")
    }
    else if (comparison == "PALB2"){
        if (!dir.exists(paste0(dir, celltype, '/', comparison))){
        dir.create(paste0(dir, celltype, '/', comparison), showWarnings = FALSE)}
        print('PALB2')
        design_df <- data.frame(colData(milo)) %>%
        filter((Parity == "NP") & (Age < 70) & (Genotype %in% c('PALB2', 'WT'))) %>%
        select(Mouse_ID, NewGenotype)
        design_df <- distinct(design_df)
        design_df <- droplevels(design_df)
        rownames(design_df) = design_df$Mouse_ID
        print('count')
        milo <- Milo(milo)
        milo <- milo[, colData(milo)$Mouse_ID %in% design_df$Mouse_ID]
        colData(milo) <- droplevels(colData(milo))
        mouse_ids <- unique(colData(milo)$Mouse_ID)
        print(mouse_ids)
        names(assays(milo))=c('logcounts2', 'logcounts','raw')
        milo <- buildGraph(milo, k = k_value, d = 20, reduced.dim = "X_scVI", transposed = TRUE)
        milo <- makeNhoods(milo, prop = 0.1, k = k_value, d=20, refined = TRUE, reduced_dims = "X_scVI", refinement_scheme = "graph")
        plotNhoodSizeHist(milo)
        k_value_str = as.character(k_value)
        ggsave(paste0(dir,celltype,'/', celltype, '_', k_value_str, comparison, ' hist_new.png'),height=6, width=5.15)
        milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="Mouse_ID")
        milo.res <- testNhoods(milo, design = ~ NewGenotype, design.df = design_df, reduced.dim='X_scVI', fdr.weighting="graph-overlap")
    }
    else if (comparison == 'Treat'){
        print('Treat')
        design_df <- data.frame(colData(milo)) %>%
        filter(Dataset == 'Senolytics') %>%
        select(Mouse_ID, Treat)
        design_df <- distinct(design_df)
        design_df <- droplevels(design_df)
        rownames(design_df) = design_df$Mouse_ID
        print(design_df)
        print('count')
        milo <- Milo(milo)
        milo <- milo[, colData(milo)$Mouse_ID %in% design_df$Mouse_ID]
        colData(milo) <- droplevels(colData(milo))
        names(assays(milo))=c('logcounts2', 'logcounts','raw')
        milo <- buildGraph(milo, k = k_value, d = 20, reduced.dim = "X_scVI", transposed = TRUE)
        milo <- makeNhoods(milo, prop = 0.1, k = k_value, d=20, refined = TRUE, reduced_dims = "X_scVI", refinement_scheme = "graph")
        plotNhoodSizeHist(milo)
        k_value_str = as.character(k_value)
        ggsave(paste0(dir,celltype,'/', celltype, '_', k_value_str, comparison, ' hist_new.png'),height=6, width=5.15)
        milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="Mouse_ID")
        milo.res <- testNhoods(milo, design = ~ Treat, design.df = design_df, reduced.dim='X_scVI', fdr.weighting="graph-overlap")
    }
    k_value_str = as.character(k_value)
    ggplot(milo.res, aes(PValue)) + geom_histogram(bins=50)
    ggsave(paste0(dir, celltype,'/',comparison,'/', comparison, '_', celltype, '_', k_value_str, 'p_distribution_new.png'),height=6, width=5.15)
    
    ggplot(milo.res, aes(logFC, -log10(SpatialFDR))) + 
    geom_point() +
    geom_hline(yintercept = 1)
    ggsave(paste0(dir, celltype,'/',comparison,'/', comparison, '_', celltype, '_', k_value_str, ' volcano_new.png'),height=6, width=5.15)

    milo.res$Diff <- sign(milo.res$logFC)
    milo.res$Diff[milo.res$SpatialFDR > 0.1] <- 0
    
    milo <- buildNhoodGraph(milo, overlap=5)
    milo.res <- annotateNhoods(milo, milo.res, 'ct_level3')
    milo.res$CellType <- milo.res$ct_level3  # no filtering

    milo.res_new = milo.res[which(milo.res$CellType != "Doublet"),]
    # milo.res_new = milo.res

    if (celltype == 'epithelial'){
        order <- c('BMYO4', 'BMYO3', 'BMYO2', 'BMYO1', 'LHS4', 'LHS3', 'LHS2', 'LHS1', 'Il33+', 
        'LASP7', 'LASP6', 'LASP5', 'LASP4', 'LASP3', 'LASP2', 'LASP1')
        order <- order[order %in% unique(milo.res_new$CellType)]
        milo_sorted <- milo.res_new[order(match(milo.res_new$CellType, order)), ]
    }
    else if (celltype == 'lymphoid'){
        order <- c("Proliferating_T","CD8_Trm/NKT","NKT17","NK","DN_NKT","B_cells",
        "ILC2","CD4_Treg","CD4_Th2","CD4_Th1","CD8_Exhausted" ,"CD8_NKT","CD8_CTL","CD8_Trm", "CD8_Tem", "CD8_Tcm","T_naive")
        ##filter order to only include cell types present in the data
        order <- order[order %in% unique(milo.res_new$CellType)]
        milo_sorted <- milo.res_new[order(match(milo.res_new$CellType, order)), ]
    }

    else if (celltype == 'myeloid'){
        order <- c("Tam_1", "Tam_2", "Tam_3", "Mast Cell","Neutrophil","pDC","mDC","DC2","DC1","Non_Classical_Monocyte","Classical_Monocyte","Mo3_3","Mo3_2","Mo3_1","Mo2","Mo1")
        ##filter order to only include cell types present in the data
        order <- order[order %in% unique(milo.res_new$CellType)]
        milo_sorted <- milo.res_new[order(match(milo.res_new$CellType, order)), ]
    }

    else if (celltype == 'stromal'){
        order <- c("PV3", "PV2", "PV1", "LEC_2", "LEC_1", "VEV", "VEAT","VEA_2", "VEA_1", "VEC", 
         "Fb8", "Fb7", "Fb6", "Fb5", "Fb4", "Fb3", "Fb2", "Fb1")
        order <- order[order %in% unique(milo.res_new$CellType)]
        milo_sorted <- milo.res_new[order(match(milo.res_new$CellType, order)), ]
    }

    print(unique(milo_sorted$CellType))
    
    plotReducedDim(milo, dimred = "UMAP", text_by = "ct_level3", colour_by=comparison, point_size=1) + guides(fill="none") 
    ggsave(paste0(dir,celltype,'/',comparison,'/', comparison, '_', celltype,'_', k_value_str, ' umap_new.png'),height=6, width=5.15)

    plotReducedDim(milo, dimred = "UMAP", colour_by='ct_level3', point_size=1) + guides(fill="none") 
    ggsave(paste0(dir,celltype,'/',comparison,'/', comparison, '_', celltype,'_', k_value_str, ' umap_new_celltype.png'),height=6, width=5.15)
    
    plotNhoodGraphDA(milo, layout="UMAP", milo_res=milo_sorted, alpha=0.1) 
    ggsave(paste0(dir, celltype,'/',comparison,'/', comparison, '_', celltype, '_', k_value_str,' umap_DA_new.png'),height=6, width=5.15)

    plotNhoodGraph(milo, layout="UMAP", colour_by = "ct_level3")
    ggsave(paste0(dir, celltype,'/',comparison,'/', comparison, '_', celltype, '_', k_value_str,' umap_DA_celltype_new.png'),height=6, width=5.15)
    
    if (!any(milo_sorted$Diff == 1)){
        plotDAbeeswarm1(milo_sorted, group.by = "CellType")
        ggsave(paste0(dir, celltype,'/',comparison,'/', comparison, '_', celltype,'_', k_value_str, ' bee_DA_new.png'),height=height_milo, width=9)
        # plotDAbeeswarm1(milo_sorted, group.by = "CellType_filtered")
        # ggsave(paste0(dir, celltype,'/',comparison,'/', comparison, '_', celltype,'_', k_value_str, ' bee_DA_new_filtered.png'),height=height_milo, width=9)
    }
    else{
        plotDAbeeswarm(milo_sorted, group.by = "CellType")
        ggsave(paste0(dir, celltype,'/',comparison,'/', comparison, '_', celltype,'_', k_value_str, ' bee_DA_new.png'),height=height_milo, width=9)
        # plotDAbeeswarm(milo_sorted, group.by = "CellType_filtered")
        # ggsave(paste0(dir, celltype,'/',comparison,'/', comparison, '_', celltype,'_', k_value_str, ' bee_DA_new_filtered.png'),height=height_milo, width=9)
    }
    return(milo.res)
}

##plot DA beeswarm without coloring by logFC
plotDAbeeswarm1 <- function(da.res, group.by=NULL, alpha=0.1, subset.nhoods=NULL){
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", group.by,")?")
    }
    if (is.numeric(da.res[,group.by])) {
      # stop(group.by, " is a numeric variable. Please bin to use for grouping.")
    }
    da.res <- mutate(da.res, group_by = da.res[,group.by])
  } else {
    da.res <- mutate(da.res, group_by = "g1")
  }

  if (!is.factor(da.res[,"group_by"])) {
    message("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, levels=unique(group_by)))
    # anno_vec <- factor(anno_vec, levels=unique(anno_vec))
  }

  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods,]
  }

  # Get position with ggbeeswarm
  beeswarm_pos <- ggplot_build(
    da.res %>%
      mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
      arrange(group_by) %>%
      ggplot(aes(group_by, logFC)) +
      geom_quasirandom()
  )

  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y

  n_groups <- unique(da.res$group_by) %>% length()

  da.res %>%
    mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
    mutate(logFC_color = ifelse(is_signif==1, logFC, 0)) %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
    mutate(pos_x = pos_x, pos_y=pos_y) %>%
    ggplot(aes(pos_x, pos_y)) +
    geom_point(color="darkgrey") +
    #scale_color_gradient2() +
    #guides(color="none") +
    xlab(group.by) + ylab("Log Fold Change") +
    scale_x_continuous(
      breaks = seq(1,n_groups),
      labels = setNames(levels(da.res$group_by), seq(1,n_groups))
      ) +
    #eom_point() +
    coord_flip() +
    theme_bw(base_size=22) +
    theme(strip.text.y =  element_text(angle=0))}