library(infercnv)
library(data.table)


setwd("/mnt/home5/pharmacology/yz2071/scratch/projects/Breast-mouse")



expr_file = './data/visium/infercnv_br_wt_region5_expr.tsv'
anno_file  = './data/visium/infercnv_br_wt_region5_labels.tsv'
outdir = './results/BR/BR_WT_region5_inferCNV_allEpi'


gene_pos  = './data/ref/mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions.txt'



anno <- fread(anno_file, header = FALSE)
colnames(anno) <- c("cell_id", "group")

ref_groups <- "WT"

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = expr_file,
  annotations_file  = anno_file,
  delim             = "\t",
  gene_order_file   = gene_pos,
  ref_group_names   = ref_groups
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,         
  out_dir = outdir,
  cluster_by_groups = TRUE,  
  denoise = TRUE,
  HMM = TRUE,              
  #HMM_type = "i6",
  #no_plot = TRUE,
  #no_prelim_plot = TRUE,
  output_format = 'pdf'
)

saveRDS(infercnv_obj, file = file.path(outdir, "infercnv_obj.rds"))





