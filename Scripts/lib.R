## HELPER FUNCTIONS ##

## libraries ##
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
#suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(patchwork))
#suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(infotheo))
suppressPackageStartupMessages(library(viridis))

## 10x matrix fun:
#suppressPackageStartupMessages(library(rhdf5))
#suppressPackageStartupMessages(library(readr))
#suppressPackageStartupMessages(library(Matrix))
#suppressPackageStartupMessages(library(tidyr))


## MATURATION GENES ##

## create smooth expression matrix for seurat object
## this only considers highly variable genes (which need to be the genes in the scaled data)
add_smooth_expression <- function(s) {
  ## smoothed expression should be applied to scaled data, doing otherwise throws an error:
  scaled_mtx <- GetAssayData(object = s, slot = "scale.data")
  pt_vec <- as.numeric(s$Pseudotime)
  pt_pred <- seq(min(pt_vec), max(pt_vec), length.out=100)
  n <- length(pt_vec)
  nd <- data.frame(x = pt_pred)
  worker <- function(i) {
    fit <- loess(scaled_mtx[i, ] ~ pt_vec, span = 0.5, degreee = 2, control=loess.control(surface="interpolate", statistics = "none"))
    predict(fit, se = FALSE)
  }
  tmp <- lapply(1:nrow(scaled_mtx), worker)
  expr_fit <- t(matrix(unlist(tmp), n))
  rownames(expr_fit) <- rownames(scaled_mtx)
  colnames(expr_fit) <- colnames(scaled_mtx)
  
  s[["smoothData"]] <- CreateAssayObject(data = expr_fit)
  return(s)
}


plot_cell_cycle_score <- function(s, th = NULL, main = "") {
  fit <- loess(abs(s$CC.Difference) ~ as.numeric(s$Pseudotime), span = 0.2, degree = 2)
  s$CC_phase_fit <- fit$fitted
  plot_df <- data.frame(pseudotime = rep(as.numeric(s$Pseudotime),2),
                        CC_score = c(s$CC.Difference, s$CC_phase_fit),
                        type = rep(c("cc_diff", "cc_fit"), each = length(s$Pseudotime)))
  if(is.null(th)) {
    g <- ggplot(plot_df, aes(x = pseudotime, y = CC_score, color = type)) +
      geom_point() +
      ylim(min(s$CC.Difference), max(s$CC.Difference)) +
      ggtitle(main)
  } else {
    g <- ggplot(plot_df, aes(x = pseudotime, y = CC_score, color = type)) +
      geom_point() +
      geom_vline(xintercept = th)  +
      ylim(min(s$CC.Difference), max(s$CC.Difference)) +
      ggtitle(main)
  }
  return(g)
}

only_keep_firsts <- function(vec) {
  firsts <- c(); l <- c()
  for(el in vec) {
    if(el %in% firsts) {
      l <- c(l, FALSE)
    } else {
      firsts <- c(firsts, el)
      l <- c(l, TRUE)
    }
  }
  return(l)
}

get_high_genes_per_bin <- function(s_traj, nbins = 10) {
  ## get scaled and smoothed data
  #mtx <- s_traj@assays$smoothData@data
  mtx <- s_traj@assays$RNA@scale.data
  ## order stuff according to pseudotime
  pt_vec <- s_traj$Pseudotime; names(pt_vec) <- colnames(s_traj)
  ordered_pt_vec <- pt_vec[order(pt_vec, decreasing = FALSE)]
  ordered_mtx <- mtx[, names(ordered_pt_vec)]
  ## prepare ordered matrix:
  bin_starts <- sapply(seq(1,ncol(ordered_mtx), length.out = nbins + 1), floor)
  bin_size <- ceiling(ncol(ordered_mtx)/nbins)
  ## this is a hard parameter:
  expr_mean <- 0.5
  ## calculate genes that are specific for each bin:
  res <- lapply(seq(1, nbins), function(idx) {
    r <- apply(ordered_mtx, 1, function(row) {
      mean(row[bin_starts[idx]:bin_starts[idx+1]]) - mean(row)
    })
    gene_abundance <- apply(ordered_mtx, 1, function(row) {
      sum(row[bin_starts[idx]:bin_starts[idx+1]] > expr_mean) / bin_size - sum(row > expr_mean)/length(row)
    })
    ## fit normal distribution to r:
    fit <- MASS::fitdistr(r, "normal")
    th <- fit$estimate[1] + 2*fit$estimate[2]
    r_names <- names(r)[r > th]
    ## also fit normal distribution to gene_abundance:
    ga_fit <- MASS::fitdistr(gene_abundance, "normal")
    ga_th <- ga_fit$estimate[1] + 2*ga_fit$estimate[2]
    ga_names <- names(gene_abundance)[gene_abundance > ga_th]
    ## keep genes that are in r_names and ga_names:
    res_names <- r_names[r_names %in% ga_names]
    
    return(res_names)
  })
  res <- res[sapply(res, function(el){length(el) > 0})]
  res <- unlist(res)
  res <- res[only_keep_firsts(res)]
  ## order according to maximum expression:
  res <- res[order(sapply(res, function(gene) {
    which.max(ordered_mtx[gene, ])
  }), decreasing = FALSE)]
  return(res)
}

order_genes_according_to_pseudotime <- function(genes, s) {
  mtx <- s@assays$smoothData@data
  #mtx <- s@assays$RNA@scale.data
  ## order stuff according to pseudotime
  pt_vec <- s$Pseudotime; names(pt_vec) <- colnames(s)
  ordered_pt_vec <- pt_vec[order(pt_vec, decreasing = FALSE)]
  ordered_mtx <- mtx[, names(ordered_pt_vec)]
  res <- genes[order(sapply(genes, function(gene) {
    which.max(ordered_mtx[gene, ])
  }), decreasing = FALSE)]
  return(res)
}


get_maturation_genes <- function(s, traj_list = c(),nbins = 10, num_trajectories = NULL) {
  if(length(traj_list) > 0) {
    if(is.null(num_trajectories)) {
      num_trajectories <- length(traj_list)
    }
    ## get maturation genes for each trajectory:
    maturation_gene_list <- lapply(names(traj_list), function(traj_name) {
      # subset for trajectory:
      traj_cellIDs <- traj_list[[traj_name]]
      s$in_traj <- as.numeric(colnames(s) %in% traj_cellIDs)
      s_traj <- subset(s, subset = in_traj == 1)
      # get maturation genes:
      traj_maturation_genes <- get_high_genes_per_bin(s_traj, nbins = nbins)
    })
    putative_maturation_genes <- levels(as.factor(unlist(maturation_gene_list)))
    ## only consider genes that are present in all trajectories
    maturation_genes <- putative_maturation_genes[sapply(putative_maturation_genes, function(gene) {
      sum(sapply(maturation_gene_list, function(traj_genes) {gene %in% traj_genes})) >= num_trajectories
    })]
    maturation_genes <- order_genes_according_to_pseudotime(maturation_genes, s)
  } else {
    ## treat whole dataset as one trajectory:
    maturation_genes <- get_high_genes_per_bin(s, nbins = nbins)
  }
  return(maturation_genes)
}


## main plotting function:
plot_mtx_with_ggplot <- function(mat_genes, s, mpm_thres = -1, legend = TRUE, ysize = 5, sort = FALSE, use_smoothed = TRUE) {
  if(sort) {
    mat_genes <- order_genes_according_to_pseudotime(mat_genes, s)
  }
  if(use_smoothed == TRUE) {
    mtx <- s@assays$smoothData@data
  } else {
    mtx <- s@assays$RNA@scale.data
  }
  
  mtx <- mtx[mat_genes, ]
  pt_vec<- as.numeric(s$Pseudotime); names(pt_vec) <- colnames(s)
  
  vmm <- expand.grid(genes = rownames(mtx), cellIDs = colnames(mtx))
  vmm$pseudotime <- rep(pt_vec, each = nrow(mtx))
  vmm$pseudotime <- sapply(vmm$pseudotime, function(el) {sprintf("%0.2f", el)})
  vmm$expression <- as.numeric(mtx)
  vmm$scaled_expression <- scales::rescale(vmm$expression,to = c(-2,2))
  pt_levels <- levels(as.factor(vmm$pseudotime))
  if(legend) {
    g <- ggplot(data = vmm, aes(x = pseudotime, y = genes, fill = scaled_expression)) +
      geom_tile() +
      scale_fill_viridis(option = "D") +
      #scale_fill_gradient2(high = "blue", mid = "yellow", low = "white") +
      scale_y_discrete(limits=rev, expand=c(0,0)) +
      scale_x_discrete(breaks = vmm$pseudotime[seq(1,length(vmm$pseudotime), length.out = 20)], expand = c(0,0),
                       limits = pt_levels[order(as.numeric(pt_levels))]) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(axis.text.y = element_text(size = ysize))
  } else {
    g <- ggplot(data = vmm, aes(x = pseudotime, y = genes, fill = scaled_expression)) +
      geom_tile(show.legend = FALSE) +
      scale_fill_viridis(option = "D") +
      #scale_fill_gradient2(high = "blue", mid = "yellow", low = "white") +
      scale_y_discrete(limits=rev, expand=c(0,0)) +
      scale_x_discrete(breaks = vmm$pseudotime[seq(1,length(vmm$pseudotime), length.out = 20)], expand=c(0,0),
                       limits = pt_levels[order(as.numeric(pt_levels))]) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(axis.text.y = element_text(size = ysize))
  }
  if(mpm_thres != -1) {
    g <- g + geom_vline(xintercept = sprintf("%0.2f", mpm_thres), color = "red")
  }
  return(g)
}


get_pseudotime_of_max_expression <- function(s, gene) {
  cell_idx <- which.max(s@assays$smoothData@data[gene, ])
  return(s$Pseudotime[cell_idx])
}

make_maturation_plots <- function(dorsal_seurat, ventral_seurat, maturation_gene_list, M_PM_border_list, plot_out, dorsal_sample, ventral_sample) {
  ## get data specific for comparison:
  ventral_M_PM_border <- M_PM_border_list[[ventral_sample]]
  dorsal_M_PM_border <- M_PM_border_list[[dorsal_sample]]
  ventral_mat_genes <- maturation_gene_list[[ventral_sample]]
  dorsal_mat_genes <- maturation_gene_list[[dorsal_sample]]
  
  ## plot maturation genes along pseudotime:
  p1 <- plot_mtx_with_ggplot(ventral_mat_genes, ventral_seurat, mpm_thres = ventral_M_PM_border, legend = T, sort = T)
  p2 <- plot_mtx_with_ggplot(dorsal_mat_genes, dorsal_seurat, mpm_thres = dorsal_M_PM_border, legend = T, sort = T)
  
  ## how many maturation genes are overlapping?
  gene_lists <- list(ventral = ventral_mat_genes, dorsal = dorsal_mat_genes)
  p3 <- ggvenn::ggvenn(gene_lists)
  ## plot overlapping maturation genes:
  overlap_genes <- ventral_mat_genes[ventral_mat_genes %in% dorsal_mat_genes]
  if(length(overlap_genes) > 0) {
    pp1 <- plot_mtx_with_ggplot(overlap_genes, ventral_seurat, mpm_thres = ventral_M_PM_border, legend = F)
    pp2 <- plot_mtx_with_ggplot(overlap_genes, dorsal_seurat, mpm_thres = dorsal_M_PM_border, legend = T)
    p4 <- wrap_plots(pp1,pp2,ncol = 2)
  } else {
    p4 <- ggplot()
  }
  
  ## quantify difference between overlapping genes (in terms of when do their expression peak?)
  if(length(overlap_genes) > 0) {
    gene_pt_diff <- sapply(overlap_genes, function(gene) {
      vp <- get_pseudotime_of_max_expression(ventral_seurat, gene)
      dp <- get_pseudotime_of_max_expression(dorsal_seurat, gene)
      return(vp - dp)
    })
    diff_df <- data.frame("gene" = overlap_genes, "diff" = gene_pt_diff)
    p5 <- ggplot(diff_df, aes(diff, gene)) +
      geom_point() +
      geom_vline(xintercept = 0.0, color = "red") +
      #xlim(-1,1) +
      theme(axis.text.y = element_text(size = 5))
  } else {
    p5 <- ggplot()
  }
  
  ## MISSING: eucledian distance between rows/cols ##
  
  # plot results to pdf:
  pdf(plot_out, width = 6, height = 6)
  plot(p1); plot(p2); plot(p3)
  plot(p4); plot(p5)
  dev.off()
}


####################### random utility functions #########################

## convert human gene names to mouse gene names:
convert_to_mouse <- function(vec) {
  sapply(vec, function(el) {
    el <- tolower(el)
    substr(el,1,1) <- toupper(substr(el,1,1))
    return(el)
  })
}

## bed file is 0-indexed.. usually we don't want that, so I add +1 here:
archr_peak_file_to_GRanges <- function(peak_file) {
  bed_df <- read.table(peak_file, sep = "\t", h=T)
  gr <- GRanges(seqnames = bed_df$seqnames,
                ranges = IRanges(start = bed_df$starts + 1, end = bed_df$ends + 1, names = bed_df$names),
                scores = bed_df$scores,
                groupScoreQuantiles = bed_df$groupScoreQuantiles,
                distToGeneStart = bed_df$distToGeneStart,
                nearestGene = bed_df$nearestGene,
                peakType = bed_df$peakType)
  return(gr)
}

encode_bed_to_GRanges <- function(bed_file) {
  bed_df <- read.table(bed_file, sep = "\t", h=F)
  gr <- GRanges(seqnames = bed_df[,1],
                ranges = IRanges(start = bed_df[,2] + 1, end = bed_df[,3] + 1, names = bed_df[,4]),
                strand = Rle(strand(c("*")), c(nrow(bed_df))),
                scores = bed_df[,5],
                signalValue = bed_df[,7],
                pValue = bed_df[,8],
                qValue = bed_df[,9],
                peak = bed_df[,10]
  )
  return(gr)
}


plot_smooth_gene_expression <- function(s, gene) {
  d <- data.frame("Pseudotime" = s$Pseudotime, "smoothed_expression" = s@assays$smoothData@data[gene, ])
  p <- ggplot(d, aes(Pseudotime, smoothed_expression)) +
    geom_point()
  return(p)
}

## functions for network analysis ##

binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}


get_sub_network <- function(sba, expr_th = 0.5, pct_th = 0.1) {
  eRegulon_AUC_sub <- eRegulon_AUC[rownames(cell_metadata)[cell_metadata$stage_broad_annotation2 == sba], ]
  
  eRegulon_AUC_sub_bin <- sapply(1:ncol(eRegulon_AUC_sub), function(j) {
    eRegulon_AUC_sub[, j] > eRegulon_AUC_thresholds[j, ]
  })
  rownames(eRegulon_AUC_sub_bin) <- rownames(eRegulon_AUC_sub); colnames(eRegulon_AUC_sub_bin) <- colnames(eRegulon_AUC_sub)
  eRegulon_AUC_sub_bin[1:10,1:5]
  
  num_regulons_per_cell <- apply(eRegulon_AUC_sub_bin, 1, sum)
  num_cells_per_regulon <- apply(eRegulon_AUC_sub_bin, 2, sum) / nrow(eRegulon_AUC_sub)
  
  hist(num_regulons_per_cell)
  hist(num_cells_per_regulon)
  
  active_regulons <- names(num_cells_per_regulon)[num_cells_per_regulon > 0.5]
  active_regulons
  
  active_TFs <- sapply(active_regulons, function(el) {strsplit(el, "_")[[1]][1]})
  
  eRegulon_md_sub <- eRegulon_md[eRegulon_md$TF %in% active_TFs, ]
  target_genes <- unique(eRegulon_md_sub$Gene)
  
  ## filter target genes for expression and pct-expressed:
  target_gene_avg_expr <- apply(cfse_seurat@assays$RNA@data[target_genes, cfse_seurat$stage_broad_annotation == sba], 1, mean)
  hist(target_gene_avg_expr, breaks = 30)
  
  target_gene_pct_expr <- apply(cfse_seurat@assays$RNA@data[target_genes, cfse_seurat$stage_broad_annotation == sba], 1, function(gv) {
    sum(gv > expr_th) / length(gv)
  })
  hist(target_gene_pct_expr)
  
  active_target_genes <- target_genes[target_gene_avg_expr > expr_th & target_gene_pct_expr > pct_th]
  
  eRegulon_md_sub_filtered <- eRegulon_md_sub[eRegulon_md_sub$Gene %in% active_target_genes, ]
  ## only keep positive interactions:
  eRegulon_md_sub_filtered <- eRegulon_md_sub_filtered[eRegulon_md_sub_filtered$TF2G_regulation == 1, ]
  return(eRegulon_md_sub_filtered)
}


get_complete_network <- function(eRegulon_md_sub, sba, only_TFs = FALSE) {
  edge_df <- eRegulon_md_sub[, c("TF","Gene","TF2G_importance_x_rho","TF2G_rho","TF2G_regulation")]
  colnames(edge_df) <- c("to","from","TF2G_importance_x_rho","TF2G_rho","TF2G_regulation")
  if(only_TFs) {
    edge_df <- edge_df[edge_df$from %in% mm10_tfs, ]
  }
  vertex_df <- data.frame("name" = unique(c(edge_df$from, edge_df$to)))
  vertex_df$type <- sapply(vertex_df$name, function(el) {
    if(el %in% edge_df$to) {"TF"} else {"gene"}
  })
  vertex_df$avg_expression <- apply(cfse_seurat@assays$RNA@data[vertex_df$name, cfse_seurat$stage_broad_annotation == sba], 1, mean)
  
  dg <- graph_from_data_frame(edge_df, directed = T, vertices = vertex_df)
  
  dg <- igraph::simplify(dg, remove.multiple = T, remove.loops = T, 
                         edge.attr.comb = c(TF2G_importance_x_rho="mean", TF2G_rho="mean", TF2G_regulation="mean"))
  
  E(dg)$width <- E(dg)$TF2G_rho*4
  edge_color_vec <- c("1"= "blue","-1" = "red", "0" = "grey")
  E(dg)$color <- edge_color_vec[as.character(E(dg)$TF2G_regulation)]
  E(dg)$edge.arrow.size <- 1 + E(dg)$TF2G_rho
  
  expr_vec <- V(dg)$avg_expression
  expr_vec <- round(expr_vec, digits = 1)
  expr_vec_levels <- unique(expr_vec)[order(unique(expr_vec))]
  colfunc <- colorRampPalette(c("white", "green"))
  col_vec <- colfunc(length(expr_vec_levels))
  names(col_vec) <- expr_vec_levels
  V(dg)$color <- col_vec[as.character(expr_vec)]
  
  vertex_shape_vec <- c("TF"="circle","gene"="square")
  V(dg)$shape <- vertex_shape_vec[V(dg)$type]
  return(dg)
}

## function for plotting one network per cell - type:
combine_networks <- function(e12_dg, e16_dg, celltype, differential_edge_colors = T) {
  ct1 <- paste0(celltype, "_e12"); ct2 <- paste0(celltype, "_e16")
  e12_df <- as_long_data_frame(e12_dg)
  e12_df$edge_id <- paste0(e12_df$from_name, "_", e12_df$to_name)
  e16_df <- as_long_data_frame(e16_dg)
  e16_df$edge_id <- paste0(e16_df$from_name, "_", e16_df$to_name)
  
  edge_df <- rbind(e12_df[,c("from_name","to_name","TF2G_importance_x_rho","TF2G_rho","TF2G_regulation","width","color","edge.arrow.size")],
                   e16_df[,c("from_name","to_name","TF2G_importance_x_rho","TF2G_rho","TF2G_regulation","width","color","edge.arrow.size")])
  e12_from_vertex_df <- e12_df[, c("from_name","from_type","from_avg_expression","from_color","from_shape")]
  colnames(e12_from_vertex_df) <- c("name","type","avg_expression","color","shape")
  e12_to_vertex_df <- e12_df[,c("to_name","to_type","to_avg_expression","to_color","to_shape")]
  colnames(e12_to_vertex_df) <- c("name","type","avg_expression","color","shape")
  e16_from_vertex_df <- e16_df[, c("from_name","from_type","from_avg_expression","from_color","from_shape")]
  colnames(e16_from_vertex_df) <- c("name","type","avg_expression","color","shape")
  e16_to_vertex_df <- e16_df[,c("to_name","to_type","to_avg_expression","to_color","to_shape")]
  colnames(e16_to_vertex_df) <- c("name","type","avg_expression","color","shape")
  
  vertex_df <- rbind(e12_from_vertex_df, e12_to_vertex_df, e16_from_vertex_df, e16_to_vertex_df)
  vertex_df <- vertex_df[!duplicated(vertex_df$name), ]
  
  vertex_df$expr_diff <- apply(cfse_seurat@assays$RNA@data[vertex_df$name, cfse_seurat$stage_broad_annotation == ct2], 1, mean) -
    apply(cfse_seurat@assays$RNA@data[vertex_df$name, cfse_seurat$stage_broad_annotation == ct1], 1, mean)
  
  merged_dg <- graph_from_data_frame(edge_df, directed = T, vertices = vertex_df)
  merged_dg <- igraph::simplify(merged_dg, edge.attr.comb = list(TF2G_importance_x_rho = "mean", TF2G_rho = "mean", TF2G_regulation = "mean",
                                                                 width = "mean", color = "ignore", edge.arrow.size = "mean"))
  
  
  merged_edge_list <- as_edgelist(merged_dg)
  edge_col_vec <- c()
  for(i in 1:nrow(merged_edge_list)) {
    v1 <-  merged_edge_list[i,1]; v2 <- merged_edge_list[i,2]
    if(!v1 %in% V(e12_dg)$name | !v2 %in% V(e12_dg)$name) {
      cc <- "#5776A7"
    } else if(!v1 %in% V(e16_dg)$name | !v2 %in% V(e16_dg)$name) {
      cc <- "#E16B8F"
    } else {
      if(are.connected(e12_dg,v1,v2) & are.connected(e16_dg,v1,v2)) {
        cc <- "grey"
      } else if(are.connected(e12_dg,v1,v2) & !are.connected(e16_dg,v1,v2)) {
        cc <- "#E16B8F"
      } else if(are.connected(e16_dg, v1,v2) & !are.connected(e12_dg, v1, v2)) {
        cc <- "#5776A7"
      } else {
        "grey"
      }
    }
    edge_col_vec <- c(edge_col_vec, cc)
  }
  table(edge_col_vec)
  
  E(merged_dg)$color <- edge_col_vec 
  
  expr_vec <- V(merged_dg)$expr_diff
  expr_vec <- round(expr_vec, digits = 1)
  expr_vec_levels <- unique(expr_vec)[order(unique(expr_vec))]
  colfunc <- colorRampPalette(c("#E16B8F", "white", "#5776A7"))
  col_vec <- colfunc(length(expr_vec_levels))
  names(col_vec) <- expr_vec_levels
  V(merged_dg)$color <- col_vec[as.character(expr_vec)]
  
  V(merged_dg)$shape <- "circle"
  V(merged_dg)$size <- 5
  
  if(!differential_edge_colors) {E(merged_dg)$color <- "black"}
  return(merged_dg)
}


get_complete_network_per_stage <- function(eRegulon_md_sub, stage, only_TFs = FALSE) {
  edge_df <- eRegulon_md_sub[, c("TF","Gene","TF2G_importance_x_rho","TF2G_rho","TF2G_regulation")]
  colnames(edge_df) <- c("to","from","TF2G_importance_x_rho","TF2G_rho","TF2G_regulation")
  if(only_TFs) {
    edge_df <- edge_df[edge_df$from %in% mm10_tfs, ]
  }
  vertex_df <- data.frame("name" = unique(c(edge_df$from, edge_df$to)))
  vertex_df$type <- sapply(vertex_df$name, function(el) {
    if(el %in% edge_df$to) {"TF"} else {"gene"}
  })
  vertex_df$avg_expression <- apply(cfse_seurat@assays$RNA@data[vertex_df$name, cfse_seurat$Experiment == stage], 1, mean)
  
  dg <- graph_from_data_frame(edge_df, directed = T, vertices = vertex_df)
  
  dg <- igraph::simplify(dg, remove.multiple = T, remove.loops = T, 
                         edge.attr.comb = c(TF2G_importance_x_rho="mean", TF2G_rho="mean", TF2G_regulation="mean"))
  
  E(dg)$width <- E(dg)$TF2G_rho*4
  edge_color_vec <- c("1"= "blue","-1" = "red", "0" = "grey")
  E(dg)$color <- edge_color_vec[as.character(E(dg)$TF2G_regulation)]
  E(dg)$edge.arrow.size <- 1 + E(dg)$TF2G_rho
  
  expr_vec <- V(dg)$avg_expression
  expr_vec <- round(expr_vec, digits = 1)
  expr_vec_levels <- unique(expr_vec)[order(unique(expr_vec))]
  colfunc <- colorRampPalette(c("white", "green"))
  col_vec <- colfunc(length(expr_vec_levels))
  names(col_vec) <- expr_vec_levels
  V(dg)$color <- col_vec[as.character(expr_vec)]
  
  vertex_shape_vec <- c("TF"="circle","gene"="square")
  V(dg)$shape <- vertex_shape_vec[V(dg)$type]
  return(dg)
}

## subset networks for gene lists or TFs ##
subset_network <- function(dg, foi_vec, lev = 2, just_upstream = FALSE) {
  foi_vec <- foi_vec[foi_vec %in% V(dg)$name]
  if(length(foi_vec) == 0) {
    print("Network doesn't contain any of the providwd features")
    return(NA)
  }
  if(just_upstream) {
    dist_from_foi <- distances(dg, mode = "in")[foi_vec,]
  } else {
    dist_from_foi <- distances(dg)[foi_vec,]
  }
  neighbors <- unique(unlist(apply(dist_from_foi, 1, function(r) {names(r)[r <= lev]})))
  tfoi_sub <- subgraph(dg, vids = neighbors)
  return(tfoi_sub)
}
