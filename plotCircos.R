require(circlize)
require(Matrix)
aml.data <- readRDS('~/Research/aml/aml_v2/data/aml_data_v4.rds')
ids <- aml.data$sample.id
mut_mat <- data.matrix(aml.data[,-1])
rownames(mut_mat) <- ids
mut_mat[mut_mat > 1] <- 1
clust_res <- read.csv('~/Research/aml/aml_v2/results/aml_cluster_v4/ClusterResults_V4.csv')

for(c in 1:4){
  local_mut_mat <- mut_mat[rownames(mut_mat) %in% clust_res$Id[clust_res$Cluster == paste0('Cluster-', c)],]  
  
  ## Extract Co-Occurrence ##
  all.ids <- sort(colnames(local_mut_mat))
  freq <- Matrix(0, ncol=ncol(local_mut_mat), nrow=ncol(local_mut_mat))
  for(i in 1:nrow(local_mut_mat)){
    mut_genes <- sort(colnames(local_mut_mat)[local_mut_mat[i,]==1])
    if(length(mut_genes) < 2)
      next
    mut_pair <- t(combn(mut_genes, 2))
    freq <- freq + sparseMatrix(i=match(mut_pair[,1], table=all.ids),
                                j=match(mut_pair[,2], table=all.ids),
                                x=1,
                                dims = c(length(all.ids), length(all.ids)))
  }
  freq <- freq/nrow(local_mut_mat)
  colnames(freq) <- all.ids
  rownames(freq) <- all.ids
  freq_ft <- reshape2::melt(as.matrix(freq))
  scale_ft <- function(x=NULL){
    return((x-min(x))/(max(x)-min(x)))
  }
  freq_ft$value <- scale_ft(freq_ft$value)
  freq_ft <- subset(freq_ft, freq_ft$value > 0)
  
  cluster_colors <- setNames(c('#c0c000',
                               '#04bf78',
                               '#ff7f00',
                               '#b956d7'), nm = paste0('Cluster-', 1:4))
  
  temp <- col2rgb(as.vector(cluster_colors))
  cluster_idx <- c
  freq_ft_cols <- sapply(freq_ft$value, function(alpha_val){
    rgb(temp[1,cluster_idx]/255, temp[2,cluster_idx]/255, temp[3,cluster_idx]/255, alpha = alpha_val)
  })
  cyto_map <- setNames(c('inv(3)/t(3;3)',
                         't(6;9)',
                         '-5/del(5q)', 
                         '-6/del(6q)', 
                         '-7/(del7q)', 
                         '-9/del(9q)', 
                         'del(12p)', 
                         'del(13q)', 
                         'del(16q)', 
                         '-17/del(17p)', 
                         'del(20q)', 
                         '+8', 
                         '-X', 
                         '-Y', 
                         expression(italic('FLT3')^'ITD'), 
                         expression(italic('FLT3')^'TKD'), 
                         expression(italic('CEBPA')^'Mo'), 
                         expression(italic('CEBPA')^'Bi'), 
                         expression(italic('IDH2')^'R140'),
                         expression(italic('IDH2')^'R172'), 
                         expression(italic('BCOR/L1')),
                         expression(italic('ASXL1')),
                         expression(italic('CBL')),
                         expression(italic('DNMT3A')),
                         expression(italic('ETV6')),
                         expression(italic('EZH2')),
                         expression(italic('GATA2')),
                         expression(italic('IDH1')),
                         expression(italic('KIT')),
                         expression(italic('KRAS')),
                         expression(italic('NPM1')),
                         expression(italic('NRAS')),
                         expression(italic('RUNX1')),
                         expression(italic('SF3B1')),
                         expression(italic('SRSF2')),
                         expression(italic('TET2')),
                         expression(italic('TP53')),
                         expression(italic('U2AF1')),
                         expression(italic('WT1')),
                         expression(italic('ZRSR2'))),
                       c('inv3',
                         'tr6',
                         'del5q', 
                         'del6q', 
                         'del7q', 
                         'del9q', 
                         'del12p', 
                         'del13q', 
                         "del16q", 
                         "del17p", 
                         "del20q",
                         "trisomy8",
                         "X",
                         "Y", 
                         'FLT3.ITD', 
                         'FLT3.TKD', 
                         'CEBPA.Mono', 
                         'CEBPA.Bi', 
                         'IDH2.140', 
                         'IDH2.172', 
                         'BCOR.L1',
                         'ASXL1',
                         'CBL',
                         'DNMT3A',
                         'ETV6',
                         'EZH2',
                         'GATA2',
                         'IDH1',
                         'KIT',
                         'KRAS',
                         'NPM1',
                         'NRAS',
                         'RUNX1',
                         'SF3B1',
                         'SRSF2',
                         'TET2',
                         'TP53',
                         'U2AF1',
                         'WT1',
                         'ZRSR2'))
  
  require(Cairo)
  CairoPDF(file=sprintf('~/Research/aml/aml_v2/plots/CircosPlot_Cluster-%d.pdf', cluster_idx), width = 8, height = 8)
  circos.par("track.height" = 0.05, 
             "track.margin" = c(0.25, 0.25))
  circos.initialize(factors = all.ids, xlim=c(0,1))
  circos.track(factors = all.ids,
               ylim=c(0,1),
               bg.col = cluster_colors[c],
               bg.border = 'black',
               panel.fun = function(x,y){
                 circos.text(CELL_META$xcenter, 
                             CELL_META$cell.ylim[2] + mm_y(1.75), ifelse(CELL_META$sector.index %in% names(cyto_map),
                                                                         cyto_map[which(names(cyto_map) == CELL_META$sector.index)],
                                                                         CELL_META$sector.index),
                             cex = 1,
                             facing = 'clockwise',
                             niceFacing = TRUE,
                             adj = c(0, 0.5))
               })
  for(i in 1:nrow(freq_ft)){
    circos.link(sector.index1 = freq_ft$Var1[i], c(0,1), sector.index2 = freq_ft$Var2[i], c(0,1), col = freq_ft_cols[i], rou = 0.675)  
  }
  dev.off()
  circos.clear()
}


