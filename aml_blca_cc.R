library(BayesLCA)
library(ggplot2)
library(survminer)
library(survival)

aml_data <- readRDS('data/aml_data_v3.rds')
mut_mat <- aml_data$mut_mat
annot_ft <- aml_data$annot

## EM ##
require(parallel)
cl <- makeCluster(20)
clusterExport(cl, 'mut_mat')
blca_res_em <- parSapply(cl, 1:1000, simplify = FALSE, function(r){
  require(BayesLCA)
  local_mut_mat <- mut_mat[sample(1:nrow(mut_mat), size = nrow(mut_mat)*0.75, replace = FALSE),
                           sample(1:ncol(mut_mat), size = ncol(mut_mat)*0.9, replace = FALSE)]
  local_blca_res <- sapply(1:15, simplify = FALSE, function(k){
    res <- blca(X=local_mut_mat, G=k, method='em', restarts=1, start.vals = 'across')
    return(res)
  })
  best_k <- which.max(sapply(local_blca_res, function(x){x$BIC}))
  best_res <- local_blca_res[[best_k]]
  cluster_assign <- setNames(apply(best_res$Z, 1, which.max), rownames(local_mut_mat))
  return(cluster_assign)
})
stopCluster(cl)
saveRDS(blca_res_em, file = 'data/LCA_CC.rds')

blca_res_em <- readRDS('data/LCA_CC.rds')

library(Matrix)
library(gRbase)

sorted_ids <- sort(rownames(mut_mat))
pair_mat <- Matrix(0, ncol = length(sorted_ids), nrow = length(sorted_ids), dimnames = list(sorted_ids, sorted_ids))
count_mat <- Matrix(0, ncol = length(sorted_ids), nrow = length(sorted_ids), dimnames = list(sorted_ids, sorted_ids))

for(r in 1:length(blca_res_em)){
  message(sprintf("Processing Run Index: %d", r))
  local_res <- blca_res_em[[r]]
  local_e <- as.data.frame(t(combnPrim(sort(match(names(local_res), table = sorted_ids)), 2)))
  local_e <- data.frame('index_a' = local_e$V1,
                        'index_b' = local_e$V2,
                        'group_a' = as.numeric(local_res[match(sorted_ids[local_e$V1], table = names(local_res))]),
                        'group_b' = as.numeric(local_res[match(sorted_ids[local_e$V2], table = names(local_res))]))
  count_mat <- count_mat + sparseMatrix(i = local_e$index_a,
                                        j = local_e$index_b,
                                        x = 1,
                                        dims = c(length(sorted_ids), length(sorted_ids)),
                                        dimnames = list(sorted_ids, sorted_ids))
  
  local_e <- subset(local_e, local_e$group_a == local_e$group_b)
  pair_mat <- pair_mat + sparseMatrix(i = local_e$index_a,
                                      j = local_e$index_b,
                                      x = 1,
                                      dims = c(length(sorted_ids), length(sorted_ids)),
                                      dimnames = list(sorted_ids, sorted_ids))
  
}

freq_mat <- pair_mat / (count_mat + 1)
freq_mat <- Matrix::forceSymmetric(freq_mat)
freq_mat <- as.matrix(freq_mat)

saveRDS(freq_mat, 'data/LCA_CC_Mat.rds')


## Generate Table ##
### Read Clusters ##
temp_clust <- read.table('~/aml/results/results_20_01_20/MLL_Removed_CombinedTable_LCA_CC_K4.tsv',
                         header = TRUE,
                         sep = '\t')
clust_assign <- setNames(gsub(pattern = 'Cluster-', replacement = '', temp_clust$Cluster),
                         temp_clust$Id)
annot_idx <- match(names(clust_assign), table = rownames(aml_data$annot))
combined_data <- data.frame('Id' = names(clust_assign),
                            'Cluster' = paste0('Cluster-', clust_assign),
                            'Type' = aml_data$annot$Type[annot_idx],
                            'OS' = aml_data$annot$OS[annot_idx],
                            'Event' = aml_data$annot$Event[annot_idx],
                            'Age' = aml_data$annot$Age[annot_idx],
                            'Sex'= aml_data$annot$Sex[annot_idx],
                            'WBC' = aml_data$annot$WBC[annot_idx],
                            'Leukopenia' = aml_data$annot$Leukopenia[annot_idx],
                            'HB' = aml_data$annot$HB[annot_idx],
                            'Anemia' = aml_data$annot$Anemia[annot_idx],
                            'Platelet' = aml_data$annot$Platelet[annot_idx],
                            'Thrombocytopenia' = aml_data$annot$Thrombocytopenia[annot_idx],
                            'BM.Blasts' = aml_data$annot$BM_Blasts[annot_idx],
                            stringsAsFactors = FALSE)
## Summarize ##
## Discrete ##
require(boot)
frq_boot <- function(data, indices){
  d <- na.omit(data[indices])
  uniq_vals <- unique(d)
  if(uniq_vals[1] == 'M' || uniq_vals[1] == 'F'){
    return(sum(d == 'F') / length(d))
  }else{
    return(sum(d == 1) / length(d))
  }
}

med_boot <- function(data, indices){
  d <- na.omit(data[indices])
  return(median(d))
}


disc_tab <- sapply(paste0('Cluster-', 1:4), simplify = FALSE, function(g){
  local_disc_tab <- sapply(c('Sex', 'Leukopenia', 'Anemia', 'Thrombocytopenia'), simplify = FALSE, function(v){
    local_data <- droplevels(subset(combined_data, combined_data$Cluster == g))
    local_data <- as.vector(local_data[,which(colnames(local_data) == v)])
    temp <- boot.ci(boot(local_data, frq_boot, R = 500), type = 'perc')
    return(sprintf("%.3f (%.3f - %.3f)", temp$t0, temp$percent[4], temp$percent[5]))
  })
  return(do.call('rbind', local_disc_tab))
})
disc_tab <- do.call('cbind', disc_tab)
colnames(disc_tab) <- paste0('Cluster-', 1:4)


cont_tab <- sapply(paste0('Cluster-', 1:4), simplify = FALSE, function(g){
  local_cont_tab <- sapply(c('Age', 'WBC', 'HB', 'Platelet', 'BM.Blasts'), simplify = FALSE, function(v){
    local_data <- droplevels(subset(combined_data, combined_data$Cluster == g))
    local_data <- as.vector(local_data[,which(colnames(local_data) == v)])
    temp <- boot.ci(boot(local_data, med_boot, R = 500), type = 'perc')
    return(sprintf("%.3f (%.3f - %.3f)", temp$t0, temp$percent[4], temp$percent[5]))
  })
  return(do.call('rbind', local_cont_tab))
})
cont_tab <- do.call('cbind', cont_tab)
colnames(cont_tab) <- paste0('Cluster-', 1:4)

clin_tab <- rbind(disc_tab, cont_tab)

mut_tab <- sapply(paste0('Cluster-', 1:4), simplify = FALSE, function(g){
  local_mut_mat <- mut_mat[match(combined_data$Id[combined_data$Cluster == g],
                                 table = rownames(mut_mat)),]
  local_mut_tab <- sapply(colnames(local_mut_mat), simplify = FALSE, function(v){
    #print(sprintf("Processing: %s - %s", g, v))
    local_mut_mat <- as.vector(local_mut_mat[,which(colnames(local_mut_mat) == v)])
    temp <- boot.ci(boot(local_mut_mat, frq_boot, R = 500), type = 'perc')
    if(is.null(temp))
      return(sprintf("%.3f (NA - NA)", 1))
    return(sprintf("%.3f (%.3f - %.3f)", temp$t0, temp$percent[4], temp$percent[5]))
  })
  local_mut_tab <- do.call('rbind', local_mut_tab)
  return(local_mut_tab)
})
mut_tab <- do.call('cbind', mut_tab)
colnames(mut_tab) <- paste0('Cluster-', 1:4)

clin_tab <- rbind(clin_tab, mut_tab)

write.csv(clin_tab, file = '~/aml/results/results_4_02_20/MLL_Removed_Table_LCA_CC_K4.csv', quote = FALSE, row.names = TRUE)
write.table(combined_data, file = '~/aml/plots/MLL_Removed_CombinedTable_LCA_CC_K4.tsv', col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, quote = FALSE)

