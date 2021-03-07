library(BayesLCA)
library(Matrix)
library(gRbase)

setwd('~/aml/data/')

aml_data <- readRDS('aml_data.rds')
ids <- aml_data$sample.id
mut_mat <- data.matrix(aml_data[,-1])
rownames(mut_mat) <- ids
mut_mat[mut_mat > 1] <- 1

## Separate Training ##
train_idx <- sample(1:nrow(mut_mat), size = nrow(mut_mat)*0.8, replace = FALSE)
train_mat <- mut_mat[train_idx,]
test_mat <- mut_mat[-train_idx,]

## Filter only pAML or sAML for test purposes ##
#train_mat <- mut_mat
#annot_ft <- annot_ft[match(rownames(train_mat), table=annot_ft$Id),]
#idx <- annot_ft$Type %in% c('N-pAML', 'A-pAML')
#train_mat <- train_mat[idx,]

## Filter only CC Cohort ##
#train_mat <- mut_mat
#cohort_ft <- readxl::read_xlsx('cohort_annot.xlsx')
#cc_ids <- paste0('sample.', cohort_ft$ID[cohort_ft$Cohort == 'MLL'])
#train_mat <- train_mat[rownames(train_mat) %in% cc_ids,]

## EM ##
require(parallel)
cl <- makeCluster(20)
clusterExport(cl, 'train_mat')
blca_res_em <- parSapply(cl, 1:10000, simplify = FALSE, function(r){
  require(BayesLCA)
  local_mut_mat <- train_mat[sample(1:nrow(train_mat), size = ceiling(nrow(train_mat)*0.75), replace = FALSE),
                             sample(1:ncol(train_mat), size = ncol(train_mat)*1.0, replace = FALSE)]
  tryCatch(expr={
    local_blca_res <- sapply(1:10, simplify = FALSE, function(k){
      res <- blca(X=local_mut_mat, G=k, method='em', restarts=1, start.vals = 'across')
      return(res)
    })
    best_k <- which.max(sapply(local_blca_res, function(x){x$BIC}))
    best_res <- local_blca_res[[best_k]]
    if(best_k == 1){
      cluster_assign <- setNames(rep(1, nrow(local_mut_mat)), rownames(local_mut_mat))
    }else{
      cluster_assign <- setNames(apply(best_res$Z, 1, which.max), rownames(local_mut_mat))  
    }
    return(cluster_assign)
  }, error=function(cond){
    return(NULL)
  }, warning=function(cond){
    return(NULL)
  })
})
stopCluster(cl)

sorted_ids <- sort(Reduce('union', lapply(blca_res_em, names)))
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

saveRDS(freq_mat, 'LCA_CC_Mat.rds')
