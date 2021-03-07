library(data.table)
library(dplyr)
library(ggplot2)
library(survminer)
library(survival)
library(ggsci)
library(cluster)
library(ComplexHeatmap)
library(circlize)
require(randomForest)

setwd('~/aml/')

## READ Data ##
aml.data <- readRDS('data/aml_data.rds')
ids <- as.character(aml.data$sample.id)
mut_mat <- data.matrix(aml.data[,-1])
rownames(mut_mat) <- ids
mut_mat[mut_mat > 1] <- 1
clust_res <- read.csv('data/ClusterResults.csv')

## RF ##
mut_mat <- mut_mat[match(clust_res$Id, table=rownames(mut_mat)),]

## Hyperparameter Selection ##
mtry_vals <- sapply(2:10, simplify = TRUE, function(m){
  local_auc <- sapply(1:50, simplify = TRUE, function(i){
    cat(sprintf('Processing: %d-%d\n', m, i))
    rand.idx <- sample(1:nrow(mut_mat), size = ceiling(nrow(mut_mat)*0.8), replace = FALSE)
    rf <- randomForest(x = mut_mat[rand.idx,], 
                       y = as.factor(clust_res$Cluster[rand.idx]),
                       ntree = 1500,
                       mtry = m,
                       replace = TRUE, 
                       importance = FALSE)
    local.pred <- predict(rf, newdata = mut_mat[-rand.idx,])
    local.pred <- table(local.pred, clust_res$Cluster[-rand.idx])
    return(sum(diag(local.pred))/sum(local.pred))
  })
})
write.csv(mtry_vals, file = 'data/Hyperparameter.csv', col.names = FALSE, row.names = FALSE, append = FALSE, quote = FALSE)

library(Cairo)
CairoPDF(file = 'plots/MLL_HyperparameterSelection.pdf', width = 6, height = 8)
boxplot(mtry_vals, main = 'Hyperparameter Selection', xlab = 'Number of Variables', ylab = 'Validation AUC', xaxt='n')
axis(side = 1, at = seq(1, 10, 2), labels = paste0('M-', seq(2,10,2)))
dev.off()

# Importance
m <- 3
imp.res <- sapply(1:100, simplify = FALSE, function(i){
  rf <- randomForest(x = mut_mat, 
                     y = as.factor(clust_res$Cluster),
                     ntree = 1500,
                     mtry = m,
                     replace = TRUE, 
                     importance = TRUE)
  return(rf$importance)
})
saveRDS(imp.res, file = 'results/imp_res.rds')
imp.res.c <- do.call('cbind', sapply(imp.res, simplify = FALSE, function(x){
  x[,3]
}))

temp <- reshape2::melt(t(imp.res.c))
temp$Var2 <- factor(temp$Var2,
                    levels = rev(names(sort(rowMeans(imp.res.c)))))
p <- ggplot(temp, aes(x = Var2, y = value)) +
  geom_boxplot() +
  theme_classic() +
  ylab('Mean Decrease Accuracy') + xlab('') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_x_discrete(labels = c("IDH1" = expression(italic("IDH1")),
                              "IDH2" = expression(italic("IDH2")),
                              "TET2" = expression(italic("TET2")),
                              "NRAS" = expression(italic("NRAS")),
                              "CEBPA" = expression(italic("CEBPA")),
                              "WT1" = expression(italic("WT1")),
                              "NPM1" = expression(italic("NPM1")),
                              "DNMT3A" = expression(italic("DNMT3A")),
                              "FLT3(TKD or ITD)" = expression(italic("FLT3 (TKD/ITD)")),
                              "abnormal.cyto" = "Abnormal Cytogenetics",
                              "SRSF2" = expression(italic("SRSF2")),
                              "ASXL1" = expression(italic("ASXL1")),
                              "trisomy8" = "Trisomy 8",
                              "del7q" = "-7/del(7q)",
                              "RUNX1" = expression(italic("RUNX1")),
                              "del5q" = "-5/del(5q)",
                              "other.cyto" = "Other Abnormal Cytogenetics",
                              "TP53" = expression(italic("TP53")),
                              "KRAS" = expression(italic("KRAS")),
                              "BCOR" = expression(italic("BCOR")),
                              "SF3B1" = expression(italic("SF3B1")),
                              "U2AF1" = expression(italic("U2AF1")),
                              "EZH2" = expression(italic("EZH2")),
                              "CBL" = expression(italic("CBL")),
                              "KIT" = expression(italic("KIT")),
                              "ZRSR2" = expression(italic("ZRSR2")),
                              "Y" = "-Y",
                              "del17p" = "del(17p)",
                              "del9p" = "del(9p)",
                              "del20q" = "del(20q)",
                              "del12p" = "del(12p)",
                              "del13q" = "del(13q)",
                              "del6q" = "-6/del(6q)",
                              "del9q" = "-9/del(9q)",
                              "X" = "-X",
                              "del16q" = "del(16q)"))
ggsave(p, file = 'plots/GlobalImportance.pdf', width = 8, height = 6)
