library(data.table)
library(dplyr)
library(ggplot2)
library(survminer)
library(survival)
library(ggsci)
library(cluster)
library(ComplexHeatmap)
library(circlize)


setwd('~/aml/')

## READ Data ##
aml.data <- readRDS('data/aml_data.rds')
ids <- as.character(aml.data$sample.id)
mut_mat <- data.matrix(aml.data[,-1])
rownames(mut_mat) <- ids
mut_mat[mut_mat > 1] <- 1
freq_mat <- readRDS('data/LCA_CC_Mat.rds')
diag(freq_mat) <- 1


## CC Metrics ##
sil.vals <- numeric(20)
for(i in 1:20){
  message('Processing K: ', i)
  
  ## Cluster ##
  hc.res <- hclust(as.dist(1-freq_mat), method = 'ward.D2')
  clust.assign <- cutree(hc.res, k=i)
  message('Minimum Sample Size: ', min(table(clust.assign)))
  
  ## Silhouette Values ##
  sil.res <- silhouette(clust.assign, dist = as.dist(1-cc_mat))
  if(i == 1){
    sil.vals[i] <- 1  
  }else{
    sil.vals[i] <- summary(sil.res)$avg.width
  }
}

pdf(file = 'plots/SilhouetteValues.pdf')
plot(sil.vals, type = 'o', col = 'black', ylim = c(0,1), ylab = '', xlab = '', xaxt = 'n')
axis(side = 1, at = seq(1, 20, 2), labels = paste0('K-', seq(1, 20, 2)))
legend('topright', 
       legend = c('Silhouette'), 
       col = 'black',
       bty = 'n',
       lty = 1, 
       cex = 0.8)
dev.off()

cluster_colors <- setNames(c('#c0c000',
                             '#04bf78',
                             '#ff7f00',
                             '#b956d7'), nm = paste0('Cluster-', 1:4))

library(Cairo)
clust.assign <- cutree(hc.res, k=4)
CairoPDF('plots/SilhouettePlot.pdf', width = 8, height = 8)
plot(silhouette(clust.assign, 
                dist = as.dist(1-freq_mat)), 
     col = cluster_colors, 
     border=NA, 
     main = 'Silhouette Plot')
dev.off()

clust_res <- data.frame('Id' = names(clust.assign),
                        'Cluster' = paste0('Cluster-', clust.assign))
write.csv(clust_res, file = 'results/ClusterResults.csv', quote = FALSE, row.names = FALSE)