library(ggplot2)
library(survminer)
library(survival)
library(ggsci)
library(dplyr)
library(randomForest)

setwd('~/Research/aml/aml_v2/')

## READ Data ##
aml_data <- readRDS('data/aml_data_v3.rds')
mut_mat <- aml_data$mut_mat
annot_ft <- aml_data$annot


## CC Metrics ##
require(cluster)
cc_mat <- readRDS('results/LCA_CC_Mat.rds')

sil.vals <- numeric(19)
for(i in 1:19){
  message('Processing K: ', i+1)
  
  ## Cluster ##
  hc.res <- hclust(as.dist(1-cc_mat), method = 'ward.D2')
  clust.assign <- cutree(hc.res, k=i+1)
  message('Minimum Sample Size: ', min(table(clust.assign)))
  
  ## Silhouette Values ##
  sil.res <- silhouette(clust.assign, dist = as.dist(1-cc_mat))
  sil.vals[i] <- summary(sil.res)$avg.width
}

pdf(file = '~/Research/aml/aml_v2/plots/SilhouetteValues.pdf')
plot(sil.vals, type = 'o', col = 'black', ylim = c(0,1), ylab = '', xlab = '', xaxt = 'n')
axis(side = 1, at = seq(1, 19, 2), labels = paste0('K-', seq(2, 20, 2)))
legend('topright', 
       legend = c('Silhouette'), 
       col = 'black',
       bty = 'n',
       lty = 1, 
       cex = 0.8)
dev.off()

# cluster_colors <- setNames(c('#ff7f00',
#                              '#c0c000',
#                              '#04bf78',
#                              '#b956d7'), nm = paste0('Cluster-', 1:4))

cluster_colors <- setNames(c('#c0c000',
                             '#04bf78',
                             '#ff7f00',
                             '#b956d7'), nm = paste0('Cluster-', 1:4))

library(Cairo)
clust.assign <- cutree(hc.res, k=4)
CairoPDF('plots/SilhouettePlot.pdf', width = 8, height = 8)
plot(silhouette(clust.assign, 
                dist = as.dist(1-cc_mat)), 
     col = cluster_colors, 
     border=NA, 
     main = 'Silhouette Plot')
dev.off()

## READ Results ##
clust_res <- read.csv('results/MLL_Removed_CombinedTable_LCA_CC_K4.tsv', stringsAsFactors = FALSE, sep = '\t')
# clust_res <- data.frame('Id' = clust_res$id,
#                         'Cluster' = clust_res$Cluster,
#                         'Type' = annot_ft$Type[match(clust_res$id, table = annot_ft$Id)],
#                         'OS' = annot_ft$OS[match(clust_res$id, table = annot_ft$Id)],
#                         'Event' = annot_ft$Event[match(clust_res$id, table = annot_ft$Id)])

# clust_res <- read.table('results/MLL_Removed_CombinedTable_LCA_CC_K4.tsv',
#                         header = TRUE,
#                         sep = '\t',
#                         quote = '',
#                         stringsAsFactors = FALSE)


clust_res$IdOrig <- gsub(pattern='^sample.', replacement = '', clust_res$Id)
write.csv(clust_res, file = '~/Research/aml/aml_v2/results/for_plots/clust_res.csv', row.names = FALSE, quote = FALSE)

## Survival ##
surv.annot <- data.frame('Cluster' = clust_res$Cluster,
                         'SampleId' = clust_res$Id,
                         'OS' = annot_ft$OS[match(clust_res$Id,
                                                  table=rownames(annot_ft))],
                         'Event' = annot_ft$Event[match(clust_res$Id,
                                                        table=rownames(annot_ft))])
surv.annot$IdOrig <- gsub(pattern='^sample.', replacement = '', surv.annot$SampleId)
write.csv(surv.annot, file = '~/Research/aml/aml_v2/results/for_plots/surv_res.csv', row.names = FALSE, quote = FALSE)

fit <- survfit(Surv(OS, Event) ~ Cluster, data = surv.annot)
p <- ggsurvplot(
  fit, 
  data = surv.annot, 
  size = 0.6,
  conf.int = FALSE,
  pval = TRUE,
  pval.coord = c(125, 0.75),
  risk.table = FALSE,    
  palette = cluster_colors,
  ggtheme = theme_bw(),
  legend.labs = paste0('Cluster-', 1:4),
  legend = 'right',
  legend.title='',
  xlab = 'Time (Months)',
  size = 11)

ggsave(print(p, newpage = FALSE),
       filename = 'plots/KaplanMeierCrossValidation.pdf',
       width = 7,
       height = 5,
       device = cairo_pdf)

## Write Survival Stats ##
for(y in c(12,24,36,60)){
  temp <- summary(fit, times=y)
  local_ft <- data.frame('Cluster' = as.character(temp$strata),
                         'Prob' = temp$surv,
                         'low' = temp$lower,
                         'upper' = temp$upper)
  write.csv(local_ft, sprintf('~/Research/aml/aml_v2/results/SurvivalProb_%dmonths.csv', y), col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE)
}

res <- pairwise_survdiff(Surv(OS, Event) ~ Cluster, data = surv.annot)
#res <- -log10(res$p.value)
#res <- round(res, digits = 2)
res <- reshape2::melt(res$p.value)
colnames(res) <- c('GroupA', 'GroupB', 'Value')

res <- rbind(res,
             data.frame('GroupA' = paste0('Cluster-', 1:4),
                        'GroupB' = paste0('Cluster-', 1:4),
                        'Value' = rep(NA, 4)))
res$Tag <- as.character(formatC(res$Value, format='e', digits=2))
res$Tag[res$GroupA == res$GroupB] <- as.character(res$GroupA[res$GroupA == res$GroupB])
res$GroupA <- factor(res$GroupA, levels=paste0('Cluster-', 1:4))
res$GroupB <- factor(res$GroupB, levels=paste0('Cluster-', 1:4))
res$Tag[res$Tag == ' NA'] <- NA

p <- ggplot(res, aes(x=GroupA, y=GroupB, fill=1-Value)) +
  geom_tile() +
  xlab('') +
  ylab('') +
  theme_void() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.direction = 'horizontal',
        legend.position = c(0.3, 0.8)) +
  geom_text(aes(label=Tag), size=6) +
  scale_fill_viridis_c(limits = c(0,1), na.value = 'white') +
  guides(fill=FALSE)
  guides(fill=guide_colourbar(title = 'P.val', 
                              barheight = 1.5, 
                              barwidth = 10, 
                              title.position = 'top', 
                              title.hjust = 0.5,
                              label.theme = element_text(size=12),
                              title.theme = element_text(size=16, face='bold')))

ggsave(p, filename = 'plots/MLL_Removed_Survival_LCA_CC_K4_Pairwise_Updated.pdf', width = 6, height = 6)

## Primary vs Secondary Survivals ##
clust_assign <- setNames(clust_res$Cluster, clust_res$Id)
col_annot <- data.frame('Cluster' = clust_assign,
                        'Id' = names(clust_assign),
                        'Type' = aml_data$annot$Type[match(names(clust_assign), table = rownames(aml_data$annot))],
                        'OS' = aml_data$annot$OS[match(names(clust_assign),table=rownames(aml_data$annot))],
                        'Event' = aml_data$annot$Event[match(names(clust_assign),table=rownames(aml_data$annot))])
col_annot$ClusterNew <- paste(substr(col_annot$Type,3,3),
                              substr(col_annot$Cluster,9,10),
                              sep = '.')
col_annot_v2 <- data.frame('Id' = rownames(col_annot),
                           'IdOrig' = gsub(pattern='^sample.', replacement = '', rownames(col_annot)),
                           'Type' = ifelse(col_annot$Type %in% c('A-pAML', 'N-pAML'), 'pAML', 'sAML'),
                           'OS' = col_annot$OS,
                           'Event' = col_annot$Event,
                           'Cluster' = col_annot$Cluster)
for(c in paste0('Cluster-', 1:4)){
  local_col_annot <- droplevels(subset(col_annot_v2, col_annot_v2$Cluster == c))
  write.csv(local_col_annot, sprintf('~/Research/aml/aml_v2/results/for_plots/%s.csv', c), row.names = FALSE, quote = FALSE)
}

local_col_annot <- droplevels(subset(col_annot, col_annot$Cluster == 'Cluster-1'))
local_col_annot$ClusterNew <- factor(local_col_annot$ClusterNew, levels=c('p.1', 's.1'))

fit <- survfit(Surv(OS, Event) ~ ClusterNew, data = local_col_annot)
c.1 <- ggsurvplot(
  fit, 
  data = local_col_annot, 
  size = 0.6,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = FALSE,    
  palette = c('red', 'blue'),
  ggtheme = theme_bw(),
  legend.title='',
  legend.labs=c('Primary AML', 'Secondary AML'), 
  xlab = 'Time (Months)')

ggsave(print(c.1, newpage = FALSE),
       filename = 'plots/SurvivalSeparate_Cluster1.pdf',
       width = 8,
       height = 8)

local_col_annot <- droplevels(subset(col_annot, col_annot$Cluster == 'Cluster-2'))
local_col_annot$ClusterNew <- factor(local_col_annot$ClusterNew, levels=c('p.2', 's.2'))

fit <- survfit(Surv(OS, Event) ~ ClusterNew, data = local_col_annot)
c.2 <- ggsurvplot(
  fit, 
  data = local_col_annot, 
  size = 0.6,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = FALSE,    
  palette = c('red', 'blue'),
  ggtheme = theme_bw(),
  legend.title='',
  legend.labs=c('Primary AML', 'Secondary AML'),
  xlab = 'Time (Months)')

ggsave(print(c.2, newpage = FALSE),
       filename = 'plots/SurvivalSeparate_Cluster2.pdf',
       width = 8,
       height = 8)


local_col_annot <- droplevels(subset(col_annot, col_annot$Cluster == 'Cluster-3'))
local_col_annot$ClusterNew <- factor(local_col_annot$ClusterNew, levels=c('p.3', 's.3'))

fit <- survfit(Surv(OS, Event) ~ ClusterNew, data = local_col_annot)
c.3 <- ggsurvplot(
  fit, 
  data = local_col_annot, 
  size = 0.6,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = FALSE,    
  palette = c('red', 'blue'),
  ggtheme = theme_bw(),
  legend.title='',
  legend.labs=c('Primary AML', 'Secondary AML'),
  xlab = 'Time (Months)')

ggsave(print(c.3, newpage = FALSE),
       filename = 'plots/SurvivalSeparate_Cluster3.pdf',
       width = 8,
       height = 8)


local_col_annot <- droplevels(subset(col_annot, col_annot$Cluster == 'Cluster-4'))
local_col_annot$ClusterNew <- factor(local_col_annot$ClusterNew, levels=c('p.4', 's.4'))

fit <- survfit(Surv(OS, Event) ~ ClusterNew, data = local_col_annot)
c.4 <- ggsurvplot(
  fit, 
  data = local_col_annot, 
  size = 0.6,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = FALSE,    
  palette = c('red', 'blue'),
  ggtheme = theme_bw(),
  legend.title='',
  legend.labs=c('Primary AML', 'Secondary AML'),
  xlab = 'Time (Months)')

ggsave(print(c.4, newpage = FALSE),
       filename = 'plots/SurvivalSeparate_Cluster4.pdf',
       width = 8,
       height = 8)


res <- pairwise_survdiff(Surv(OS, Event) ~ ClusterNew, data = col_annot)
#res <- -log10(res$p.value)
#res <- round(res, digits = 2)
res <- reshape2::melt(res$p.value)
colnames(res) <- c('GroupA', 'GroupB', 'Value')
res <- rbind(res,
             data.frame('GroupA' = c(paste0('p.', 1:4), paste0('s.', 1:4)),
                        'GroupB' = c(paste0('p.', 1:4), paste0('s.', 1:4)),
                        'Value' = rep(NA, 8)))
res$Tag <- as.character(formatC(res$Value, format='e', digits=2))
res$Tag[res$GroupA == res$GroupB] <- as.character(res$GroupA[res$GroupA == res$GroupB])
res$GroupA <- factor(res$GroupA, levels = c(paste0('p.', 1:4), paste0('s.', 1:4)))
res$GroupB <- factor(res$GroupB, levels = c(paste0('p.', 1:4), paste0('s.', 1:4)))
res$Tag[res$Tag == ' NA'] <- NA

res_p <- subset(res, grepl(pattern = '^p', res$GroupA) & grepl(pattern = '^p', res$GroupB))
res_s <- subset(res, grepl(pattern = '^s', res$GroupA) & grepl(pattern = '^s', res$GroupB))
p <- ggplot(res_s, aes(x=GroupA, y=GroupB, fill=1-Value)) +
  geom_tile() +
  xlab('') +
  ylab('') +
  theme_void() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.direction = 'horizontal',
        legend.position = c(0.3, 0.8)) +
  geom_text(aes(label=Tag), size=6) +
  scale_fill_viridis_c(limits = c(0,1), na.value = 'white') +
  guides(fill=FALSE)
  guides(fill=guide_colourbar(title = '-log10(P.val)', 
                              barheight = 1.5, 
                              barwidth = 10, 
                              title.position = 'top', 
                              title.hjust = 0.5,
                              label.theme = element_text(size=12),
                              title.theme = element_text(size=16, face='bold')))

ggsave(p, filename = 'plots/MLL_Removed_Survival_LCA_CC_K4_Separate_Pairwise_ComparisonS.pdf', width = 8, height = 8)

pheatmap(p.val.mat.ft, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         breaks = seq(from = 0, to = max(p.val.mat.ft, na.rm = TRUE) + 1, 0.1),
         color = colorRampPalette(c('white', 'firebrick'))(length(seq(from = 0, to = max(p.val.mat.ft, na.rm = TRUE) + 1, 0.1))),
         na_col = 'white', 
         angle_col = 45,
         display_numbers = temp,
         fontsize = 10,
         legend = FALSE,
         filename = '~/aml/results/results_4_02_20/MLL_Removed_Survival_LCA_CC_K4_Separate_Pairwise.pdf',
         width = 4,
         height = 4, 
         fontsize_row = 10, 
         fontsize_col = 10)
graphics.off()
gc()


## Plot Consensus Matrix ##
freq_mat <- readRDS('results/LCA_CC_Mat.rds')
diag(freq_mat) <- 1

## Heatmap ##
require(ComplexHeatmap)
require(circlize)

hc_res <- hclust(as.dist(1-freq_mat), method = 'ward.D2')
clust_assign <- setNames(clust_res$Cluster, clust_res$Id)
clust_assign <- clust_assign[match(colnames(freq_mat), table=names(clust_assign))]

## Modify cluster order ##
temp.clust.1 <- which(clust_assign == 'Cluster-1')
clust_assign[which(clust_assign == 'Cluster-2')] <- 'Cluster-1'
clust_assign[temp.clust.1] <- 'Cluster-2'
temp.clust.2 <- which(clust_assign == 'Cluster-2')
clust_assign[which(clust_assign == 'Cluster-3')] <- 'Cluster-2'
clust_assign[temp.clust.2] <- 'Cluster-3'

col_fun <- colorRamp2(c(0, 0.5, 1), c("dodgerblue", "white", "firebrick"))
column_ha <- HeatmapAnnotation(Type = ifelse(annot_ft$Type[match(rownames(freq_mat), table = rownames(annot_ft))] %in% c('N-pAML', 'A-pAML'), 'pAML', 'sAML'),
                               Cluster = clust_assign[match(rownames(freq_mat), table = names(clust_assign))],
                               col = list(Cluster = cluster_colors,
                                          Type = setNames(c('red', 'dodgerblue'), c('pAML', 'sAML'))),
                               annotation_legend_param = list(Type = list(grid_height = unit(0.5, "cm"),
                                                                          grid_width = unit(0.2, "cm"),
                                                                          labels_gp = gpar(fontsize = 8),
                                                                          title_gp = gpar(fontsize = 10, fontface = 'bold')),
                                                              Cluster = list(grid_height = unit(0.5, 'cm'),
                                                                             grid_width = unit(0.2, 'cm'),
                                                                             labels_gp = gpar(fontsize = 8),
                                                                             title_gp = gpar(fontsize = 10, fontface = 'bold'))))
ht <- Heatmap(freq_mat,
              col = col_fun,
              cluster_rows = hc_res,
              cluster_columns = hc_res,
              show_row_names = FALSE,
              show_column_names = FALSE,
              show_heatmap_legend = TRUE,
              heatmap_legend_param = list(title = 'Frequency', legend_height = unit(4, 'cm'), title_position = 'leftcenter-rot'),
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              use_raster = TRUE,
              raster_device = 'png',
              top_annotation = column_ha)

library(Cairo)
CairoPDF('plots/LCA_CC_Heatmap_Updated.pdf', width = 5, height = 4)
#png('results/LCA_CC_Heatmap_Updated.png', width = 380, height = 360)
draw(ht)
dev.off()

## RF ##
require(randomForest)
mut_mat <- mut_mat[match(clust_res$Id, table=rownames(mut_mat)),]

## Hyperparameter Selection ##
mtry_vals <- sapply(2:10, simplify = TRUE, function(m){
  local_auc <- sapply(1:50, simplify = TRUE, function(i){
    cat(sprintf('Processing: %d-%d\n', m, i))
    rand.idx <- sample(1:nrow(mut_mat), size = 0.8*nrow(mut_mat), replace = FALSE)
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
write.csv(mtry_vals, file = 'results/hyperparameter_selection_CrossValidation.csv', col.names = FALSE, row.names = FALSE, append = FALSE, quote = FALSE)

library(Cairo)
CairoPDF(file = 'plots/HyperparameterSelectionCrossValidation.pdf', width = 6, height = 8)
boxplot(mtry_vals, main = 'Hyperparameter Selection', xlab = 'Number of Variables', ylab = 'Validation AUC', xaxt='n')
axis(side = 1, at = seq(1, 10, 2), labels = paste0('M-', seq(2,10,2)))
dev.off()
#rf_model <- readRDS(file='data/lca_cc_rf_model.rds')
m <- 6

## Survival Crossvalidation ##
train_mat <- mut_mat[match(clust_res$Id, table = rownames(mut_mat)),]
test_mat <- mut_mat[!(rownames(mut_mat) %in% rownames(train_mat)),]

rf_cv <- randomForest(x = train_mat, 
                      y = as.factor(clust_res$Cluster),
                      ntree = 1500,
                      mtry = m,
                      sampsize = c(50,50,50,50),
                      replace = TRUE, 
                      importance = FALSE)
## Plot Combined ##
library(ggfortify)
library(survival)
library(survminer)

pred_res <- predict(rf_cv, newdata = test_mat)
temp_ft <- rbind(data.frame('id' = clust_res$Id,
                            'Cluster' = clust_res$Cluster,
                            'OS' = clust_res$OS,
                            'Event' = clust_res$Event,
                            'Data' = rep('train', nrow(clust_res)),
                            stringsAsFactors = FALSE),
                 data.frame('id' = names(pred_res),
                            'Cluster' = pred_res,
                            'OS' = annot_ft$OS[match(names(pred_res), table = annot_ft$Id)],
                            'Event' = annot_ft$Event[match(names(pred_res), table = annot_ft$Id)],
                            'Data' = rep('test', length(pred_res)),
                            stringsAsFactors = FALSE))
c_plots <- sapply(paste0('Cluster-', 1:4), simplify = FALSE, function(s){
  local_ft <- subset(temp_ft, temp_ft$Cluster == s)
  l_fit <- survfit(Surv(OS, Event) ~ Data, data = local_ft)
  c_title = s
  if(s == 'Cluster-2'){
    c_title <- 'Cluster-1'
  }else if(s == 'Cluster-1'){
    c_title <- 'Cluster-2'
  }
  
  p <- autoplot(l_fit, main=c_title, grid=FALSE) + 
    theme_minimal() + 
    xlim(0, 175) +
    xlab('Time (Months)') +
    ylab('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_discrete(name = 'Data', labels = setNames(c('Test', 'Train'), c('test', 'train'))) +
    guides(fill=FALSE)
  ggsave(p, filename = sprintf('plots/KaplanMeier_CrossValidation_%s.pdf', c_title), device = cairo_pdf, width = 6, height = 4)
  return(p)
})

library(cowplot)
plot_grid(plotlist = c_plots, nrow = 2, ncol = 2)

## Survival Validation ##
require(rms)
require(Hmisc)
aml_data <- readRDS('data/aml_data_v3.rds')
mut_mat <- aml_data$mut_mat
annot_ft <- aml_data$annot
clust_res <- read.csv('results/LCA_CC_CrossValidation_K4.csv', stringsAsFactors = FALSE)
clust_res <- data.frame('Id' = clust_res$id,
                        'Cluster' = clust_res$Cluster,
                        'Type' = annot_ft$Type[match(clust_res$id, table = annot_ft$Id)],
                        'OS' = annot_ft$OS[match(clust_res$id, table = annot_ft$Id)],
                        'Event' = annot_ft$Event[match(clust_res$id, table = annot_ft$Id)])
train_mat <- mut_mat[match(clust_res$Id, table = rownames(mut_mat)),]
test_mat <- mut_mat[!(rownames(mut_mat) %in% rownames(train_mat)),]

rf_cv <- randomForest(x = train_mat, 
                      y = as.factor(clust_res$Cluster),
                      ntree = 1500,
                      mtry = 6,
                      sampsize = c(50,50,50,50),
                      replace = TRUE, 
                      importance = FALSE)
pred_res <- predict(rf_cv, newdata = test_mat)

train_ft <- data.frame('id' = clust_res$Id,
                       'Cluster' = clust_res$Cluster,
                       'OS' = clust_res$OS,
                       'Event' = clust_res$Event,
                       'Data' = rep('train', nrow(clust_res)))
id_tag <- train_ft$Cluster == 'Cluster-1'
train_ft$Cluster[train_ft$Cluster == 'Cluster-2'] <- 'Cluster-1'
train_ft$Cluster[id_tag] <- 'Cluster-2'

test_ft <- data.frame('id' = names(pred_res),
                      'Cluster' = pred_res,
                      'OS' = annot_ft$OS[match(names(pred_res), table = annot_ft$Id)],
                      'Event' = annot_ft$Event[match(names(pred_res), table = annot_ft$Id)],
                      'Data' = rep('test', length(pred_res)))

id_tag <- test_ft$Cluster == 'Cluster-1'
test_ft$Cluster[test_ft$Cluster == 'Cluster-2'] <- 'Cluster-1'
test_ft$Cluster[id_tag] <- 'Cluster-2'

write.csv(train_ft, file = '~/Research/aml/aml_v2/results/Training_Survival.csv', append = FALSE, row.names = FALSE, col.names = TRUE, quote=FALSE)
write.csv(test_ft, file = '~/Research/aml/aml_v2/results/Test_Survival.csv', append = FALSE, row.names = FALSE, col.names = TRUE, quote=FALSE)


fit <- survfit(Surv(OS, Event) ~ Cluster, data = train_ft)
logrank_res <- survdiff(Surv(OS, Event) ~ Cluster, data = train_ft)
p_val_res <- pchisq(150, df=3, lower.tail = FALSE)

p <- autoplot(fit, grid=FALSE, conf.int = FALSE) + 
  theme_minimal() + 
  xlim(0, 175) +
  xlab('Time (Months)') +
  ylab('') +
  labs(color='') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = cluster_colors) +
  guides(fill=FALSE) +
  annotate("text", x = 150, y = 0.75, label = 'p-value < 2e-16', parse=TRUE)
ggsave(p, filename = 'plots/KaplanMeierTrainData.pdf', width = 6, height = 4, device = cairo_pdf)


fit <- survfit(Surv(OS, Event) ~ Cluster, data = test_ft)
logrank_res <- survdiff(Surv(OS, Event) ~ Cluster, data = test_ft)
p_val_res <- pchisq(30.5, df=3, lower.tail = FALSE)

p <- autoplot(fit, grid=FALSE, conf.int = FALSE) + 
  theme_minimal() + 
  xlim(0, 175) +
  xlab('Time (Months)') +
  ylab('') +
  labs(color='') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = cluster_colors) +
  guides(fill=FALSE) +
  annotate("text", x = 150, y = 0.75, label = 'p-value < 2e-06', parse=TRUE)
ggsave(p, filename = 'plots/KaplanMeierTestData.pdf', width = 6, height = 4, device = cairo_pdf)


l_fit <- rms::cph(Surv(OS, Event) ~ Cluster, data=train_ft, x=TRUE, y=TRUE, surv=TRUE)
test_est <- survest(l_fit, newdata = test_ft, times = 175)$surv
rcorr.cens(x=test_est, S=Surv(test_ft$OS, test_ft$Event))


for(s in paste0('Cluster-', 1:4)){
  train_mat_clust <- train_mat[match(clust_res$Id[clust_res$Cluster == s], table = rownames(train_mat)),]
  test_mat_clust <- test_mat[match(names(pred_res)[pred_res == s], table = rownames(test_mat)),]
  
  feat_idx <- apply(train_mat_clust, 2, function(x){min(table(x)) > 0.1*length(x) && length(unique(x)) == 2})
  train_mat_clust <- train_mat_clust[,feat_idx]
  test_mat_clust <- test_mat_clust[,feat_idx]
  
  train_ft <- cbind(data.frame('OS' = annot_ft$OS[match(rownames(train_mat_clust), table = annot_ft$Id)],
                               'Event' = annot_ft$Event[match(rownames(train_mat_clust), table = annot_ft$Id)],
                               stringsAsFactors = FALSE),
                    as.data.frame(train_mat_clust))
  test_ft <- cbind(data.frame('OS' = annot_ft$OS[match(rownames(test_mat_clust), table = annot_ft$Id)],
                              'Event' = annot_ft$Event[match(rownames(test_mat_clust), table = annot_ft$Id)],
                              stringsAsFactors = FALSE),
                   as.data.frame(test_mat_clust))
  ## Model ##
}

require(randomForest)
mut_mat <- mut_mat[match(clust_res$Id, table = rownames(mut_mat)),]
imp.res <- sapply(1:100, simplify = FALSE, function(i){
  rf <- randomForest(x = mut_mat, 
                     y = as.factor(clust_res$Cluster),
                     ntree = 1500,
                     mtry = m,
                     sampsize = c(50,50,50,50),
                     replace = TRUE, 
                     importance = TRUE)
  return(rf$importance)
})
saveRDS(imp.res, file = 'data/imp_res_CrossValidation.rds')
imp.res.c <- do.call('cbind', sapply(imp.res, simplify = FALSE, function(x){
  x[,5]
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
                              "X" = "-X",
                              "del16q" = "del(16q)"))
ggsave(p, file = 'plots/LCA_CC_RandomForest_GlobalImportance.pdf', width = 8, height = 6)

## Write Median ##
require(dplyr)
imp_ft <- temp %>%
  group_by(Var2) %>%
  summarise(Importance=median(value))
write.csv(imp_ft, 
          file = '~/Research/aml/aml_v2/results/GlobalImportance.csv', 
          col.names = TRUE, 
          row.names = FALSE,
          append = FALSE, 
          quote = FALSE)

## Class Level Importance ##
imp.res.c <- do.call('rbind',
                     sapply(imp.res, simplify = FALSE, function(x){
                       x <- x[,1:4]
                       colnames(x) <- paste0('Cluster-', 1:4)
                       temp <- reshape2::melt(x)
                       colnames(temp) <- c('Variable', 'Cluster', 'MeanDecreaseAccuracy')
                       return(temp)
                     }))
imp.res.c.ft <- imp.res.c %>%
  group_by(Cluster, Variable) %>%
  summarise(Importance=median(MeanDecreaseAccuracy))

write.csv(imp.res.c.ft, file = '~/Research/aml/aml_v2/results/ClassImportance.csv', row.names = FALSE, quote = FALSE)

p <- ggplot(imp.res.c, aes(x = Variable, y = MeanDecreaseAccuracy)) +
  geom_boxplot() +
  facet_wrap(~Cluster) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        strip.text = element_text(face='bold', size=12)) +
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
                              "X" = "-X",
                              "del16q" = "del(16q)"))
ggsave(p, file = 'plots/LCA_CC_RandomForest_Importance.pdf', width = 12, height = 10)  

## Cluster Profiles  && Class Level Importance ##
mut_mat_ft <- reshape2::melt(mut_mat)
colnames(mut_mat_ft) <- c('Sample', 'Variable', 'Status')
mut_mat_ft$Cluster <- clust_assign[match(mut_mat_ft$Sample, table = names(clust_assign))]
#mut_mat_ft$Status <- ifelse(mut_mat_ft$Status == 1, 'Mut', 'No-Mut')

require(dplyr)
require(ggplot2)
require(ggsci)
for(i in 1:4){
  c <- paste0('Cluster-', i)
  local_imp_res_c <- droplevels(subset(imp.res.c, imp.res.c$Cluster == c))
  local_mut_mat_ft <- droplevels(subset(mut_mat_ft, mut_mat_ft$Cluster == c))
  
  local_mut_mat_ft <- local_mut_mat_ft %>% 
    group_by(Cluster, Variable) %>% 
    summarise(Freq=sum(Status)/n())
  
  local_imp_res_stat <- local_imp_res_c %>%
    group_by(Cluster, Variable) %>%
    summarise(Med=median(MeanDecreaseAccuracy))
  cluster_color <- pal_jco()(4)[i]
  
  write.csv(local_mut_mat_ft, file = sprintf('~/Research/aml/aml_v2/results/for_plots/%s_mut_profile.csv', c), row.names = FALSE, quote = FALSE)
  write.csv(local_imp_res_stat, file = sprintf('~/Research/aml/aml_v2/results/for_plots/%s_acc_profile.csv', c), row.names = FALSE, quote = FALSE)
  next
  alpha_res <- ifelse(local_imp_res_stat$Med >= 0.02, 1, 0.2)
  
  ## Plot ##
  p_col <- ggplot(local_mut_mat_ft, aes(x = Variable, y = Freq)) +
    geom_col(fill = I(cluster_color), alpha = alpha_res) +
    xlab('') +
    ylab('Frequency') +
    ggtitle('') +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          strip.text = element_text(face='bold', size=12))

  p_imp <- ggplot(local_imp_res_c, aes(x = Variable, y = MeanDecreaseAccuracy)) +
    geom_boxplot() +
    theme_minimal() +
    ggtitle('') +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
          strip.text = element_text(face='bold', size=12)) +
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
                                "X" = "-X",
                                "del16q" = "del(16q)"))
  
  title <- cowplot::ggdraw() +
    cowplot::draw_label(sprintf('%s', c),
               fontface = 'bold',
               x = 0,
               hjust = 0)
  p <- cowplot::plot_grid(title, p_col, p_imp, align='v', ncol = 1, rel_heights = c(0.1, 1, 1))
  cowplot::save_plot(filename = sprintf('~/aml/results/results_4_02_20/MergedImportance_%s.pdf', c),
                     p, base_height = 6, base_width = 6)
  
}


## Pie Chart ##
col_annot <- data.frame('Cluster' = clust_assign,
                        'Id' = names(clust_assign),
                        'Type' = aml_data$annot$Type[match(names(clust_assign), table = rownames(aml_data$annot))],
                        'OS' = aml_data$annot$OS[match(names(clust_assign),table=rownames(aml_data$annot))],
                        'Event' = aml_data$annot$Event[match(names(clust_assign),table=rownames(aml_data$annot))])

p_aml <- col_annot[col_annot$Type == 'N-pAML' | col_annot$Type == 'A-pAML',]
s_aml <- col_annot[col_annot$Type == 'N-sAML' | col_annot$Type == 'A-sAML',]

p_aml <- as.data.frame(table(p_aml$Cluster)/nrow(p_aml))
s_aml <- as.data.frame(table(s_aml$Cluster)/nrow(s_aml))

write.csv(p_aml, file = '~/Research/aml/aml_v2/results/for_plots/Pie_pAML.csv', row.names = FALSE, quote = FALSE)
write.csv(s_aml, file = '~/Research/aml/aml_v2/results/for_plots/Pie_sAML.csv', row.names = FALSE, quote = FALSE)

require(ggplot2)
require(dplyr)
require(ggsci)


p_aml <- p_aml %>% 
  arrange(desc(Var1)) %>%
  mutate(prop = Freq  *100) %>%
  mutate(ypos = cumsum(prop)-0.5*prop)
p_aml$Lab <- sprintf('%d%%', round(p_aml$Freq*100))

p_paml <- ggplot(p_aml, aes(x="", y=prop, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  ggtitle('pAML') +
  theme(legend.position="none",
        plot.title = element_text(size=18, hjust = 0.5),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  geom_text(aes(y = ypos, label = Lab), x = 1.2, color = "black", size=4) +
  scale_fill_jco()

s_aml <- s_aml %>% 
  arrange(desc(Var1)) %>%
  mutate(prop = Freq  *100) %>%
  mutate(ypos = cumsum(prop)-0.5*prop)
s_aml$Lab <- sprintf('%d%%', round(s_aml$Freq*100))

p_saml <- ggplot(s_aml, aes(x="", y=prop, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  ggtitle('sAML') +
  theme(plot.title = element_text(size=18, hjust = 0.5),
        legend.position = 'none',
        plot.margin = unit(c(0,0,0,0), "cm")) +
  geom_text(aes(y = ypos, label = Lab), x = 1.2, color = "black", size=4) +
  scale_fill_jco() +
  guides(fill = guide_legend(title='Cluster', title.theme = element_text(face='bold', size=12)))

temp <- ggplot(s_aml, aes(x="", y=prop, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  ggtitle('sAML') +
  theme(plot.title = element_text(size=18, face='bold', hjust = 0.5),
        legend.box.margin = unit(c(0,0,-1,0), "cm"),
        legend.position = 'top') +
  geom_text(aes(y = ypos, label = Lab), color = "black", size=4) +
  scale_fill_jco() +
  guides(fill = guide_legend(title='', title.theme = element_text(face='bold', size=12)))

legend <- cowplot::get_legend(temp)

p_pie <- cowplot::plot_grid(p_paml, p_saml, align='h', ncol = 2, scale = 1)
p <- cowplot::plot_grid(legend, p_pie, align='v', ncol=1, rel_heights = c(0.1,1), axis='t')


cowplot::save_plot(filename = 'plots/MergedPie.pdf',
                   p, base_height = 5, base_width = 7)

## Barplots ##
prob.res.all <- sapply(1:10000, function(i){
  pSam <- sample(rownames(col_annot)[col_annot$Type == 'A-pAML' | col_annot$Type == 'N-pAML'], size = 150, replace = TRUE)
  sSam <- sample(rownames(col_annot)[col_annot$Type == 'A-sAML' | col_annot$Type == 'N-sAML'], size = 150, replace = TRUE)
  rSam <- union(pSam, sSam)
  local_annot <- data.frame('Samples' = rSam,
                            'Clusters' = col_annot$Cluster[match(rSam, table = rownames(col_annot))],
                            'Type' = col_annot$Type[match(rSam, table = rownames(col_annot))])
  local_annot$Type <- ifelse(local_annot$Type %in% c('A-pAML', 'N-pAML'), 'pAML', 'sAML')
  prob.res <- tapply(local_annot$Type, local_annot$Clusters, function(x){
    sum(x == 'pAML') / length(x)
  })
  return(prob.res)
})

plot.data.ft <- data.frame('Cluster' = c(paste0('Cluster-', 1:4), paste0('Cluster-', 1:4)),
                           'Type' = c(rep('pAML', 4), rep('sAML', 4)),
                           'Prob' = c(apply(prob.res.all, 1, mean), 1 - apply(prob.res.all, 1, mean)),
                           'Count' = paste0('n=', c(tapply(col_annot$Type, col_annot$Cluster, function(s){sum(s %in% c('N-pAML', 'A-pAML'))}),
                                                    tapply(col_annot$Type, col_annot$Cluster, function(s){sum(s %in% c('N-sAML', 'A-sAML'))}))))
write.csv(plot.data.ft, file = '~/Research/aml/aml_v2/results/for_plots/BarPlot.csv', row.names = FALSE, quote = FALSE)
write.csv(t(prob.res.all), file = '~/Research/aml/aml_v2/results/for_plots/BarPlotHistogram.csv', row.names = FALSE, quote = FALSE)

p <- ggplot(plot.data.ft, aes(x = Cluster, y = Prob, fill = Type)) +
  geom_col(position = 'dodge') +
  geom_text(aes(label = Count), position = position_dodge(width = 1), vjust = -0.5) +
  xlab('Identified Clusters') + ylab('Frequency') +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 16)) +
  scale_fill_manual(values=setNames(c('red', 'blue'), c('pAML', 'sAML')))

## Histograms ##
clust.1 <- data.frame('Prob' = prob.res.all[1,])
p.1 <- ggplot(clust.1) + 
  geom_histogram(aes(x=Prob), bins = 500, fill = 'black', alpha = 0.6) +
  ylab('') +
  xlab('') +
  xlim(0, 1) +
  geom_vline(xintercept = 0.5, linetype='dashed') +
  theme(panel.background = element_rect(fill=NA),
        axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.1),
        plot.margin = margin(0,0.5,0,0.5,'cm'),
        axis.text.x = element_text(size=4)) +
  scale_y_continuous(breaks=NULL)


clust.2 <- data.frame('Prob' = prob.res.all[2,])
p.2 <- ggplot(clust.2) + 
  geom_histogram(aes(x=Prob), bins = 500, fill = 'black', alpha = 0.6) +
  ylab('') +
  xlab('') +
  xlim(0, 1) +
  geom_vline(xintercept = 0.5, linetype='dashed') +
  theme(panel.background = element_rect(fill=NA),
        axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.1),
        plot.margin = margin(0,0.5,0,0.5,'cm'),
        axis.text.x = element_text(size=4)) +
  scale_y_continuous(breaks=NULL)

clust.3 <- data.frame('Prob' = prob.res.all[3,])
p.3 <- ggplot(clust.3) + 
  geom_histogram(aes(x=Prob), bins = 500, fill = 'black', alpha = 0.6) +
  ylab('') +
  xlab('') +
  xlim(0, 1) +
  geom_vline(xintercept = 0.5, linetype='dashed') +
  theme(panel.background = element_rect(fill=NA),
        axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.1),
        plot.margin = margin(0,0.5,0,0.25,'cm'),
        axis.text.x = element_text(size=4)) +
  scale_y_continuous(breaks=NULL)

clust.4 <- data.frame('Prob' = prob.res.all[4,])
p.4 <- ggplot(clust.4) + 
  geom_histogram(aes(x=Prob), bins = 500, fill = 'black', alpha = 0.6) +
  ylab('') +
  xlab('') +
  xlim(0, 1) +
  geom_vline(xintercept = 0.5, linetype='dashed') +
  theme(panel.background = element_rect(fill=NA),
        axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.1),
        plot.margin = margin(0,2,0,0.5,'cm'),
        axis.text.x = element_text(size=4)) +
  scale_y_continuous(breaks=NULL)

p_hist <- cowplot::plot_grid(p.1, p.2, p.3, p.4, align='h', ncol=4)
save_plot('plots/barplot.pdf', ggdraw(insert_xaxis_grob(p, p_hist, position = 'top')), base_height = 8, base_width = 12)


p_comb <- cowplot::plot_grid(p, p_hist, ncol = 1, align = 'v', axis='t', rel_heights = c(1, 0.2), rel_widths = c(1, 0.7))

## Logistic Regression ##
library(logistf)
setwd('~/Research/aml/aml_v2/')

## READ Data ##
aml_data <- readRDS('data/aml_data_v3.rds')
mut_mat <- aml_data$mut_mat
annot_ft <- aml_data$annot
annot_ft$TypePrimary <- ifelse(annot_ft$Type %in% c('N-pAML', 'A-pAML'), 1, 0)

## Multivariate ##
require(pROC)
auc.vals <- sapply(1:100, function(i){
  p.samples <- sample(rownames(annot_ft)[which(annot_ft$TypePrimary == 1)], size = 370, replace = FALSE)
  s.samples <- sample(rownames(annot_ft)[which(annot_ft$TypePrimary == 0)], size = 370, replace = FALSE)
  c.samples <- union(p.samples, s.samples)
  c.samples <- c.samples[sample(1:length(c.samples), size = length(c.samples), replace = FALSE)]
  train.data <- as.data.frame(mut_mat[match(c.samples, table = rownames(mut_mat)),])
  train.labels <- annot_ft$TypePrimary[match(rownames(train.data), table = rownames(annot_ft))]
  test.data <- as.data.frame(mut_mat[!(rownames(mut_mat) %in% c.samples),])
  test.labels <- annot_ft$TypePrimary[match(rownames(test.data), table = rownames(annot_ft))]
  
  
  ## Logistic Fit
  train.data$Type <- train.labels
  test.data$Type <- test.labels
  
  ## Check
  if(sum(apply(train.data[,-ncol(train.data)], 2, function(x){length(unique(x)) < 2})) > 0){
    return(NA)
  }
  
  logit.fit <- logistf(Type~., data = train.data)
  coef_ft <- logit.fit$coefficients
  
  ## Predict
  pred <- as.vector(exp((as.matrix(test.data[,-ncol(test.data)]) %*% as.matrix(coef_ft)[-1,]) + coef_ft[1]) / (1+exp((as.matrix(test.data[,-ncol(test.data)]) %*% as.matrix(coef_ft)[-1,]) + coef_ft[1])))
  
  ## AUC
  auc.val <- as.numeric(auc(roc(test.data$Type, pred)))
  return(auc.val)
})
pdf('~/aml/results/Logistic_AUC.pdf', width=8, height = 6)
hist(auc.vals, xlab = 'AUC Value', main = 'Multivariate Logistic Regression', breaks = 50)
dev.off()

write.csv(as.data.frame(auc.vals), file = '~/Research/aml/aml_v2/data/multivariate_auc.csv', row.names = FALSE, col.names = FALSE, append = FALSE, quote = FALSE)

clust_res <- read.table('results/MLL_Removed_CombinedTable_LCA_CC_K4.tsv',
                        header = TRUE,
                        sep = '\t',
                        quote = '',
                        stringsAsFactors = FALSE)
clust_assign <- setNames(clust_res$Cluster, clust_res$Id)
clust_assign <- clust_assign[match(rownames(mut_mat), table = names(clust_assign))]

m <- 6

require(randomForest)
imp.res <- sapply(1:50, simplify = FALSE, function(i){
  cat(sprintf('Processing Idx: %d', i))
  rand.idx <- sample(1:nrow(mut_mat), size = 0.8*nrow(mut_mat), replace = FALSE)
  rf <- randomForest(x = mut_mat[rand.idx,], 
                     y = as.factor(clust_assign[rand.idx]),
                     ntree = 1500,
                     mtry = m,
                     replace = TRUE, 
                     importance = FALSE)
  local.pred <- predict(rf, newdata = mut_mat[-rand.idx,])
  local.pred <- table(local.pred, clust_assign[-rand.idx])
  return(sum(diag(local.pred))/sum(local.pred))
})
write.csv(as.data.frame(imp.res), file = '~/Research/aml/aml_v2/data/rf_auc.csv', row.names = FALSE, col.names = FALSE, append = FALSE, quote = FALSE)

## For Plotting ##
library(party)
sample_idx <- colnames(mut_mat) %in% c('ASXL1', 'IDH2', 'NPM1', 'RUNX1', 'SRSF2', 'TP53', 'NPM1', 'DNMT3A', 'del5q', 'del17p', 'FLT3(TKD or ITD)')
local_mut_mat <- as.data.frame(mut_mat[,sample_idx])
local_mut_mat$Type <- as.factor(clust_assign)

rf_res <- ctree(Type~., data=local_mut_mat)
library(Cairo)
CairoPDF('~/Research/aml/aml_v2/plots/SampleTree.pdf', width=32, height = 8)
plot(rf_res, type='extended')
dev.off()

aml.data <- readRDS('~/Research/aml/aml_v2/data/aml_data_v4.rds')
ids <- aml.data$sample.id
mut_mat <- data.matrix(aml.data[,-1])
rownames(mut_mat) <- ids

clust_res <- read.csv('~/Research/aml/aml_v2/results/aml_cluster_v4/ClusterResults_V4.csv')
for(c in unique(clust_res$Cluster)){
  local_mut_mat <- mut_mat[rownames(mut_mat) %in% clust_res$Id[clust_res$Cluster == c],]
  require(circlize)
  circos.initialize(factors=c('Chr1', 'Chr2', 'Chr3'), xlim=c(0,1))
  circos.track(ylim=c(0,1), panel.fun=function(x,y){
    xlim=c(0,1)
    ylim=c(0,1)
    circos.rect(xlim[1], 0, xlim[2], 1, col=rand_color(1))
  }, track.height=0.15, bg.border=NA)
}