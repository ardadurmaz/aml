library(readxl)
library(randomForest)
require(rms)
require(Hmisc)

## Survival Validation ##
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
                       'Age' = annot_ft$Age[match(clust_res$Id, table=annot_ft$Id)],
                       'Data' = rep('train', nrow(clust_res)))
id_tag <- train_ft$Cluster == 'Cluster-1'
train_ft$Cluster[train_ft$Cluster == 'Cluster-2'] <- 'Cluster-1'
train_ft$Cluster[id_tag] <- 'Cluster-2'

test_ft <- data.frame('id' = names(pred_res),
                      'Cluster' = pred_res,
                      'OS' = annot_ft$OS[match(names(pred_res), table = annot_ft$Id)],
                      'Event' = annot_ft$Event[match(names(pred_res), table = annot_ft$Id)],
                      'Age' = annot_ft$Age[match(names(pred_res), table=annot_ft$Id)],
                      'Data' = rep('test', length(pred_res)))

id_tag <- test_ft$Cluster == 'Cluster-1'
test_ft$Cluster[test_ft$Cluster == 'Cluster-2'] <- 'Cluster-1'
test_ft$Cluster[id_tag] <- 'Cluster-2'

combined_ft <- rbind(test_ft, train_ft)
rownames(combined_ft) <- 1:nrow(combined_ft)
combined_ft$Sex <- annot_ft$Sex[match(combined_ft$id, table=annot_ft$Id)]
combined_ft$WBC <- annot_ft$WBC[match(combined_ft$id, table=annot_ft$Id)]
combined_ft$Leukopenia <- annot_ft$Leukopenia[match(combined_ft$id, table=annot_ft$Id)]
combined_ft$HB <- annot_ft$HB[match(combined_ft$id, table=annot_ft$Id)]
combined_ft$Anemia <- annot_ft$Anemia[match(combined_ft$id, table=annot_ft$Id)]
combined_ft$Platelet <- annot_ft$Platelet[match(combined_ft$id, table=annot_ft$Id)]
combined_ft$Thrombocytopenia <- annot_ft$Thrombocytopenia[match(combined_ft$id, table=annot_ft$Id)]
combined_ft$BM.Blasts <- annot_ft$BM_Blasts[match(combined_ft$id, table=annot_ft$Id)]

prep_mat <- function(x=NULL, filter=TRUE){
  x <- as.data.frame(ifelse(x>=1, 'Mut', 'NonMut'))
  if(filter)
    x <- x[,apply(x, 2, function(k){min(c(sum(k=='Mut', na.rm = TRUE),sum(k=='NonMut', na.rm = TRUE))) > 15})]
  for(i in 1:ncol(x)){
    x[,i] <- factor(x[,i], 
                    levels = c('NonMut', 'Mut'))
  }
  return(x)
}
require(survAUC)

## Cluster 1 ##
local_ft <- subset(combined_ft, combined_ft$Cluster == 'Cluster-1')
local_ft <- na.omit(cbind(as.data.frame(local_ft[,c(3,4,5,6,7,8,10,12)]), as.data.frame(prep_mat(local_ft[,c(9,11,13)]))))
local_ft$Age.Cat <- ifelse(local_ft$Age > 60, 'High', 'Low')
local_ft_train <- subset(local_ft, local_ft$Data == 'train')
local_ft_test <- subset(local_ft, local_ft$Data == 'test')

library(MASS)
coxph_model <- coxph(Surv(OS, Event)~strata(Age.Cat)+Sex+WBC+HB+Platelet+Leukopenia+Anemia+Thrombocytopenia, data=local_ft_train)
stepAIC(coxph_model, direction = 'both', k = 2)
coxph_model_step <- coxph(Surv(OS, Event)~strata(Age.Cat)+Sex+WBC+Leukopenia, data=local_ft_train,x=TRUE,y=TRUE)
cox.zph(coxph_model_step)

## Skip AUC for now ##
require(survAUC)
pred_res_train <- predict(coxph_model_step, newdata = local_ft_train, type='survival')
pred_res_test <- predict(coxph_model_step, newdata = local_ft_test, type='survival')
Surv.rsp <- Surv(local_ft_train$OS, local_ft_train$Event)
Surv.rsp.new <- Surv(local_ft_test$OS, local_ft_test$Event)
times <- seq(10, 50, 100)
AUC_cd <- AUC.cd(Surv.rsp, Surv.rsp.new, lp=pred_res_train, lpnew=pred_res_test, times=30)


feat_mapping <- setNames(colnames(mut_mat), paste0('feat.', 1:ncol(mut_mat)))
surv_pred_res <- sapply(paste0('Cluster-', 1:4), simplify = FALSE, function(clust){
  train_ft_local <- subset(train_ft, train_ft$Cluster == clust)
  test_ft_local <- subset(test_ft, test_ft$Cluster == clust)
  
  train_surv_ft <- cbind(data.frame('OS' = train_ft_local$OS,
                                    'Event' = train_ft_local$Event,
                                    'Age' = train_ft_local$Age),
                         prep_mat(as.data.frame(ifelse(mut_mat[match(train_ft_local$id, table=rownames(mut_mat)),]==1,'Mut','NonMut')), feat_map = feat_mapping))
  test_surv_ft <- cbind(data.frame('OS' = test_ft_local$OS,
                                   'Event' = test_ft_local$Event,
                                   'Age' = test_ft_local$Age),
                        prep_mat(as.data.frame(ifelse(mut_mat[match(test_ft_local$id, table=rownames(mut_mat)),]==1,'Mut','NonMut')), feat_map = feat_mapping, filter = FALSE))
  ## Cox-PH ##
  cox_fit <- rms::cph(Surv(OS, Event)~., data=train_surv_ft, x=TRUE, y=TRUE, surv=TRUE)
  sink(file = sprintf('~/Research/aml/aml_v2/results/survival_%s.tsv', clust))
  print(cox_fit)
  sink()
  surv_pred <- survest(cox_fit, newdata=test_surv_ft, times=175)
  surv_est <- surv_pred$surv
  return(rcorr.cens(x=surv_est, S=Surv(test_surv_ft$OS, test_surv_ft$Event)))
  
})
surv_pred_res <- do.call('cbind', surv_pred_res)
write.csv(surv_pred_res, file='~/Research/aml/aml_v2/results/survival_predictions_InternalTest.csv', quote = FALSE, row.names = TRUE, col.names = TRUE, append = FALSE)


## External Data ##
setwd('~/Research/aml/aml_v2/')
valid_data <- read_xlsx('data/MDACC.xlsx', sheet=1)
rf_model <- readRDS('data/aml_rf_model_v4.rds')

valid_mut <- data.frame('inv3' = as.integer(trimws(valid_data$`inv(3)/t(3;3)`)),
                        'tr6' = as.integer(trimws(valid_data$`t(6;9)`)),
                        'del5q' = valid_data$`–5/del(5q)`,
                        'del6q' = valid_data$`del6q/-6`,
                        'del7q' = valid_data$`–7/del(7q)`,
                        'del9q' = valid_data$`–9/del(9q)`,
                        'del12p' = valid_data$`del(12p)`,
                        'del13q' = valid_data$`del(13)(q)`,
                        'del16q' = valid_data$`del(16q)`,
                        'del17p' = valid_data$`–17/del(17p)`,
                        'del20q' = valid_data$`del(20q)`,
                        'trisomy8' = valid_data$`trisomoy 8`,
                        'X' = valid_data$`–X`,
                        'Y' = valid_data$`–Y`,
                        'ASXL1' = valid_data$ASXL1,
                        'BCOR.L1' = valid_data$`BCOR/L1`,
                        'CBL' = valid_data$CBL,
                        'CEBPA.Mono' = valid_data$`CEBPA-Mono`,
                        'CEBPA.Bi' = valid_data$`CEBPA-Bi`,
                        'DNMT3A' = valid_data$DNMT3A,
                        'ETV6' = valid_data$ETV6,
                        'EZH2' = valid_data$EZH2,
                        'FLT3.TKD' = valid_data$`FLT-3 TKD`,
                        'FLT3.ITD' = valid_data$FLT3_ITD,
                        'GATA2' = valid_data$GATA2,
                        'IDH1' = valid_data$IDH1,
                        'IDH2.140' = valid_data$`IDH2-140`,
                        'IDH2.172' = valid_data$`IDH2-172`,
                        'KIT' = valid_data$KIT,
                        'KRAS' = valid_data$KRAS,
                        'NPM1' = valid_data$NPM1,
                        'NRAS' = valid_data$NRAS,
                        'RUNX1' = valid_data$RUNX1,
                        'SF3B1' = valid_data$SF3B1,
                        'SRSF2' = valid_data$SRSF2,
                        'TET2' = valid_data$TET2,
                        'TP53' = valid_data$TP53,
                        'U2AF1' = valid_data$U2AF1,
                        'WT1' = valid_data$WT1,
                        'ZRSR2' = valid_data$ZRSR2, 
                        check.names = FALSE)
rownames(valid_mut) <- paste0('sample.',valid_data$Accession)
valid_mut[valid_mut>1] <- 1

valid_mut <- valid_mut[rowSums(is.na(valid_mut)) < 5,]
valid_mut <- na.roughfix(valid_mut)
pred_res <- predict(rf_model, newdata = valid_mut)

surv_ft <- data.frame('id' = rownames(valid_mut),
                      'OS' = valid_data$`OS (Months)`[match(rownames(valid_mut), table=paste0('sample.', valid_data$Accession))],
                      'Event' = valid_data$Died[match(rownames(valid_mut), table=paste0('sample.', valid_data$Accession))],
                      'Cluster' = pred_res)
## Kaplan-Meier Plot ##
write.csv(surv_ft, file = '~/Research/aml/aml_v2/results/MDACC_Survival.csv', append = FALSE, row.names = FALSE, col.names = TRUE, quote=FALSE)

fit <- survfit(Surv(OS, Event) ~ Cluster, data = surv_ft)
logrank_res <- survdiff(Surv(OS, Event) ~ Cluster, data = surv_ft)
p_val_res <- pchisq(31.4, df=3, lower.tail = FALSE)
cluster_colors <- setNames(c('#ff7f00',
                             '#c0c000',
                             '#04bf78',
                             '#b956d7'), nm = paste0('Cluster-', 1:4))


## CoxPH ##
surv_ft_train <- na.omit(surv_ft_train)
surv_ft_test <- na.omit(surv_ft_test)

cox_fit <- coxph(Surv(OS, Event) ~ Cluster, data=surv_ft_train, x=TRUE, y=TRUE, method='efron')
lp <- predict(cox_fit)
lpnew <- predict(cox_fit, newdata = surv_ft_test)
Surv.rsp <- Surv(surv_ft_train$OS, surv_ft_train$Event)
Surv.rsp.new <- Surv(surv_ft_test$OS, surv_ft_test$Event)
times <- c(24)

AUC_hc <- AUC.hc(Surv.rsp, Surv.rsp.new, lpnew, times)
AUC_sh <- AUC.sh(Surv.rsp, Surv.rsp.new, lp, lpnew, times)

require(survminer)
require(survival)
hr_res <- sapply(1:250, simplify = FALSE, function(i){
  local_surv_ft <- surv_ft[sample(1:nrow(surv_ft), size=nrow(surv_ft)*0.8, replace = TRUE),]
  cox_fit <- coxph(Surv(OS, Event)~Cluster, data = local_surv_ft)
  return(exp(cox_fit$coefficients))
})
hr_res <- do.call('rbind', hr_res)
hr_res_validation <- hr_res
hr_ft <- rbind(colMeans(hr_res),
               apply(hr_res, 2, sd))
rownames(hr_ft) <- c('HR', 'SD')
write.table(hr_ft, file = '~/Research/aml/aml_v2/results/SurvivalValidation.tsv', sep = '\t', append = FALSE, quote = FALSE, row.names = TRUE, col.names = TRUE)

## Training Data ##
surv_ft <- read.table('results/SurvivalData.tsv', sep = '\t', header = TRUE)
surv_ft <- na.omit(surv_ft)
surv_ft$Cluster <- factor(surv_ft$Cluster, levels = paste0('Cluster-', c(4,3,2,1)))

hr_res <- sapply(1:250, simplify = FALSE, function(i){
  local_surv_ft <- surv_ft[sample(1:nrow(surv_ft), size=nrow(surv_ft)*0.8, replace = TRUE),]
  if(sum(tapply(local_surv_ft$Event, local_surv_ft$Cluster, function(x){sum(x, na.rm = TRUE)}) > 0) == 4){
    cox_fit <- coxph(Surv(OS, Event)~Cluster, data = local_surv_ft)
    return(exp(cox_fit$coefficients))
  }else{
    return(NULL)
  }
})
hr_res <- do.call('rbind', hr_res)
hr_res_training <- hr_res
hr_ft <- rbind(colMeans(hr_res),
               apply(hr_res, 2, sd))
rownames(hr_ft) <- c('HR', 'SD')
write.table(hr_ft, file = '~/Research/aml/aml_v2/results/SurvivalValidationTrainingData.tsv', sep = '\t', append = FALSE, quote = FALSE, row.names = TRUE, col.names = TRUE)

## Boxplot ##
train_res <- reshape2::melt(hr_res_training)
colnames(train_res) <- c('Run', 'Comp', 'HR')
train_res$type <- 'Training'

valid_res <- reshape2::melt(hr_res_validation)
colnames(valid_res) <- c('Run', 'Comp', 'HR')
valid_res$type <- 'Validation'

combined_res <- rbind(train_res, valid_res)

p <- ggplot(combined_res, aes(x=Comp, y=HR, fill=type)) +
  geom_boxplot(notch=TRUE) +
  ylab('Hazard Ratio') +
  xlab('') +
  labs(fill='') +
  theme_classic() +
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=12, face='bold')) +
  scale_x_discrete(labels=c('ClusterCluster-3' = 'Cluster-3', 
                            'ClusterCluster-2' = 'Cluster-2',
                            'ClusterCluster-1' = 'Cluster1'))
ggsave(p, filename = '~/Research/aml/aml_v2/plots/SurvivalValidation.pdf')