library(logistf)

prep_mat <- function(x=NULL, filter=TRUE){
  x <- as.data.frame(ifelse(x>=1, 'Mut', 'NonMut'))
  if(filter)
    x <- x[,apply(x, 2, function(k){min(c(sum(k=='Mut', na.rm = TRUE),sum(k=='NonMut', na.rm = TRUE))) > 0.00*nrow(x)})]
  for(i in 1:ncol(x)){
    x[,i] <- factor(x[,i], 
                    levels = c('NonMut', 'Mut'))
  }
  return(x)
}

library(readxl)
aml.data <- read_excel('~/Research/aml/aml_v2/data/ULR-MLR.xlsx',
                       sheet = 1,
                       na = c('', '-', 'NA', 'NaN'))
aml.data.annot <- data.frame('id' = paste0('sample.', aml.data$ID),
                             'OS' = as.numeric(trimws(aml.data$OS)),
                             'Event' = as.integer(trimws(aml.data$`D=1/A=0`)),
                             'Type' = trimws(aml.data$Type))
aml.data.ft <- data.frame('inv3' = as.integer(trimws(aml.data$`inv(3)/ t(3;3)`)),
                          'tr6' = as.integer(trimws(aml.data$`t(6;9)(p23;q34)`)),
                          'del5q' = as.integer(trimws(aml.data$`–5/del(5q)`)),
                          'del6q' = as.integer(trimws(aml.data$`del6q/-6`)),
                          'del7q' = as.integer(trimws(aml.data$`–7/del(7q)`)),
                          'del9q' = as.integer(trimws(aml.data$`–9/del(9q)`)),
                          'del12p' = as.integer(trimws(aml.data$`del(12p)`)),
                          'del13q' = as.integer(trimws(aml.data$`del(13)(q)`)),
                          'del16q' = as.integer(trimws(aml.data$`del(16q)`)),
                          'del17p' = as.integer(trimws(aml.data$`–17/del(17p)`)),
                          'del20q' = as.integer(trimws(aml.data$`del(20q)`)),
                          'trisomy8' = as.integer(trimws(aml.data$`trisomoy 8`)),
                          'Complex' = as.integer(trimws(aml.data$Complex)),
                          'X' = as.integer(trimws(aml.data$`–X`)),
                          'Y' = as.integer(trimws(aml.data$`–Y`)),
                          'ASXL1' = as.integer(trimws(aml.data$ASXL1)),
                          'BCOR.L1' = as.integer(trimws(aml.data$`BCOR/L1`)),
                          'CBL' = as.integer(trimws(aml.data$CBL)),
                          'CEBPA.Mono' = as.integer(trimws(aml.data$`M-CEBPA`)),
                          'CEBPA.Bi' = as.integer(trimws(aml.data$`Bi-CEBPA`)),
                          'DNMT3A' = as.integer(trimws(aml.data$DNMT3A)),
                          'ETV6' = as.integer(trimws(aml.data$ETV6)),
                          'EZH2' = as.integer(trimws(aml.data$EZH2)),
                          'FLT3.ITD' = as.integer(trimws(aml.data$`FLT3-ITD`)),
                          'FLT3.TKD' = as.integer(trimws(aml.data$`FLT3-TKD`)),
                          'GATA2' = as.integer(trimws(aml.data$GATA2)),
                          'IDH2.140' = as.integer(trimws(aml.data$`IDH2-140`)),
                          'IDH2.172' = as.integer(trimws(aml.data$`IDH2-172`)),
                          'IDH1' = as.integer(trimws(aml.data$IDH1)),
                          'KIT' = as.integer(trimws(aml.data$KIT)),
                          'KRAS' = as.integer(trimws(aml.data$KRAS)),
                          'NPM1' = as.integer(trimws(aml.data$NPM1)),
                          'NRAS' = as.integer(trimws(aml.data$NRAS)),
                          'RUNX1' = as.integer(trimws(aml.data$RUNX1)),
                          'SF3B1' = as.integer(trimws(aml.data$SF3B1)),
                          'SRSF2' = as.integer(trimws(aml.data$SRSF2)),
                          'TET2' = as.integer(trimws(aml.data$TET2)),
                          'TP53' = as.integer(trimws(aml.data$TP53)),
                          'U2AF1' = as.integer(trimws(aml.data$U2AF1)),
                          'WT1' = as.integer(trimws(aml.data$WT1)),
                          'ZRSR2' = as.integer(trimws(aml.data$ZRSR2)))

aml_annot <- aml.data.annot
aml_mut <- prep_mat(aml.data.ft)
rownames(aml_mut) <- aml_annot$id

aml_mut <- aml_mut[,colSums(is.na(aml_mut)) < nrow(aml_mut)*0.6]
aml_mut <- aml_mut[rowSums(is.na(aml_mut)) < 0.6*ncol(aml_mut),]
aml_mut <- aml_mut[, apply(aml_mut, 2, function(x){min(sum(x=='Mut', na.rm = TRUE), sum(x=='NonMut', na.rm = TRUE))}) > 50]
aml_mut <- aml_mut[rowSums(is.na(aml_mut)) < 0.6*ncol(aml_mut),]
aml_annot <- aml_annot[match(rownames(aml_mut), table=aml_annot$id),]

## Impute ##
require(missForest)
aml_mut_imputed <- missForest(aml_mut, variablewise = TRUE, ntree = 2500, maxiter = 10)
aml_mut_imputed <- aml_mut_imputed$ximp
aml_annot$TypePrimary <- ifelse(aml_annot$Type %in% c('N-pAML', 'A-pAML'), 1, 0)

## Multivariate ##
require(pROC)
auc.vals <- sapply(1:100, function(i){
  p.samples <- sample(aml_annot$id[which(aml_annot$TypePrimary == 1)], size = 250, replace = FALSE)
  s.samples <- sample(aml_annot$id[which(aml_annot$TypePrimary == 0)], size = 250, replace = FALSE)
  c.samples <- union(p.samples, s.samples)
  c.samples <- c.samples[sample(1:length(c.samples), size = length(c.samples), replace = FALSE)]
  train.data <- as.data.frame(aml_mut_imputed[match(c.samples, table = rownames(aml_mut_imputed)),])
  train.labels <- aml_annot$TypePrimary[match(rownames(train.data), table = aml_annot$id)]
  test.data <- as.data.frame(aml_mut_imputed[!(rownames(aml_mut_imputed) %in% c.samples),])
  test.labels <- aml_annot$TypePrimary[match(rownames(test.data), table = aml_annot$id)]
  
  ## Logistic Fit
  train.data$Type <- train.labels
  test.data$Type <- test.labels
  
  # temp <- as.matrix(train.data)
  # temp[temp == 'Mut'] <- 1
  # temp[temp == 'NonMut'] <- 0
  # temp <- as.matrix(test.data)
  # temp[temp == 'Mut'] <- 1
  # temp[temp == 'NonMut'] <- 0
  # test.data <- as.data.frame(temp)
  
  logit.fit <- glm(Type~., data = train.data, family = binomial(link='logit'))
  
  ## Predict
  pred <- predict(logit.fit, newdata = test.data, type='response')
  
  ## AUC
  auc.val <- as.numeric(auc(roc(test.data$Type, pred)))
  return(auc.val)
})
pdf('~/Research/aml/aml_v2/results/Logistic_AUC_V4.pdf', width=8, height = 6)
hist(auc.vals, xlab = 'AUC Value', main = 'Multivariate Logistic Regression', breaks = 50)
dev.off()


## Plot Model effects using all data ##
p.samples <- aml_annot$id[which(aml_annot$TypePrimary == 1)]
s.samples <- aml_annot$id[which(aml_annot$TypePrimary == 0)]
c.samples <- union(p.samples, s.samples)
c.samples <- c.samples[sample(1:length(c.samples), size = length(c.samples), replace = FALSE)]
train.data <- as.data.frame(aml_mut_imputed[match(c.samples, table = rownames(aml_mut_imputed)),])
train.labels <- aml_annot$TypePrimary[match(rownames(train.data), table = aml_annot$id)]

## Logistic Fit
train.data$Type <- train.labels
logit.fit <- glm(Type~., 
                 data = train.data, 
                 family = binomial(link='logit'), 
                 weights = ifelse(train.data$Type == 1, 0.15, 0.85))
conf.res <- confint(logit.fit)
temp <- summary(logit.fit)
combined <- data.frame('Estimate' = temp$coefficients[,1],
                       'ciLow' = conf.res[,1],
                       'ciHigh' = conf.res[,2],
                       'pVal' = temp$coefficients[,4])
combined <- combined[-1,]

write.table(combined, '~/Research/aml/aml_v2/results/Logistic_Multivariate_V4.tsv', col.names = TRUE, row.names = TRUE, sep = '\t', append = FALSE, quote = FALSE)

plot.data <- data.frame('Estimate' = combined[,1],
                        'ciLow' = combined[,2],
                        'ciHigh' = combined[,3])

plot.data$Variable <- gsub(pattern = '`', replacement = '', rownames(combined))
plot.data$Variable <- factor(plot.data$Variable,
                             levels = names(sort(setNames(plot.data$Estimate, plot.data$Variable), decreasing = FALSE)))

require(rcompanion)
require(ggplot2)
p <- ggplot(plot.data, aes(x = Estimate, y = Variable)) +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = 0.5, height = 1, color = "gray50") +
  geom_point(size = 2, color = "orange") +
  xlim(-6,6) +
  ylab("") +
  xlab("Odds ratio (log scale)") +
  #annotate(geom = "text", y =1.1, x = 3.5, label ="Model p < 0.001\nPseudo R^2 = 0.24", size = 3.5, hjust = 0) + 
  ggtitle("Primary vs Secondary AML") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'))

ggsave(p, filename = 'results/Logistic-Effects-Multivariate.pdf', width = 8, height = 6, device = cairo_pdf)

## Univariate ##
uni_res <- sapply(colnames(aml_mut_imputed), simplify = FALSE, function(x){
  local.train.data <- data.frame('variable' = aml_mut_imputed[,colnames(aml_mut_imputed) %in% x],
                                 'type' = aml_annot$Type,
                                 'type.primary' = ifelse(aml_annot$Type %in% c('N-pAML', 'A-pAML'), 'pAML', 'sAML'))
  local.train.data <- subset(local.train.data, local.train.data$type %in% c('A-pAML', 'A-sAML'))
  local.train.data$type <- factor(local.train.data$type, levels=c('A-sAML', 'A-pAML'))
  #local.train.data$type.primary <- factor(local.train.data$type.primary, levels=c('sAML', 'pAML'))
  tryCatch(expr={
    logit.fit <- glm(type~variable, 
                     data = local.train.data, 
                     family = binomial(link='logit'), 
                     weights = ifelse(local.train.data$type == 'N-pAML', 0.1, 0.9))
    conf.res <- confint(logit.fit)
    temp <- summary(logit.fit)
    combined <- data.frame('Estimate' = temp$coefficients[,1],
                           'ciLow' = conf.res[,1],
                           'ciHigh' = conf.res[,2],
                           'pVal' = temp$coefficients[,4])
    return(combined[-1,])
  }, error=function(err){
    return(NULL)
  })
})
uni_res <- do.call('rbind', uni_res)
write.table(uni_res, '~/Research/aml/aml_v2/results/Logistic_Univariate_APvsAS.tsv', col.names = TRUE, row.names = TRUE, sep = '\t', append = FALSE, quote = FALSE)
