library(shiny)
library(shinycssloaders)
library(ggplot2)
library(readxl)
library(randomForest)
library(FactoMineR)
library(cluster)
library(ggplot2)
library(shinythemes)
library(survminer)
library(survival)
library(coxed)
library(ggsci)

# Define UI for data upload app ----
ui <- fluidPage(
  
  theme = shinytheme('lumen'),
  # App title ----
  titlePanel("AML Genomic Subclassification"),
  tags$hr(),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      selectInput("cyto", "Select Cytogenetic Abnormalities",
                  c('inv(3)/t(3;3)','t(6;9)','-5/del(5q)', '-6/del(6q)', '-7/(del7q)', '-9/del(9q)', 'del(12p)', 'del(13q)', 'del(16q)', '-17/del(17p)', 'del(20q)', '+8', '-X', '-Y'), 
                  multiple = TRUE),
      selectInput('mut', 'Select Gene Mutations',
                  c('ASXL1','BCOR/L1','CBL','CEBPA (Monoallelic)','CEBPA (Biallelic)','DNMT3A','ETV6','EZH2','FLT3 (ITD)', 'FLT3 (TKD)', 'GATA2','IDH1','IDH2 (R140)','IDH2 (R172)', 'KIT','KRAS','NPM1','NRAS','RUNX1','SF3B1','SRSF2','TET2','TP53','U2AF1','WT1','ZRSR2'),
                  multiple = TRUE),
      actionButton(inputId = 'predict_rf', label = "Predict")
    ),
    mainPanel(
      width = 10,
      plotOutput('rf_results')
    )
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  plot_data <- eventReactive(input$predict_rf, {
    rf_model <- readRDS('app_data/aml_rf_model_v4.rds')
    surv_ft <- read.csv('app_data/SurvivalData_V4.tsv')
    
    ## Prediction ##
    x <- union(input$mut, input$cyto)
    if(is.null(x)){
      x <- numeric(nrow(rf_model$importance))
    }else{
      cyto_map <- setNames(c('inv(3)/t(3;3)','t(6;9)','-5/del(5q)', '-6/del(6q)', '-7/(del7q)', '-9/del(9q)', 'del(12p)', 'del(13q)', 'del(16q)', '-17/del(17p)', 'del(20q)', '+8', '-X', '-Y', 'FLT3 (ITD)', 'FLT3 (TKD)','CEBPA (Monoallelic)','CEBPA (Biallelic)','IDH2 (R140)','IDH2 (R172)', 'BCOR/L1'),
                           c('inv3','tr6','del5q', 'del6q', 'del7q', 'del9q', 'del12p', 'del13q', "del16q", "del17p", "del20q","trisomy8","X","Y", 'FLT3.ITD', 'FLT3.TKD', 'CEBPA.Mono', 'CEBPA.Bi', 'IDH2.140', 'IDH2.172', 'BCOR.L1'))
      if(any(x %in% cyto_map)){
        cyto.idx <- which(x %in% cyto_map)
        x[cyto.idx] <- names(cyto_map)[match(x[cyto.idx], table=cyto_map)]
      }
      x <- ifelse(rownames(rf_model$importance) %in% x, 1, 0)
    }
    cluster_colors <- setNames(c('#c0c000',
                                 '#04bf78',
                                 '#ff7f00',
                                 '#b956d7'), nm = paste0('Cluster-', 1:4))
    
    local_predict <- predict(rf_model, newdata = matrix(x, nrow = 1), type = 'prob')
    plot_ft <- reshape2::melt(local_predict)
    colnames(plot_ft) <- c('Input', 'Cluster', 'Prob')
    plot_ft$Input <- paste0('Input-', plot_ft$Input)
    p_bar <- ggplot(plot_ft, aes(x = Input, y = Prob, fill = Cluster)) +
      geom_col(stat='identity', position='dodge') +
      xlab('') +
      ylab('Probability') +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      scale_fill_manual(values = cluster_colors,
                        labels = c('Genomic Cluster 1', 'Genomic Cluster 2', 'Genomic Cluster 3', 'Genomic Cluster 4')) +
      guides(fill=guide_legend(title=''))
    
    ## Survival ##
    pred_cluster <- plot_ft$Cluster[which.max(plot_ft$Prob)]
    fit <- survfit(Surv(OS, Event) ~ Cluster, data = surv_ft)
    
    ## Mod Colors ##
    temp <- col2rgb(as.vector(cluster_colors))
    temp <- sapply(1:ncol(temp), function(c){rgb(temp[1,c]/255, temp[2,c]/255, temp[3,c]/255, alpha = ifelse(paste0('Cluster-', c) == pred_cluster, 1, 0.05))})
    
    p_surv <- ggsurvplot(
      fit, 
      data = surv_ft, 
      size = 0.6,
      conf.int = TRUE,
      pval = FALSE,
      risk.table = FALSE,    
      palette = temp,
      conf.int.style = 'step',
      legend = 'none',
      xlab='Time (Months)',
      ylab='Survival Probability',
      ggtheme = theme_classic(base_size = 16))
    
    temp <- summary(fit, times=c(12,24,36,60))
    local.lab <- data.frame('Group' = temp$strata,
                            'Times' = temp$time,
                            'Survival' = temp$surv,
                            'loCI' = temp$lower,
                            'hiCI' = temp$upper)
    med <- data.frame('Group' = rownames(temp$table),
                      'Median' = temp$table[,7],
                      'loCI' = temp$table[,8],
                      'hiCI' = temp$table[,9])
    local.lab <- subset(local.lab, local.lab$Group == sprintf('Cluster=%s', pred_cluster))
    med <- subset(med, med$Group == sprintf('Cluster=%s', pred_cluster))
    label_ft <- sprintf("Median Survival: %.2f [%.2f-%.2f]\nSurvival Probability (1-Year): %.2f [%.2f-%.2f]\nSurvival Probability (2-Year): %.2f [%.2f-%.2f]\nSurvival Probability (3-Year): %.2f [%.2f-%.2f]\nSurvival Probability (5-Year): %.2f [%.2f-%.2f]", med$Median, med$loCI, med$hiCI,local.lab$Survival[1], local.lab$loCI[1], local.lab$hiCI[1],local.lab$Survival[2], local.lab$loCI[2], local.lab$hiCI[2],local.lab$Survival[3], local.lab$loCI[3], local.lab$hiCI[3],local.lab$Survival[4], local.lab$loCI[4], local.lab$hiCI[4])
    p_surv$plot <- p_surv$plot +
      ggplot2::annotate("text", x = 125, y = 0.8, label = label_ft, size = 4)
    p_comb <- gridExtra::arrangeGrob(p_surv$plot, p_bar, nrow = 1)
    p_comb
  })
  output$rf_results <- renderPlot({
    local_p <- plot_data()
    grid::grid.draw(local_p)
  })
}

# Create Shiny ##
shinyApp(ui, server)