library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(GEOquery)
library(R.utils)
library(reshape2)
library(limma)
library(dplyr)
library(tidyr)
library(caret)
library(MLmetrics)

options(shiny.maxRequestSize = 200*1024^2)

# Load and preprocess data
set.seed(123)
gse_15852 <- getGEO("GSE15852")[[1]]
exprs_data <- log2(exprs(gse_15852))
pdata <- pData(gse_15852)
grade_vec <- pdata$`grade:ch1`

# PCA
expr_for_pca <- t(exprs_data)
pca_res <- prcomp(expr_for_pca, center = TRUE, scale. = TRUE)
grade_fac <- factor(grade_vec, levels = c("normal", "grade 1", "grade 2", "grade 3"))
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Grade = grade_fac)

# Classification
train_idx <- createDataPartition(grade_vec, p = 0.8, list = FALSE)
train_samples <- colnames(exprs_data)[train_idx]
test_samples  <- colnames(exprs_data)[-train_idx]
expr_train <- exprs_data[, train_samples]
expr_test  <- exprs_data[, test_samples]
labels_train <- factor(grade_vec[train_idx], levels = c("normal", "grade 1", "grade 2", "grade 3"))
labels_test  <- factor(grade_vec[-train_idx], levels = c("normal", "grade 1", "grade 2", "grade 3"))
levels(labels_train) <- make.names(levels(labels_train))
levels(labels_test)  <- make.names(levels(labels_test))

# Trend analysis
grade_map <- c("normal"=0, "grade.1"=1, "grade.2"=2, "grade.3"=3)
grade_num <- grade_map[as.character(labels_train)]
design_trend <- model.matrix(~ grade_num)
fit_trend <- lmFit(expr_train, design_trend)
fit_trend <- eBayes(fit_trend)
res_trend <- topTable(fit_trend, coef = "grade_num", number = Inf, adjust.method = "BH")
sig_genes <- rownames(subset(res_trend, adj.P.Val < 0.05 & abs(logFC) > 1))

# Model training and evaluation
x_train <- t(expr_train[sig_genes, , drop = FALSE])
x_test  <- t(expr_test [sig_genes, , drop = FALSE])
y_train <- labels_train
y_test  <- labels_test

trctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = multiClassSummary, savePredictions = "final")
svm_mod <- train(x = x_train, y = y_train, method = "svmLinear", trControl = trctrl)
pred_test <- predict(svm_mod, newdata = x_test)
predicted_probs <- predict(svm_mod, newdata = x_test, type = "prob")
cm <- confusionMatrix(pred_test, y_test)
accuracy_val <- cm$overall["Accuracy"]
kappa_val <- cm$overall["Kappa"]

results_df <- data.frame(SampleID = rownames(x_test), TrueGrade = y_test, PredictedGrade = pred_test, round(predicted_probs * 100, 1))

ui <- fluidPage(
  theme = shinytheme("cerulean"),
  navbarPage(
    id = "main_tabs",
    "Breast Cancer Tool",
    
    tabPanel("About",
             mainPanel(
               p("Breast cancer is the most common cancer among females, and one in seven women are diagnosed with it in their lifetime. This is why we created (Insert app name). This is a tool that medical researchers can use to enhance knowledge that facilitates early breast cancer detection, and understanding the interaction between genes and breast cancer diagnoses. But also to detect it when it matters, and know what stage the cancer is at for further action.This Shiny web application uses machine learning models to help predict and visualise breast cancer diagnoses. The application, which was developed using R and is driven by interactive data visualisation, gives users an easy-to-use way to examine patient data and prediction outcomes.

Sources:

National Breast Cancer Foundation (2021). Breast Cancer Statistics in Australia. [online] National Breast Cancer Foundation. Available at: https://nbcf.org.au/about-breast-cancer/breast-cancer-stats/.
"))) ,
    
    tabPanel("How to use", "Insert how to use it here"),
    
    tabPanel("Data",
             sidebarLayout(position = "right",
                           sidebarPanel(
                             h4("Upload your CSV file"),
                             fileInput("file1", "Choose CSV File", accept = ".csv"),
                             p("File must be csv format and max size is 200 MB.")),
                           mainPanel(
                             h4("Preview of uploaded dataset"),
                             DT::dataTableOutput("dataPreview")
                           ))),
    
    tabPanel("Analysis",
             tabsetPanel(
               id = "analysis_tabs",
               tabPanel("Uploaded Data", 
                        h4("Uploaded Data Table (Select multiple rows below)"),
                        actionButton("goToGene", "Analyse Selected Genes"),
                        br(), br(),
                        DT::dataTableOutput("predictionTable")),
               
               tabPanel("Gene Analysis", 
                        h4("Selected Rows"),
                        DT::dataTableOutput("selectedRowTable")),
               
               tabPanel("Real predictions",
                        h4("Model Prediction Results"),
                        DT::dataTableOutput("resultsTable")),
               
               tabPanel("Compare Samples?", "Comments section goes here")
             )
    ),
    
    tabPanel("Evaluation",
             fluidRow(
               column(6, plotOutput("pcaPlot")),
               column(6,
                      h4("Model Performance Metrics"),
                      verbatimTextOutput("accStats")))
    )
  )
)

server <- function(input, output, session) {
  data <- reactive({ req(input$file1); read.csv(input$file1$datapath) })
  selected_rows <- reactiveVal(data.frame())
  
  output$dataPreview <- DT::renderDataTable({ data() })
  output$predictionTable <- DT::renderDataTable({ datatable(data(), selection = "multiple") })
  
  observeEvent(input$goToGene, {
    req(input$predictionTable_rows_selected)
    selected <- data()[input$predictionTable_rows_selected, , drop = FALSE]
    selected_rows(selected)
    updateTabsetPanel(session, "main_tabs", selected = "Analysis")
    updateTabsetPanel(session, "analysis_tabs", selected = "Gene Analysis")
  })
  
  output$selectedRowTable <- DT::renderDataTable({ req(nrow(selected_rows()) > 0); selected_rows() })
  output$resultsTable <- DT::renderDataTable({ datatable(results_df) })
  output$pcaPlot <- renderPlot({
    ggplot(pca_df, aes(x = PC1, y = PC2, color = Grade)) +
      geom_point(size = 3) +
      theme_bw() +
      labs(title = "PCA of Samples by Grade")
  })
  output$accStats <- renderPrint({
    cat("Test set Accuracy:", round(accuracy_val, 3), "\n")
    cat("Test set Kappa:", round(kappa_val, 3), "\n")
  })
}

shinyApp(ui = ui, server = server)

