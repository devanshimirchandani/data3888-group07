cat(">>> Sourcing ui.R <<<\n")

# Test Preva
library(shiny)
library(GEOquery)
library(DT)
library(ggplot2)
library(bslib)
library(reshape2)
library(dplyr)
library(pheatmap)
library(limma)
library(caret)
library(MLmetrics)
library(caret)
library(MLmetrics)
library(pROC)
library(bslib)

ui <- fluidPage(
  theme = bs_theme(bg = "#fcf8f3", fg = "#2d2d2d", primary = "#e3a5b2"),
  tags$head(tags$style(HTML("
    .topbar-container {
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding: 10px 30px;
      background-color: #fcf8f3;
      border-bottom: 2px solid #e3a5b2;
    }
    .logo-title {
      height: 70px;
    }
    .toolbar-group {
      display: flex;
      gap: 24px;
      margin-left: 30px;
    }
    .toolbar-group .nav-tab {
      background: none !important;
      border: none !important;
      font-size: 17px;
      font-weight: 500;
      padding: 10px 0;
      color: #555;
    }
    .toolbar-group .nav-tab:hover {
      color: #e3a5b2 !important;
      border-bottom: 2px solid #e3a5b2;
    }
    .toolbar-group .nav-tab.active {
      color: black !important;
      font-weight: bold;
      border-bottom: 2px solid black;
    }
    .dataTable tbody tr.selected {
      background-color: #f3d4dc !important;
    }
  "))),
  
  div(class = "topbar-container",
      img(src = "Preva.jpg", class = "logo-title"),
      div(class = "toolbar-group",
          actionButton("about_tab", "About", class = "nav-tab"),
          actionButton("how_tab", "How to Use", class = "nav-tab"),
          actionButton("data_tab", "Data", class = "nav-tab"),
          actionButton("analysis_tab", "Analysis", class = "nav-tab"),
          actionButton("eval_tab", "Evaluation", class = "nav-tab")
      )
  ),
  
  uiOutput("conditional_selectors"),
  
  mainPanel(uiOutput("main_ui"), width = 12)
)





prepare_comparison_data <- function() {
  gse <- getGEO("GSE15852", GSEMatrix = TRUE)[[1]]
  M <- log2(exprs(gse) + 1)
  p <- pData(gse)
  grades <- factor(p$`grade:ch1`, levels = c("normal", "grade 1", "grade 2", "grade 3"))
  
  set.seed(123)
  train_idx <- createDataPartition(grades, p = 0.8, list = FALSE)
  expr_train <- M[, train_idx]
  labels_train <- make.names(grades[train_idx])
  
  labels_factor <- factor(labels_train, levels = c("normal", "grade.1", "grade.2", "grade.3"))
  design <- model.matrix(~ 0 + labels_factor)
  colnames(design) <- levels(labels_factor)
  fit <- lmFit(expr_train, design)
  fit <- eBayes(fit)
  res <- topTable(fit, number = Inf, adjust.method = "BH")
  
  sig_genes <- rownames(subset(res, adj.P.Val < 0.05))
  group_means <- sapply(levels(labels_factor), function(lv) {
    rowMeans(expr_train[, labels_train == lv, drop = FALSE])
  })
  expr_diff <- apply(group_means, 1, function(x) max(x) - min(x))
  sig_genes_strong <- sig_genes[expr_diff[sig_genes] > 1]
  res_sig <- res[sig_genes_strong, ]
  res_sig <- res_sig[order(res_sig$adj.P.Val), ]
  
  top_gene_sets <- list(
    top30 = rownames(res_sig)[1:30],
    top50 = rownames(res_sig)[1:50],
    top80 = rownames(res_sig)[1:80],
    top110 = rownames(res_sig)[1:110]
  )
  
  results <- list()
  cv_folds <- c(5, 10)
  
  for (gene_set_name in names(top_gene_sets)) {
    gene_set <- top_gene_sets[[gene_set_name]]
    x_train <- t(expr_train[gene_set, , drop = FALSE])
    y_train <- factor(labels_train, levels = c("normal", "grade.1", "grade.2", "grade.3"))
    
    results[[gene_set_name]] <- list()
    
    for (cv in cv_folds) {
      trctrl <- trainControl(method = "cv", number = cv, savePredictions = "final")
      
      svm_mod <- train(x = x_train, y = y_train, method = "svmLinear", trControl = trctrl)
      knn_mod <- train(x = x_train, y = y_train, method = "knn", trControl = trctrl,
                       tuneGrid = expand.grid(k = 5:10))
      rf_mod <- train(x = x_train, y = y_train, method = "rf", trControl = trctrl, ntree = 100)
      
      results[[gene_set_name]][[paste0("CV", cv)]] <- list(
        SVM = svm_mod$resample$Accuracy,
        KNN = knn_mod$resample$Accuracy,
        RF = rf_mod$resample$Accuracy
      )
    }
  }
  
  results
}
