# Preva

options(shiny.maxRequestSize = 100 * 1024^2) 

library(shiny)
library(shinyjs)
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
library(tidyr)
library(rintrojs)


load_data <- function(dataset = "GSE15852") {
  df <- data.frame(Patient = character(), Grade = factor(levels = c("0", "1", "2", "3")))
  
  
  if (dataset == "GSE15852") {
    print("Loading GSE15852.RData")
    load("GSE15852.RData")  
    df <- as.data.frame(t(exprs_15852))
    df$Patient <- rownames(df)
    if (!is.null(pdata_15852$grade)) {
      df$Grade <- pdata_15852$grade[match(df$Patient, rownames(pdata_15852))]
    }
    if ("histopathological exam:ch1" %in% colnames(pdata_15852)) {
      df$Histopathology <- pdata_15852$`histopathological exam:ch1`[match(df$Patient, rownames(pdata_15852))]
    }
    
    return(df)
  }
  else if (dataset == "GSE10810") {
    load("data_GSE10810.RData")
    df <- as.data.frame(t(expr_10810))
    df$Patient <- rownames(df)
    if (exists("pdata_10810") && "grade" %in% colnames(pdata_10810)) {
      df$Grade <- pdata_10810$grade[match(df$Patient, rownames(pdata_10810))]
    }
    
  } else if (dataset == "GSE17907") {
    load("data_GSE17907.RData")
    df <- as.data.frame(t(expr_17907))
    df$Patient <- rownames(df)
    if (exists("pdata_17907") && "grade" %in% colnames(pdata_17907)) {
      df$Grade <- pdata_17907$grade[match(df$Patient, rownames(pdata_17907))]
    }
    
  } else if (dataset == "Combined") {
    load("data_Combined.RData")  # contains expr_corrected and metadata
    df <- as.data.frame(t(expr_corrected))
    df$Patient <- rownames(df)
    
    if (exists("metadata") && all(c("sample", "grade") %in% colnames(metadata))) {
      metadata_use <- metadata[metadata$sample %in% rownames(df), ]
      
      if (nrow(metadata_use) > 0) {
        df <- merge(df, metadata_use[, c("sample", "grade")],
                    by.x = "Patient", by.y = "sample", all.x = TRUE)
        
        if ("grade" %in% colnames(df)) {
          colnames(df)[which(colnames(df) == "grade")] <- "Grade"
        }
      } else {
        warning("No matching samples between metadata and expression matrix.")
        df$Grade <- NA
      }
    } else {
      df$Grade <- NA
    }
  }
  
  
  if (!"Grade" %in% colnames(df)) {
    if (nrow(df) == 0) {
      df$Grade <- factor(levels = c("0", "1", "2", "3"))
    } else {
      df$Grade <- rep(NA, nrow(df))
    }
  }
  
  return(df)
}

make_result <- function(pred, prob, labels_test, x_test) {
  rownames(prob) <- rownames(x_test)
  df <- data.frame(
    Patient = rownames(prob),
    TrueGrade = gsub("X", "", as.character(labels_test)),
    PredictedGrade = gsub("X", "", as.character(pred)),
    `Normal Chance` = paste0(round(prob[, "X0"] * 100, 1), "%"),
    `Grade 1 Chance` = paste0(round(prob[, "X1"] * 100, 1), "%"),
    `Grade 2 Chance` = paste0(round(prob[, "X2"] * 100, 1), "%"),
    `Grade 3 Chance` = paste0(round(prob[, "X3"] * 100, 1), "%")
  )
  return(df)
}

prepare_predictions_all_models <- function() {
  load("GSE15852.RData")
  M <- exprs_15852
  p <- pdata_15852
  
  grade_vec <- tolower(p$grade)
  grade_vec <- gsub("grade ", "", grade_vec)
  grade_vec <- gsub("grade", "", grade_vec)
  grade_vec <- gsub("normal", "0", grade_vec)
  grade_vec <- trimws(grade_vec)
  
  grade_num_all <- as.numeric(grade_vec)
  grade_fac_all <- factor(make.names(grade_vec), levels = c("X0", "X1", "X2", "X3"))
  
  set.seed(123)
  train_idx <- createDataPartition(grade_fac_all, p = 0.8, list = FALSE)
  expr_train <- M[, train_idx]
  expr_test  <- M[, -train_idx]
  
  grade_num <- grade_num_all[train_idx]
  design <- model.matrix(~ grade_num)
  
  labels_train <- grade_fac_all[train_idx]
  labels_test  <- grade_fac_all[-train_idx]
  
  fit <- lmFit(expr_train, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
  sig_genes <- rownames(subset(res, adj.P.Val < 0.05 & abs(logFC) > 1))
  
  x_train <- t(expr_train[sig_genes, , drop = FALSE])
  x_test  <- t(expr_test[sig_genes, , drop = FALSE])
  rownames(x_test) <- colnames(expr_test)
  
  trctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                         summaryFunction = multiClassSummary, savePredictions = "final")
  
  svm_mod <- train(x = x_train, y = labels_train, method = "svmLinear", trControl = trctrl)
  knn_mod <- train(x = x_train, y = labels_train, method = "knn", trControl = trctrl, tuneGrid = expand.grid(k = 5:10))
  rf_mod  <- train(x = x_train, y = labels_train, method = "rf",  trControl = trctrl)
  
  pred_svm <- predict(svm_mod, newdata = x_test)
  pred_knn <- predict(knn_mod, newdata = x_test)
  pred_rf  <- predict(rf_mod,  newdata = x_test)
  
  prob_svm <- predict(svm_mod, newdata = x_test, type = "prob")
  prob_knn <- predict(knn_mod, newdata = x_test, type = "prob")
  prob_rf  <- predict(rf_mod,  newdata = x_test, type = "prob")
  
  list(
    result = list(
      SVM = make_result(pred_svm, prob_svm, labels_test, x_test),
      KNN = make_result(pred_knn, prob_knn, labels_test, x_test),
      RF  = make_result(pred_rf,  prob_rf,  labels_test, x_test)
    ),
    models = list(
      SVM = svm_mod,
      KNN = knn_mod,
      RF  = rf_mod
    )
  )
}




prepare_predictions_for <- function(dataset_name) {
  if (dataset_name == "GSE10810") {
    load("data_GSE10810.RData")
    M <- expr_10810
    p <- pdata_10810
  } else if (dataset_name == "GSE17907") {
    load("data_GSE17907.RData")
    M <- expr_17907
    p <- pdata_17907
  } else if (dataset_name == "Combined") {
    load("data_Combined.RData")
    M <- expr_corrected
    p <- metadata
  } else {
    stop("Unsupported dataset")
  }
  
  grade_vec <- tolower(p$grade)
  grade_vec <- gsub("grade ", "", grade_vec)
  grade_vec <- gsub("grade", "", grade_vec)
  grade_vec <- gsub("normal", "0", grade_vec)
  grade_vec <- trimws(grade_vec)
  
  grade_fac_all <- factor(make.names(grade_vec), levels = c("X0", "X1", "X2", "X3"))
  
  set.seed(123)
  train_idx <- createDataPartition(grade_fac_all, p = 0.8, list = FALSE)
  expr_train <- M[, train_idx]
  expr_test <- M[, -train_idx]
  labels_train <- grade_fac_all[train_idx]
  labels_test <- grade_fac_all[-train_idx]
  
  grade_num <- as.numeric(as.character(factor(labels_train, labels = 0:3)))
  design <- model.matrix(~ grade_num)
  fit <- lmFit(expr_train, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
  sig_genes <- rownames(subset(res, adj.P.Val < 0.05 & abs(logFC) > 1))
  
  x_train <- t(expr_train[sig_genes, , drop = FALSE])
  x_test  <- t(expr_test[sig_genes, , drop = FALSE])
  rownames(x_test) <- colnames(expr_test)
  
  trctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                         summaryFunction = multiClassSummary, savePredictions = "final")
  
  svm_mod <- train(x = x_train, y = labels_train, method = "svmLinear", trControl = trctrl)
  knn_mod <- train(x = x_train, y = labels_train, method = "knn", trControl = trctrl,
                   tuneGrid = expand.grid(k = 5:10))
  rf_mod  <- train(x = x_train, y = labels_train, method = "rf", trControl = trctrl)
  
  pred_svm <- predict(svm_mod, newdata = x_test)
  pred_knn <- predict(knn_mod, newdata = x_test)
  pred_rf  <- predict(rf_mod,  newdata = x_test)
  
  prob_svm <- predict(svm_mod, newdata = x_test, type = "prob")
  prob_knn <- predict(knn_mod, newdata = x_test, type = "prob")
  prob_rf  <- predict(rf_mod,  newdata = x_test, type = "prob")
  
  list(
    result = list(
      SVM = make_result(pred_svm, prob_svm, labels_test, x_test),
      KNN = make_result(pred_knn, prob_knn, labels_test, x_test),
      RF  = make_result(pred_rf,  prob_rf,  labels_test, x_test)
    ),
    models = list(
      SVM = svm_mod,
      KNN = knn_mod,
      RF  = rf_mod
    )
  )
}






tmp <- prepare_predictions_all_models()
all_results_df <- list(GSE15852 = tmp$result)
all_model_list <- list(GSE15852 = tmp$models)

for (ds in c("GSE10810", "GSE17907", "Combined")) {
  tmp <- prepare_predictions_for(ds)
  all_results_df[[ds]] <- tmp$result
  all_model_list[[ds]] <- tmp$models
}




summary_df <- readRDS("model_summary_df.rds")


predict_uploaded_data <- function(uploaded_expr, model_name = "SVM", dataset = "GSE15852") {
  mod <- all_model_list[[dataset]][[model_name]]
  
  if (is.null(mod$trainingData)) {
    showNotification("Model missing training data; cannot extract gene names.", type = "error")
    return(data.frame())
  }
  
  train_genes <- colnames(mod$trainingData)[-ncol(mod$trainingData)]
  
  if (is.null(colnames(uploaded_expr))) {
    showNotification("Uploaded data has no column names.", type = "error")
    return(data.frame())
  }
  
  matching_genes <- intersect(colnames(uploaded_expr), train_genes)
  
  if (length(matching_genes) == 0) {
    showNotification("No matching genes found between uploaded data and model!", type = "error")
    return(data.frame())
  }
  
  test_data <- tryCatch({
    tmp <- uploaded_expr[, matching_genes, drop = FALSE]
    tmp <- tmp[, train_genes, drop = FALSE]  # fill with NA for missing
    tmp
  }, error = function(e) {
    showNotification("Prediction failed due to invalid test data.", type = "error")
    return(NULL)
  })
  
  if (is.null(test_data)) return(data.frame())
  
  pred <- tryCatch({
    predict(mod, newdata = test_data)
  }, error = function(e) {
    showNotification(paste("Prediction failed:", conditionMessage(e)), type = "error")
    return(NULL)
  })
  
  prob <- tryCatch({
    predict(mod, newdata = test_data, type = "prob")
  }, error = function(e) {
    showNotification(paste("Probability prediction failed:", conditionMessage(e)), type = "error")
    return(NULL)
  })
  
  if (is.null(pred) || is.null(prob)) return(data.frame())
  
  df <- data.frame(
    Patient = rownames(test_data),
    PredictedGrade = gsub("X", "", as.character(pred))
  )
  
  for (col in colnames(prob)) {
    df[[paste0(gsub("X", "Grade.", col), ".Chance")]] <- paste0(round(prob[[col]] * 100, 1), "%")
  }
  
  return(df)
}







comparison_results <- reactiveVal(NULL)



ui <- fluidPage(
  introjsUI(), 
  theme = bs_theme(bg = "#fcf8f3", fg = "#2d2d2d", primary = "#e3a5b2"),
  tags$head(tags$style(HTML("
  
  body {
    margin-top: 10px !important;
  }
  

  
  .topbar-container {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 16px 30px 12px 30px;
    background-color: #fcf8f3;
    border-bottom: 2px solid #e3a5b2;
    box-shadow: 0 2px 6px rgba(0,0,0,0.05);
  }

  .logo-title {
    height: 90px;
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
    transition: border-bottom 0.3s ease, color 0.3s ease;
    position: relative;
  }

  .toolbar-group .nav-tab:hover {
    color: #e3a5b2 !important;
    border-bottom: 2px solid #e3a5b2;
  }

  .toolbar-group .nav-tab.active {
    color: #d07095 !important;
    font-weight: bold;
  }

  .toolbar-group .nav-tab.active::after {
    content: '';
    position: absolute;
    bottom: -4px;
    left: 0;
    width: 100%;
    height: 3px;
    background-color: #d07095;
    border-radius: 2px;
  }

    .introjs-tooltip {
    border: 2px solid #e3a5b2 !important;
    box-shadow: 0px 0px 10px rgba(227, 165, 178, 0.3) !important;
  }
  .introjs-tooltiptext, .introjs-helperNumberLayer {
    color: #e3a5b2 !important;
  }
  .introjs-button {
    background-color: #e3a5b2 !important;
    border: none !important;
  }
  
  .dataTable tbody tr.selected {
    background-color: #f3d4dc !important;
  }

  h2 {
    font-size: 28px;
    font-weight: 700;
    margin-top: 10px !important;
    margin-bottom: 10px;
  }

  h3 {
    margin-top: 15px !important;
  }

  h1 {
    margin-top: 25px !important;
  }
  
  .nav-tabs > li > a {
    color: #e3a5b2 !important;
    font-weight: 500;
    padding: 6px 16px;
    margin: 6px 6px 0 0;
    font-size: 17px;
    border-radius: 999px !important;
    border: none !important;
    background-color: transparent !important;
    transition: all 0.3s ease;
    box-shadow: none !important;
  }

  .nav-tabs > li > a:hover {
    background-color: #fce9ef !important;
    border-color: #e3a5b2 !important;
    color: #d07095 !important;
    box-shadow: 0 1px 2px rgba(227,165,178,0.25) !important;
  }

  .nav-tabs > li.active > a,
  .nav-tabs > li.active > a:focus,
  .nav-tabs > li.active > a:hover {
    background-color: #e3a5b2 !important;
    color: white !important;
    border: 1px solid #e3a5b2 !important;
    font-weight: bold;
    box-shadow: inset 0 -2px 0 #c66690;
  }

  .nav-tabs {
    border-bottom: none !important;
  }

  .tab-content {
    background-color: transparent !important;
    border: none !important;
    padding-top: 10px;
    margin-top: 0px;
  }
    .selector-row {
    margin-bottom: 0px;  
    }

"))),
  
  div(class = "topbar-container",
      
      div(style = "flex: 1; display: flex; align-items: center; justify-content: flex-start;",
          img(src = "Preva.png", class = "logo-title",
              `data-step` = 1,
              `data-intro` = "Welcome to Preva! This logo reminds you you're in a gene expression prediction tool.")
      ),
      
      div(class = "toolbar-group", style = "flex: 5; display: flex; justify-content: center; gap: 24px;",
          
          actionButton("about_tab", "About Preva", class = "nav-tab",
                       `data-step` = 2,
                       `data-intro` = "Click here to learn what Preva is about."),
          actionButton("intro_tab", "Introduction", class = "nav-tab",
                       `data-step` = 3,
                       `data-intro` = "This section introduces you to breast cancer and gene data."),
          actionButton("data_tab", "Data Processing", class = "nav-tab",
                       `data-step` = 4,
                       `data-intro` = "Select this tab to upload, browse, and prepare your datasets."),
          actionButton("analysis_tab", "Predictive Modelling", class = "nav-tab",
                       `data-step` = 5,
                       `data-intro` = "This is where you build machine learning models."),
          actionButton("eval_tab", "Evaluation", class = "nav-tab",
                       `data-step` = 6,
                       `data-intro` = "Finally, check how well your model performed!")
      ),
      
      div(style = "flex: 1; display: flex; justify-content: flex-end; align-items: center;",
          actionButton("help_modal_button", label = NULL, icon = icon("question-circle"),
                       style = "background: none; border: none; color: #e3a5b2; font-size: 22px;",
                       class = "nav-tab")
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







server <- function(input, output, session) {
  
  
  uploaded_data <- reactiveVal(NULL)
  
  uploaded_exprs <- reactiveVal(NULL)
  
  overview_df <- reactiveVal(NULL)
  

  
  observeEvent(input$convert_format, {
    req(uploaded_data())
    first_obj <- uploaded_data()
    
    if (is.data.frame(first_obj)) {
      df_tmp <- as.data.frame(first_obj)
      numeric_cols <- sapply(df_tmp, is.numeric)
      
      if (any(numeric_cols)) {
        if (is.data.frame(first_obj)) {
          numeric_cols <- sapply(first_obj, is.numeric)
          expr_try <- as.matrix(first_obj[numeric_cols])  
          if (nrow(expr_try) > 0 && ncol(expr_try) > 0) {
            expr_try <- t(expr_try)  
            uploaded_exprs(expr_try)
            

            updateSelectInput(session, "patient_col", choices = rownames(expr_try), selected = rownames(expr_try)[1])
            updateSelectInput(session, "grade_col", choices = c("None")) 
            showNotification("Conversion complete. Please manually select Patient column if needed.", type = "message")
            return()
          }
        }
        
        if (nrow(expr_try) > 0 && ncol(expr_try) > 0) {
          uploaded_exprs(expr_try)
          updateSelectInput(session, "patient_col", choices = colnames(pdata_try))
          updateSelectInput(session, "grade_col", choices = c("None", colnames(pdata_try)))
          showNotification("Auto-conversion successful. Please verify patient/grade columns.", type = "message")
          return()
        }
      }
    }
    
    showNotification("Conversion failed: no numeric expression matrix found.", type = "error")})
    
  
  
  data <- reactive({
    if (input$dataset_choice == "Uploaded" && !is.null(uploaded_exprs())) {
      expr_df <- uploaded_exprs()
      df <- as.data.frame(t(expr_df))
      df$Patient <- rownames(df)
      
      if (!is.null(uploaded_data()) && !is.null(input$pdata_object) && !is.null(input$grade_col)) {
        pdata <- uploaded_data()[[input$pdata_object]]
        if (is.data.frame(pdata) && input$grade_col %in% colnames(pdata)) {
          grade_info <- pdata[, c(input$grade_col), drop = FALSE]
          rownames(grade_info) <- rownames(pdata)
          colnames(grade_info) <- "Grade"
          df <- merge(df, grade_info, by.x = "Patient", by.y = "row.names", all.x = TRUE)
        }
      }
      
      if ("Grade" %in% colnames(df)) {
        grade_vec <- tolower(df$Grade)
        grade_vec <- gsub("grade ", "", grade_vec)
        grade_vec <- gsub("grade", "", grade_vec)
        grade_vec <- gsub("normal", "0", grade_vec)
        grade_vec <- trimws(grade_vec)
        df$Grade <- grade_vec
      }
      return(df)
    }
    
    dataset <- input$dataset_choice
    if (is.null(dataset) || dataset == "") return(data.frame())
    load_data(dataset)
  })
  
  
  # Uploaded_data
  
  uploaded_data <- reactiveVal(NULL)
  


  observeEvent(input$browse_data, {
    inFile <- input$browse_data
    if (is.null(inFile)) return(NULL)
    
    e <- new.env()
    load(inFile$datapath, envir = e)
    uploaded_data(e)
    
    obj_names <- ls(e)
    

    updateSelectInput(session, "expr_object", choices = obj_names)
    updateSelectInput(session, "pdata_object", choices = obj_names)

    pdata_candidates <- grep("^pdata", obj_names, value = TRUE)
    if (length(pdata_candidates) > 0) {
      pdata <- e[[pdata_candidates[1]]]
      if (is.data.frame(pdata)) {
        updateSelectInput(session, "patient_col", choices = colnames(pdata))
        updateSelectInput(session, "grade_col", choices = c("None", colnames(pdata)))
      }
    }
    
    showNotification("Upload complete. Please confirm expression + phenotype objects.", type = "message")
  })
  
  
  

  
  observeEvent(input$pdata_object, {
    req(uploaded_data())
    req(input$pdata_object != "") 
    pdata <- uploaded_data()[[input$pdata_object]]
    
    if (input$pdata_object %in% names(uploaded_data())) {
      pdata <- uploaded_data()[[input$pdata_object]]
      
      if (is.data.frame(pdata)) {
        updateSelectInput(session, "patient_col", choices = colnames(pdata))
        updateSelectInput(session, "grade_col", choices = c("None", colnames(pdata)))
      }
    }
  })
  
  
  observeEvent(input$check_overview, {
    req(input$expr_object)
    expr <- uploaded_data()[[input$expr_object]]
    if (!is.matrix(expr) && !is.data.frame(expr)) {
      showModal(modalDialog("Invalid expression matrix!", easyClose = TRUE))
      return()
    }
    
    expr <- as.data.frame(t(expr))
    expr$Patient <- rownames(expr)
    
    if (!is.null(input$pdata_object) && input$grade_col != "None" && input$patient_col != "") {
      pdata <- uploaded_data()[[input$pdata_object]]
      if (is.data.frame(pdata)) {
        if (input$patient_col %in% colnames(pdata)) {
          rownames(pdata) <- pdata[[input$patient_col]]
        }
        if (input$grade_col %in% colnames(pdata)) {
          pdata_use <- pdata[, input$grade_col, drop = FALSE]
          colnames(pdata_use) <- "Grade"
          expr <- merge(expr, pdata_use, by.x = "Patient", by.y = "row.names", all.x = TRUE)
        }
      }
    }
    
    rownames(expr) <- expr$Patient
    expr$Patient <- NULL
    uploaded_exprs(expr)  
    
    cols_to_show <- intersect(c("Patient", "Grade"), colnames(expr))
    preview_df <- expr[, cols_to_show, drop = FALSE]
    preview_df$Patient <- rownames(expr)
    
    output$overview_table <- renderDT({
      datatable(preview_df, options = list(pageLength = 10))
    })
    
    showModal(modalDialog(
      title = "Dataset Overview",
      DTOutput("overview_table"),
      easyClose = TRUE
    ))
    
    showNotification("Dataset has been loaded and is ready for prediction!", type = "message")
  })
  
  
  
  
  
  
  
  
  observe({
    ds <- input$dataset_choice
    model <- input$model_choice
    
    if (is.null(ds) || is.null(model)) return()
    
    if (ds == "Uploaded" && !is.null(uploaded_exprs())) {
      df <- predict_uploaded_data(uploaded_exprs(), model)
    } else {
      if (!(ds %in% names(all_results_df)) || !(model %in% names(all_results_df[[ds]]))) {
        updateCheckboxGroupInput(session, "selected_genes", choices = character(0))
        return()
      }
      df <- all_results_df[[ds]][[model]]
    }
    
    if (is.null(df) || ncol(df) < 6) {
      updateCheckboxGroupInput(session, "selected_genes", choices = character(0))
      return()
    }
    

    gene_cols <- colnames(df)[6:ncol(df)]
    gene_cols <- gene_cols[!grepl("Chance", gene_cols, ignore.case = TRUE)]
    
    updateCheckboxGroupInput(session, "selected_genes", choices = gene_cols)
  })
  
  
  
  
  
  
  page <- reactiveVal("About")
  
  observeEvent(input$about_tab, page("About"))
  observeEvent(input$intro_tab, page("Introduction"))
  observeEvent(input$data_tab, page("DataProcessing"))
  observeEvent(input$analysis_tab, page("Predictive Modelling"))
  observeEvent(input$eval_tab, page("Evaluation"))
  observeEvent(input$findings_tab, page("Findings"))
  

  observeEvent(input$evaluation_tab, {
    if (input$evaluation_tab == "Model Comparison" &&
        is.null(comparison_results())) {
      
      showModal(modalDialog(
        title = "Preparing model comparison data...",
        "This may take a moment as multiple models are being trained. Please wait.",
        easyClose = FALSE,
        footer = NULL
      ))
      
      result <- prepare_comparison_data()
      comparison_results(result)
      saveRDS(result, "comparison_cache.rds")
      
      removeModal()
    }
  })
  
  
  
  observe({
    if (is.null(input$tour_shown)) {
      showModal(modalDialog(
        title = div(style = "color: #e3a5b2; font-size: 24px; font-weight: bold;",
                    "🎀 Welcome to Preva"),
        div(style = "font-size: 15px;",
            "This guided tour will help you explore the key parts of the Preva application, from selecting datasets to evaluating predictive models."),
        easyClose = FALSE,
        footer = actionButton("start_full_tour", "Start Tour", class = "btn btn-primary", 
                              style = "background-color: #e3a5b2; border-color: #e3a5b2;")
      ))
    }
  })
  
  observeEvent(input$start_full_tour, {
    removeModal()
    introjs(session, options = list(
      "nextLabel" = "Next →",
      "prevLabel" = "← Back",
      "doneLabel" = "Finish",
      "showStepNumbers" = TRUE,
      "showBullets" = TRUE,
      "showProgress" = TRUE,
      "exitOnOverlayClick" = TRUE,
      "exitOnEsc" = TRUE
    ))
  })
  
  observeEvent(input$tour_button, {
    introjs(session, options = list(
      "nextLabel" = "Next →",
      "prevLabel" = "← Back",
      "doneLabel" = "Finish",
      "showStepNumbers" = TRUE,
      "exitOnOverlayClick" = TRUE,
      "showBullets" = TRUE,
      "showButtons" = TRUE
    ))
  })
  
  

  observeEvent(input$help_modal_button, {
    showModal(modalDialog(
      title = tags$span(style = "color: #e3a5b2; font-size: 20px; font-weight: bold;", "🔍 Help with Preva"),
      easyClose = TRUE,
      footer = modalButton("Close"),
      HTML("
      <p>This is the <b>Preva</b> application — an interactive tool for exploring gene expression data and predicting breast cancer grades using machine learning.</p>
      <ul>
        <li><b>Data Processing:</b> Upload and explore datasets.</li>
        <li><b>Predictive Modelling:</b> Select models and make predictions.</li>
        <li><b>Evaluation:</b> Assess performance using PCA, ROC, F1, etc.</li>
        <li><b>Tour:</b> Click the question icon at the top to take a guided tour.</li>
      </ul>
    ")
    ))
  })
  
  
  
  
  output$conditional_selectors <- renderUI({
    if (page() %in% c("Predictive Modelling", "Evaluation")) {
      div(class = "selector-row",
        fluidRow(
          column(4,
                 selectInput("dataset_choice", "Dataset:",
                             choices = c("GSE15852", "GSE10810", "GSE17907", "Combined (GSE10810 and GSE17907)", "Uploaded"),
                             selected = "GSE15852")
                 
          ),
          column(4,
                 selectInput("model_choice", "Model:",
                               choices = c("SVM", "KNN", "RF"),
                               selected = "SVM")
               
        )
       )
      )
    }
  })
  
  output$main_ui <- renderUI({
    req(input$dataset_choice)
    
    switch(page(),
           "About" = {
             tagList(
               h2("About Preva"),
               p("Breast cancer is the most common cancer among females, and one in seven women are diagnosed with it in their lifetime. This is why we created Preva. This is a tool that medical researchers can use to enhance knowledge that facilitates early breast cancer detection, and understanding the interaction between genes and breast cancer diagnoses. But also to detect it when it matters, and know what stage the cancer is at for further action."),
               p("This Shiny web application uses machine learning models to help predict and visualise breast cancer diagnoses. The application, which was developed using R and is driven by interactive data visualisation, gives users an easy-to-use way to examine patient data and prediction outcomes."),
               tags$b("Sources:"),
               tags$ul(tags$li(a("National Breast Cancer Foundation (2021)", href = "https://nbcf.org.au/about-breast-cancer/breast-cancer-stats/", target = "_blank")))
             )
           },
          
           "Introduction" = {
             tabsetPanel(
               tabPanel("Breast Cancer",
                        h2("Why Focus on Breast Cancer?"),
                        p("Coming soon")  
               ),
               tabPanel("Gene Expression Data",
                        h2("Understanding Gene Expression"),
                        p("Coming soon")
               )
             )
           },
           
           "DataProcessing" = {
             tabsetPanel(
               tabPanel("About Datasets",
                        h2("About the Datasets"),
                        p("This application supports several public gene expression datasets relevant to breast cancer research."),
                        tags$ul(
                          tags$li(strong("GSE15852:"), " Coming soon"),
                          tags$li(strong("GSE10810:"), " Coming soon"),
                          tags$li(strong("GSE17907:"), " Coming soon"),
                          tags$li(strong("Combined:"), " A merged dataset combining GSE10810 and GSE17907.")
                        )
               ),
               tabPanel("Browse Dataset",
                        fluidRow(
                          column(6,
                                 fileInput("browse_data", "Upload your own dataset (.RData)", accept = ".RData"),
                                 selectInput("expr_object", "Select expression matrix object (genes × samples):", choices = NULL, selected = ""),
                                 selectInput("pdata_object", "Select phenotype (metadata) object:", choices = NULL, selected = ""),
                                 actionButton("check_overview", "Check Dataset Overview"), 
                                 DTOutput("browse_table")
                          ),
                          column(6,
                                 selectInput("patient_col", "Select Patient Column", choices = NULL),
                                 selectInput("grade_col", "Select Grade Column (optional)", choices = NULL)
                          )
                        )
               ),
               tabPanel("Dataset Overview",
                        selectInput("dataset_choice_data", "Select Dataset:",
                                    choices = c("GSE15852", "GSE10810", "GSE17907", "Combined (GSE10810 and GSE17907)"),
                                    selected = "GSE15852"),
                        DTOutput("data_table"))
               
             )
           },
           
           "Predictive Modelling" = {tabsetPanel(
             id = "modelling_tab",
             tabPanel("Model Building",
                      h2("Model:"),
                      p("SVM:
                        KNN: 
                        RF:")  
             ),
             tabPanel("Prediction", 
                      h2("Prediction Table"),
                      fluidRow(
                        column(6, actionButton("select_all", "Select All")),
                        column(6, actionButton("clear_all", "Clear Selection"))
                      ),
                      DTOutput("prediction_table")
             ),
             
             
             tabPanel("Gene Analysis",
                      h3("Select Genes to Visualize"),
                      p("This section allows you to explore how different genes are expressed across breast cancer grades. ",
                        "By selecting the top variable genes, you can generate interactive boxplots that show how gene expression ",
                        "levels differ between groups (e.g., normal, grade 1, grade 2, grade 3)."),
                      p("You don’t need to write any code — just use the slider to adjust how many genes you’d like to view, ",
                        "and the app will automatically select the most variable ones (those that show the biggest changes across patients)."),
                      p(strong("Tip:"), " If you're not sure what to select, just start with the default setting and explore. ",
                        "The plots will update based on your choices!"),
                      checkboxGroupInput("selected_genes", "Choose Genes:", choices = NULL),
                      sliderInput("top_gene_slider", "Top Variable Genes:", min = 5, max = 50, value = 10, step = 5),
                      plotOutput("gene_box", height = "400px")
             )
             ,
             
             tabPanel("Performance",
                      checkboxInput("use_all_data", "Use full dataset (ignore selection)", value = FALSE),
                      
                      h3("Grade Distribution"),
                      p("This section helps you evaluate how well the model predicted breast cancer grades."),
                      
                      p("You can choose between two visualisation styles:",
                        tags$ul(
                          tags$li(strong("Bar + Line Plot:"), " Bars show predicted grade counts, and the line represents the true grade counts, making it easy to compare them directly."),
                          tags$li(strong("Confusion matrix heatmap:"), " Displays true vs predicted grades in a grid format, helping you spot correct and incorrect predictions clearly.")
                        ),
                        "These plots provide an intuitive and visually distinct summary of model performance."
                      ),
                      
                      radioButtons("bar_plot_style", "Plot Style:",
                                   choices = c("Bar + Line Plot" = "facet", "Confusion Matrix Heatmap" = "heatmap"),
                                   inline = TRUE),
                      
                      div(style = "display: flex; justify-content: center;",
                          plotOutput("bar_plot", width = "600px", height = "350px")),
                      
                      
                      h3("Heatmap"),
                      p("This heatmap visualises the expression levels of the top 10 most variable genes across selected patients."),
                      p("Each row represents a gene, and each column represents a patient. The colours reflect gene expression intensity, ",
                        "helping you identify patterns and differences between samples. It can be useful for spotting clusters of patients ",
                        "with similar expression profiles."),
                      div(style = "display: flex; justify-content: center;",
                        plotOutput("heatmap_plot", width = "750px", height = "400px"))
             )
           )
           },
           
           
           "Evaluation" = {
             tabsetPanel(
               id = "evaluation_tab",
               tabPanel("PCA",
                        h3("PCA Plot"),
                        p("This PCA plot allows you to explore how samples group by grade using selected gene expression profiles. ",
                          "You can choose to visualise the sample structure using genes that vary most (top variable genes), or the genes actually used by the selected model."),
                        
                        fluidRow(
                          column(4,
                                 radioButtons("pca_source", "Use genes from:",
                                              choices = c("Top variable genes" = "top_var", "Model-selected genes" = "model"),
                                              selected = "top_var"
                                 ),
                                 
                                 conditionalPanel(
                                   condition = "input.pca_source == 'top_var'",
                                   tagList(
                                     sliderInput("topN", "Number of top variable genes:", min = 30, max = 110, value = 80, step = 10),
                                     helpText("Top variable genes are selected by ranking all genes based on their variance across samples. ",
                                              "This is a model-independent way to explore natural clustering of samples.")
                                   )
                                 ),
                                 
                                 conditionalPanel(
                                   condition = "input.pca_source == 'model'",
                                   helpText("Model-selected genes are the features used by the machine learning model you selected. ",
                                            "This shows how well the model's selected genes group samples by grade.")
                                 ),
                                 helpText("Ellipses are only shown for grades with at least 3 samples.")
                                 
                              
                          ),
                          
                          column(8,
                                 plotOutput("pca_plot", height = "420px")
                          )
                        )
               ),
               tabPanel("F score",
                        h3("F1 Score Summary"),
                        p("The F1 score combines precision and recall for each class, and is especially useful for evaluating imbalanced datasets. ",
                          "This table shows the per-grade classification metrics for the selected model and dataset."),
                        p(strong("Precision:"), "Out of all predicted positive samples for a class, how many were actually correct."),
                        p(strong("Recall:"), "Out of all actual positive samples for a class, how many were correctly identified."),
                        p(strong("F1 Score:"), "The harmonic mean of precision and recall. A balanced measure of accuracy for each class."),
                        DT::dataTableOutput("fscore_table")
               ),
               tabPanel("ROC / AUC",
                        h3("Receiver Operating Characteristic (ROC) Curve"),
                        
                        p("The ROC curve illustrates the performance of a classification model at various threshold settings. ",
                          "It shows the trade-off between sensitivity (true positive rate) and 1 − specificity (false positive rate)."),
                        
                        p("A perfect classifier would reach the top-left corner of the plot (100% sensitivity and 100% specificity). ",
                          "The closer the curve follows this ideal path, the better the model's ability to distinguish between classes."),
                        
                        p(strong("Area Under the Curve (AUC):"), 
                          "The AUC value represents the overall ability of the model to correctly classify positive and negative cases. ",
                          "An AUC of 0.5 means random guessing, while a score close to 1.0 indicates excellent prediction performance."),
                        
                        p("You can select different models and datasets above to view their ROC performance here."),
                        p("This plot shows a separate one-vs-rest ROC curve for each grade class. Each curve evaluates how well the model distinguishes one specific grade from all others."),
                        
                        div(style = "text-align: center;",
                            plotOutput("roc_plot", height = "400px", width = "600px")
                        )
                        
               ),
               tabPanel("Confusion Matrix",
                        h3("Confusion Matrix"),
                        p("This matrix shows how often predicted grades match the actual grades. Diagonal values are correct predictions, and off-diagonal cells represent errors."),
                        # plotOutput("conf_matrix_plot", height = "400px"),
                        h4("False Negative Rate (FNR) Explanation"),
                        p(strong("False Negative Rate (FNR):"), 
                          "Also known as the miss rate, FNR is the proportion of actual positive cases that the model failed to detect. ",
                          "It is computed as 1 − Recall (Sensitivity)."),
                        
                        p("FNR is important in clinical applications, where missing a diagnosis may be more serious than over-predicting one.")
               ),
               tabPanel("Cross Validation",
                        h3("Cross Validation Overview"),
                        p("Cross-validation is a technique used to estimate the generalisation performance of a machine learning model. ",
                          "It works by splitting the data into multiple training/testing folds and averaging performance across folds."),
                        p("This helps prevent overfitting and provides a more robust measure of model accuracy."),
                        p("In future versions, this section will include visualisations of accuracy, F1 score, and other metrics across folds.")
               ),
               tabPanel("Feature Importance (Gene)",
                        h3("Feature Importance Explanation"),
                        p("Feature importance highlights which genes contributed most to the predictions made by the model."),
                        p("This helps identify potentially biologically relevant markers that distinguish between grades of breast cancer."),
                        p("In future versions, this section will include ranked gene lists and barplots showing their relative importance.")
               ),
               tabPanel("Metric Summary",
                        h3("Evaluation Summary Across Gene Sets"),
                        p("This plot shows how model performance changes across different gene sets. Currently shown for demonstration only."),
                        sliderInput("metric_topn", "Top N Genes:", min = 30, max = 110, value = 80, step = 10),
                        plotOutput("summary_plot", height = "400px")
               ),
               tabPanel("Model Comparison",
                        h3("Model Comparison"),
                        p("This section compares the predictive performance of different machine learning models across various gene selection strategies and cross-validation settings."),
                        
                        p("You can choose the number of top genes used for model building and the cross-validation scheme (5-fold or 10-fold) to explore how model performance changes."),
                        
                        p(strong("Note:"), " This comparison uses the selected dataset only. The ",
                          em("Model"), " selector at the top is not used here, as the goal is to evaluate multiple models side-by-side."),
                        
                        fluidRow(
                          column(4,
                                 selectInput("comp_topn", "Select Top N Genes:",
                                             choices = c("top30", "top50", "top80", "top110"),
                                             selected = "top30"),
                                 radioButtons("comp_cv", "CV Fold:",
                                              choices = c("5-fold" = "CV5", "10-fold" = "CV10"),
                                              selected = "CV5")
                          ),
                          column(8,
                                 plotOutput("comparison_plot"),
                                 verbatimTextOutput("mean_accuracy_text")
                          )
                        )
               )
               
             )
           },
           
           
           "Findings" = {
             tagList(
               h2("Key Findings Summary"),
               p("This section will summarize the key takeaways from the modelling and evaluation."),
               h3("Limitations"),
               p("Limitations of current study and modelling process."),
               h3("Future Directions"),
               p("Suggestions for further model validation, independent test datasets, and robustness enhancements.")
             )
           }
           
    )
  })
  
  output$pca_plot <- renderPlot({
    df <- load_data(input$dataset_choice)
    
    if (nrow(df) == 0 || !"Grade" %in% colnames(df) || all(is.na(df$Grade))) {
      plot.new()
      text(0.5, 0.5, "No data available for PCA.", cex = 1.3)
      return()
    }
    
    

    expr_mat <- df[, sapply(df, is.numeric), drop = FALSE]
    grade <- df$Grade
    rownames(expr_mat) <- df$Patient

    if (input$pca_source == "top_var") {
      gene_vars <- apply(expr_mat, 2, var, na.rm = TRUE)
      top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:input$topN])
      x <- expr_mat[, top_genes, drop = FALSE]
    }
    

    else if (input$pca_source == "model") {
      model_obj <- all_model_list[[input$dataset_choice]][[input$model_choice]]
      if (is.null(model_obj) || is.null(model_obj$trainingData)) {
        plot.new()
        text(0.5, 0.5, "Model not available for PCA.", cex = 1.3)
        return()
      }
      
      selected_genes <- colnames(model_obj$trainingData)[-ncol(model_obj$trainingData)]
      selected_genes <- intersect(selected_genes, colnames(expr_mat))
      if (length(selected_genes) < 2) {
        plot.new()
        text(0.5, 0.5, "Not enough model-selected genes found in dataset.", cex = 1.3)
        return()
      }
      
      x <- expr_mat[, selected_genes, drop = FALSE]
    }
    

    pca <- prcomp(x, scale. = TRUE)
    pc_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Grade = grade)
    pc_var <- round(summary(pca)$importance[2, 1:2] * 100, 1)
    
    
    grade_counts <- pc_df %>% count(Grade)
    valid_grades <- grade_counts$Grade[grade_counts$n >= 3]
    
    gg <- ggplot(pc_df, aes(x = PC1, y = PC2, color = Grade)) +
      geom_point(size = 3, alpha = 0.8)
    
    if (length(valid_grades) > 0) {
      gg <- gg + stat_ellipse(data = pc_df %>% filter(Grade %in% valid_grades),
                              type = "norm", linetype = "dashed")
    }
    
    gg <- gg +
      theme_minimal(base_size = 13) +
      labs(
        title = "PCA of Samples (by Grade)",
        x = paste0("PC1 (", pc_var[1], "%)"),
        y = paste0("PC2 (", pc_var[2], "% variance)")
      )
    
    gg
    
  })
  
  
  
  
  
  
  
  
  output$data_table <- renderDT({
    dataset <- input$dataset_choice_data
    if (dataset == "Combined (GSE10810 and GSE17907)") {
      dataset <- "Combined"
    }
    
    df <- load_data(dataset)
    
    target_cols <- c("Patient", "Grade", "Histopathology")
    existing_cols <- intersect(target_cols, colnames(df))
    
    datatable(df[, existing_cols, drop = FALSE],
              selection = "multiple",
              options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10),
              rownames = FALSE)
  })
  
  
  
  output$browse_table <- renderDT({
    env <- uploaded_data()
    req(env)
    req(input$expr_object)
    df <- env[[input$expr_object]]
    if (!is.data.frame(df)) df <- as.data.frame(df)
    datatable(df, options = list(pageLength = 10), rownames = TRUE)
  })


  output$prediction_table <- renderDT({
    req(input$model_choice, input$dataset_choice)
    if (input$dataset_choice == "Uploaded") {
      df <- predict_uploaded_data(uploaded_exprs(), input$model_choice)
    } else {
      ds <- ifelse(input$dataset_choice == "Combined (GSE10810 and GSE17907)", "Combined", input$dataset_choice)
      df <- all_results_df[[ds]][[input$model_choice]]
    }
    datatable(df, options = list(pageLength = 10))
  })
  
  
  

  
  
  observeEvent(input$dataset_choice, {
    if (page() == "Predictive Modelling") {
      isolate({
        updateTabsetPanel(session, "modelling_tab", selected = input$modelling_tab)
      })
    } else if (page() == "Evaluation") {
      isolate({
        updateTabsetPanel(session, "evaluation_tab", selected = input$evaluation_tab)
      })
    }
  })
  

  
  observeEvent(input$select_all, {
    req(input$dataset_choice, input$model_choice)
    df <- all_results_df[[input$dataset_choice]][[input$model_choice]]
    req(!is.null(df))
    
    proxy <- dataTableProxy("prediction_table")
    selectRows(proxy, 1:nrow(df))
  })
  observeEvent(input$clear_all, {
    proxy <- dataTableProxy("prediction_table")
    selectRows(proxy, NULL)
  })

  

  output$gene_box <- renderPlot({
    df <- load_data(input$dataset_choice)
    
    if (nrow(df) == 0 || !"Grade" %in% colnames(df)) return(NULL)
    
    expr_mat <- df[, sapply(df, is.numeric), drop = FALSE]
    
    gene_vars <- apply(expr_mat, 2, var, na.rm = TRUE)
    top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:input$top_gene_slider])
    
    expr_long <- reshape2::melt(df[, c("Grade", top_genes), drop = FALSE])
    
    ggplot(expr_long, aes(x = Grade, y = value, fill = Grade)) +
      geom_boxplot(alpha = 0.6) +
      facet_wrap(~variable, scales = "free") +
      theme_minimal() +
      labs(x = "Grade", y = "Expression", title = "Gene Expression by Grade")
  })
  
  
  
  
  
  
  output$bar_plot <- renderPlot({
    df <- all_results_df[[input$dataset_choice]][[input$model_choice]]
    req(df)

    if (nrow(df) == 0) return(NULL)
    
    
    req("TrueGrade" %in% colnames(df), "PredictedGrade" %in% colnames(df))
    
    if (!isTRUE(input$use_all_data)) {
      selected <- input$prediction_table_rows_selected
      if (length(selected) > 0) {
        df <- df[selected, ]
      } else {
        return(NULL)  # Avoid empty plot
      }
    }
    
    df$TrueGrade <- factor(df$TrueGrade, levels = c("0", "1", "2", "3"))
    df$PredictedGrade <- factor(df$PredictedGrade, levels = c("0", "1", "2", "3"))
    
    if (input$bar_plot_style == "facet") {
      # ✅ Bar + Line Combo
      library(dplyr)
      
      pred_count <- df %>%
        count(PredictedGrade) %>%
        rename(Grade = PredictedGrade, PredictedCount = n)
      
      true_count <- df %>%
        count(TrueGrade) %>%
        rename(Grade = TrueGrade, TrueCount = n)
      
      plot_df <- full_join(pred_count, true_count, by = "Grade") %>%
        mutate(across(c(PredictedCount, TrueCount), ~replace_na(., 0))) %>%
        mutate(Grade = factor(Grade, levels = c("0", "1", "2", "3")))
      
      ggplot(plot_df, aes(x = Grade)) +
        geom_col(aes(y = PredictedCount, fill = "Predicted"), width = 0.6) +
        geom_line(aes(y = TrueCount, group = 1, color = "True"), size = 1.5) +
        geom_point(aes(y = TrueCount, color = "True"), size = 3) +
        scale_fill_manual(values = c("Predicted" = "#e3a5b2")) +
        scale_color_manual(values = c("True" = "#2d2d2d")) +
        theme_minimal(base_size = 13) +
        labs(
          title = "True vs Predicted Grade Distribution",
          y = "Number of Samples",
          x = "Grade",
          fill = "",
          color = ""
        )
      
    } else {
      # ✅ Confusion Matrix Heatmap
      library(dplyr)
      conf_mat <- df %>%
        count(TrueGrade, PredictedGrade)
      
      ggplot(conf_mat, aes(x = TrueGrade, y = PredictedGrade, fill = n)) +
        geom_tile(color = "white") +
        geom_text(aes(label = n), size = 5) +
        scale_fill_gradient(low = "#fcf8f3", high = "#e3a5b2") +
        theme_minimal(base_size = 13) +
        labs(
          title = "Confusion Matrix: True vs Predicted Grade",
          x = "True Grade",
          y = "Predicted Grade",
          fill = "Count"
        )
    }
  })
  
  
  
  
  output$heatmap_plot <- renderPlot({
    df <- load_data(input$dataset_choice)
    if (nrow(df) == 0 || !"Grade" %in% colnames(df)) return(NULL)
    
    expr_mat <- df[, sapply(df, is.numeric), drop = FALSE]
    
    gene_vars <- apply(expr_mat, 2, var, na.rm = TRUE)
    top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:10])
    
    if (!isTRUE(input$use_all_data)) {
      selected <- input$prediction_table_rows_selected
      if (length(selected) > 0) {
        expr <- as.matrix(expr_mat[selected, top_genes, drop = FALSE])
        rownames(expr) <- df$Patient[selected]
      } else {
        return(NULL)
      }
    } else {
      expr <- as.matrix(expr_mat[, top_genes, drop = FALSE])
      rownames(expr) <- df$Patient
    }
    
    pheatmap(
      t(expr),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 8,
      main = "Top 10 Variable Genes",
      color = colorRampPalette(c("#e3a5b2", "#fcf8f3", "#c8b0d8"))(50)
    )
  })
  
  
  
  output$overview_table <- DT::renderDataTable({
    req(overview_df())
    datatable(overview_df(), options = list(pageLength = 10), rownames = TRUE)
  })
  
  

  
  
  output$fscore_table <- DT::renderDataTable({
    dataset <- input$dataset_choice
    model <- input$model_choice
    
    model_obj <- all_model_list[[dataset]][[model]]
    req(model_obj)
    
    pred <- predict(model_obj)
    true <- model_obj$trainingData[, ncol(model_obj$trainingData)]
    cm <- caret::confusionMatrix(pred, true)
    
    stats <- as.data.frame(cm$byClass)
    stats <- round(stats[, c("Precision", "Recall", "F1")], 3)
    stats$Grade <- rownames(stats)
    
    stats <- stats[, c("Grade", "Precision", "Recall", "F1")]
    
    stats$Grade <- gsub("Class: ", "", stats$Grade)
    
    DT::datatable(stats,
                  options = list(
                    paging = FALSE,
                    searching = FALSE,
                    info = FALSE,
                    ordering = FALSE,
                    autoWidth = TRUE
                  ),
                  rownames = FALSE
    )
  })
  
  
  
  
  
  output$roc_plot <- renderPlot({
    library(pROC)
    dataset <- input$dataset_choice
    model <- input$model_choice
    
    model_obj <- all_model_list[[dataset]][[model]]
    req(model_obj)
    
    pred_probs <- model_obj$pred
    obs <- pred_probs$obs
    classes <- levels(obs)
    
    colors <- c("#c8b0d8", "#e3a5b2", "#a2b5c4", "#b1d1c5")  # Morandi tones
    auc_list <- c()
    
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), 
         xlab = "False Positive Rate", ylab = "True Positive Rate",
         main = paste("ROC Curves -", model, "on", dataset),
         cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
    abline(a = 0, b = 1, lty = 2, col = "gray")
    
    for (i in seq_along(classes)) {
      class <- classes[i]
      bin_obs <- ifelse(obs == class, 1, 0)
      roc_obj <- roc(response = bin_obs, predictor = pred_probs[[class]], quiet = TRUE)
      auc_val <- round(auc(roc_obj), 3)
      auc_list[class] <- auc_val
      
      lines(roc_obj, col = colors[i], lwd = 2)
    }
    
    legend_labels <- paste0("Grade ", names(auc_list), " (AUC = ", auc_list, ")")
    
    legend("bottomleft", legend = legend_labels,
           col = colors[seq_along(auc_list)], lwd = 2, bty = "n", cex = 1)
    
  }, width = 600, height = 400)
  
  
  output$comparison_plot <- renderPlot({
    comp_data <- comparison_results()
    req(comp_data)
    
    top_choice <- input$comp_topn
    cv_choice <- input$comp_cv
    
    if (is.null(comp_data[[top_choice]]) || is.null(comp_data[[top_choice]][[cv_choice]])) {
      plot.new()
      text(0.5, 0.5, "No data available for this selection.", cex = 1.3)
      return()
    }
    
    model_names <- c("SVM", "KNN", "RF")
    df_list <- lapply(model_names, function(model) {
      acc <- comp_data[[top_choice]][[cv_choice]][[model]]
      if (!is.null(acc)) {
        data.frame(Model = model, Accuracy = acc)
      } else {
        NULL
      }
    })
    
    df <- do.call(rbind, df_list)
    req(nrow(df) > 0)
    
    ggplot(df, aes(x = Model, y = Accuracy, fill = Model)) +
      geom_boxplot(width = 0.6, alpha = 0.8) +
      theme_minimal(base_size = 13) +
      labs(title = paste("Model Accuracy Comparison (", top_choice, ",", cv_choice, ")"),
           x = "Model", y = "Accuracy") +
      scale_fill_manual(values = c("SVM" = "#c8b0d8", "KNN" = "#e3a5b2", "RF" = "#b1d1c5")) +
      theme(legend.position = "none")
  })
  
  
  
  output$mean_accuracy_text <- renderPrint({
    comp_data <- comparison_results()
    req(comp_data)
    
    top_choice <- input$comp_topn
    cv_choice <- input$comp_cv
    
    if (is.null(comp_data[[top_choice]]) ||
        is.null(comp_data[[top_choice]][[cv_choice]])) return("No data available.")
    
    for (model in c("SVM", "KNN", "RF")) {
      acc <- comp_data[[top_choice]][[cv_choice]][[model]]
      if (!is.null(acc)) {
        cat(paste0(model, ": Mean Accuracy = ", round(mean(acc), 3), "\n"))
      }
    }
  })
  
  
  
  output$summary_plot <- renderPlot({
    req(summary_df)
    
    eval_long <- tidyr::pivot_longer(summary_df, 
                                     cols = c("Accuracy", "MacroF1", "Kappa"),
                                     names_to = "Metric", values_to = "Value")
    
    ggplot(eval_long, aes(x = GeneSet, y = Value, fill = Model)) +
      geom_col(position = position_dodge(0.8)) +
      facet_wrap(~ Metric, scales = "free_y") +
      labs(title = "Test Set Evaluation: Accuracy, Macro F1, Kappa",
           x = "Top Gene Set", y = "Metric") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })
  
  
  
  
  
  
}  

shinyApp(ui, server)