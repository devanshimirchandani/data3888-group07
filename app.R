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
library(smotefamily)
library(pROC)

load_data <- function(dataset = "GSE15852") {
  df <- data.frame(Patient = character(), Grade = factor(levels = c("0", "1", "2", "3")))
  
  
  if (dataset == "GSE15852") {
    print("Loading GSE15852.RData")
    load("GSE15852.RData")  
    df <- as.data.frame(t(exprs_15852))
    df$Patient <- rownames(df)
    if (!is.null(metadata_15852$grade)) {
      df$Grade <- metadata_15852$grade[match(df$Patient, rownames(metadata_15852))]
    }
    if ("histopathological exam:ch1" %in% colnames(metadata_15852)) {
      df$Histopathology <- metadata_15852$`histopathological exam:ch1`[match(df$Patient, rownames(metadata_15852))]
    }
    
    return(df)
  }
  else if (dataset == "GSE10810") {
    load("data_GSE10810.RData")
    df <- as.data.frame(t(expr_10810))
    df$Patient <- rownames(df)
    if (exists("metadata_10810") && "grade" %in% colnames(metadata_10810)) {
      df$Grade <- metadata_10810$grade[match(df$Patient, rownames(metadata_10810))]
    }
    
  } else if (dataset == "GSE17907") {
    load("data_GSE17907.RData")
    df <- as.data.frame(t(expr_17907))
    df$Patient <- rownames(df)
    if (exists("metadata_17907") && "grade" %in% colnames(metadata_17907)) {
      df$Grade <- metadata_17907$grade[match(df$Patient, rownames(metadata_17907))]
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
    Patient         = rownames(prob),
    TrueGrade       = gsub("X", "", as.character(labels_test)),
    PredictedGrade  = gsub("X", "", as.character(pred)),
    X0.Chance       = prob[, "X0"],
    X1.Chance       = prob[, "X1"],
    X2.Chance       = prob[, "X2"],
    X3.Chance       = prob[, "X3"],
    check.names = FALSE,    # 保持列名不被自动改成合法 R 名称
    stringsAsFactors = FALSE
  )
  
  
  df$`Normal Chance`   <- paste0(round(df$`X0.Chance` * 100, 1), "%")
  df$`Grade 1 Chance`  <- paste0(round(df$`X1.Chance` * 100, 1), "%")
  df$`Grade 2 Chance`  <- paste0(round(df$`X2.Chance` * 100, 1), "%")
  df$`Grade 3 Chance`  <- paste0(round(df$`X3.Chance` * 100, 1), "%")
  
  
  df[, c("X0.Chance","X1.Chance","X2.Chance","X3.Chance")] <- NULL
  
  
  colnames(df) <- c(
    "Patient",
    "True Grade",
    "Predicted Grade",
    "Normal Chance",
    "Grade 1 Chance",
    "Grade 2 Chance",
    "Grade 3 Chance"
  )
  
  return(df)
}


prepare_predictions_all_models <- function() {
  load("GSE15852.RData")
  M    <- exprs_15852
  meta <- metadata_15852
  
  # 2. clean grade into "0","1","2","3"
  gv <- tolower(meta$grade)
  gv <- gsub("grade ?", "", gv)
  gv[gv=="normal"] <- "0"
  gv <- trimws(gv)
  fac_all <- factor(paste0("X", gv), levels = c("X0","X1","X2","X3"))
  
  # 3. train/test split
  set.seed(123)
  idx <- createDataPartition(fac_all, p=0.8, list=FALSE)
  expr_tr <- M[, idx]; expr_te <- M[, -idx]
  y_tr    <- fac_all[idx];   y_te    <- fac_all[-idx]
  
  # 4. SMOTE on training set
  df_tr <- as.data.frame(t(expr_tr))
  df_tr$grade <- y_tr
  sm <- SMOTE(X=df_tr[,-ncol(df_tr)], target=df_tr$grade,
              K=5, dup_size=5)$data
  expr_tr_bal <- t(as.matrix(sm[,-ncol(sm)]))
  colnames(expr_tr_bal) <- paste0("SMOTE_", seq_len(nrow(sm)))
  y_tr_bal <- factor(sm$class, levels=c("X0","X1","X2","X3"))
  
  # 5. limma feature selection
  design <- model.matrix(~0 + y_tr_bal); colnames(design) <- levels(y_tr_bal)
  fit    <- lmFit(expr_tr_bal, design) |> eBayes()
  res    <- topTable(fit, number=Inf, adjust.method="BH")
  sig1   <- rownames(res)[res$adj.P.Val < .05 ]
  gm     <- sapply(levels(y_tr_bal), function(l) rowMeans(expr_tr_bal[,y_tr_bal==l]))
  diff1  <- apply(gm[sig1,,drop=FALSE],1,function(x)max(x)-min(x))
  sig2   <- sig1[diff1 > 1]
  r2     <- res[sig2, ]; r2 <- r2[order(r2$adj.P.Val),]
  
  # 6. build top gene sets
  genesets <- list()
  if(nrow(r2)>=30)  genesets$top30  <- rownames(r2)[1:30]
  if(nrow(r2)>=50)  genesets$top50  <- rownames(r2)[1:50]
  if(nrow(r2)>=80)  genesets$top80  <- rownames(r2)[1:80]
  if(nrow(r2)>=110) genesets$top110 <- rownames(r2)[1:110]
  
  # 7. train & predict for each set
  results <- list(); models <- list()
  for(gs in names(genesets)) {
    g    <- genesets[[gs]]
    xtr  <- t(expr_tr_bal[g, ,drop=FALSE])
    xte  <- t(expr_te   [g, ,drop=FALSE]); rownames(xte)<-colnames(expr_te)
    ytrf <- y_tr_bal;          yte  <- y_te
    
    ctrl <- trainControl(
      method="cv", number=5,
      classProbs=TRUE, summaryFunction=multiClassSummary,
      savePredictions="final"
    )
    svm_mod <- train(x=xtr,y=ytrf,method="svmLinear",trControl=ctrl)
    knn_mod <- train(x=xtr,y=ytrf,method="knn",      trControl=ctrl,
                     tuneGrid=expand.grid(k=5:10))
    rf_mod  <- train(x=xtr,y=ytrf,method="rf",       trControl=ctrl)
    
    ps <- list(
      SVM = list(pred=predict(svm_mod,xte), prob=predict(svm_mod,xte,type="prob")),
      KNN = list(pred=predict(knn_mod,xte), prob=predict(knn_mod,xte,type="prob")),
      RF  = list(pred=predict(rf_mod, xte), prob=predict(rf_mod, xte,type="prob"))
    )
    
    results[[gs]] <- lapply(ps, function(z) make_result(z$pred, z$prob, yte, xte))
    models[[gs]]  <- list(SVM=svm_mod, KNN=knn_mod, RF=rf_mod)
  }
  
  list(result=results, models=models, genesets=genesets)
}

tmp <- prepare_predictions_all_models()
all_results_df <- list(GSE15852 = tmp$result)
all_model_list <- list(GSE15852 = tmp$models)
all_genesets_df <- list(GSE15852 = tmp$genesets)




prepare_predictions_for <- function(dataset_name) {
  if (dataset_name == "GSE17907") {
    load("data_GSE17907.RData")
    M <- expr_17907
    meta <- metadata_17907
  } else if (dataset_name == "Combined") {
    load("data_Combined.RData")
    M <- expr_corrected
    meta <- metadata
  } else {
    stop("Unsupported dataset")
  }
  # 2. clean grade into "0","1","2","3"
  gv <- tolower(meta$grade)
  gv <- gsub("grade ?", "", gv)
  gv[gv=="normal"] <- "0"
  gv <- trimws(gv)
  fac_all <- factor(paste0("X", gv), levels = c("X0","X1","X2","X3"))
  
  # 3. train/test split
  set.seed(123)
  idx <- createDataPartition(fac_all, p=0.8, list=FALSE)
  expr_tr <- M[, idx]; expr_te <- M[, -idx]
  y_tr    <- fac_all[idx];   y_te    <- fac_all[-idx]
  
  # 4. SMOTE on training set
  df_tr <- as.data.frame(t(expr_tr))
  df_tr$grade <- y_tr
  sm <- SMOTE(X=df_tr[,-ncol(df_tr)], target=df_tr$grade,
              K=1, dup_size=5)$data
  expr_tr_bal <- t(as.matrix(sm[,-ncol(sm)]))
  colnames(expr_tr_bal) <- paste0("SMOTE_", seq_len(nrow(sm)))
  y_tr_bal <- factor(sm$class, levels=c("X0","X1","X2","X3"))
  
  # 5. limma feature selection
  design <- model.matrix(~0 + y_tr_bal); colnames(design) <- levels(y_tr_bal)
  fit    <- lmFit(expr_tr_bal, design) |> eBayes()
  res    <- topTable(fit, number=Inf, adjust.method="BH")
  sig1   <- rownames(res)[res$adj.P.Val < .05 ]
  gm     <- sapply(levels(y_tr_bal), function(l) rowMeans(expr_tr_bal[,y_tr_bal==l]))
  diff1  <- apply(gm[sig1,,drop=FALSE],1,function(x)max(x)-min(x))
  sig2   <- sig1[diff1 > 1]
  r2     <- res[sig2, ]; r2 <- r2[order(r2$adj.P.Val),]
  
  # 6. build top gene sets
  genesets <- list()
  if(nrow(r2)>=30)  genesets$top30  <- rownames(r2)[1:30]
  if(nrow(r2)>=50)  genesets$top50  <- rownames(r2)[1:50]
  if(nrow(r2)>=80)  genesets$top80  <- rownames(r2)[1:80]
  if(nrow(r2)>=110) genesets$top110 <- rownames(r2)[1:110]
  
  # 7. train & predict for each set
  results <- list(); models <- list()
  for(gs in names(genesets)) {
    g    <- genesets[[gs]]
    xtr  <- t(expr_tr_bal[g, ,drop=FALSE])
    xte  <- t(expr_te   [g, ,drop=FALSE]); rownames(xte)<-colnames(expr_te)
    ytrf <- y_tr_bal;          yte  <- y_te
    
    ctrl <- trainControl(
      method="cv", number=5,
      classProbs=TRUE, summaryFunction=multiClassSummary,
      savePredictions="final"
    )
    svm_mod <- train(x=xtr,y=ytrf,method="svmLinear",trControl=ctrl)
    knn_mod <- train(x=xtr,y=ytrf,method="knn",      trControl=ctrl,
                     tuneGrid=expand.grid(k=5:10))
    rf_mod  <- train(x=xtr,y=ytrf,method="rf",       trControl=ctrl)
    
    ps <- list(
      SVM = list(pred=predict(svm_mod,xte), prob=predict(svm_mod,xte,type="prob")),
      KNN = list(pred=predict(knn_mod,xte), prob=predict(knn_mod,xte,type="prob")),
      RF  = list(pred=predict(rf_mod, xte), prob=predict(rf_mod, xte,type="prob"))
    )
    
    results[[gs]] <- lapply(ps, function(z) make_result(z$pred, z$prob, yte, xte))
    models[[gs]]  <- list(SVM=svm_mod, KNN=knn_mod, RF=rf_mod)
  }
  
  list(result=results, models=models, genesets=genesets)
}



for (ds in c("GSE17907", "Combined")) {
  tmp <- prepare_predictions_for(ds)
  all_results_df[[ds]] <- tmp$result
  all_model_list[[ds]] <- tmp$models
  all_genesets_df[[ds]] <- tmp$genesets
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

# Combined Comparison
comparison_results_combined <- reactiveVal(NULL)


ui <- fluidPage(
  introjsUI(), 
  useShinyjs(),
  theme = bs_theme(bg = "#fcf8f3", fg = "#2d2d2d", primary = "#e3a5b2"),
  tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Poppins&display=swap"),
  tags$style(HTML("
    body {
      font-family: Poppins, sans-serif;
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
  

    h1, h2, h3 {
      font-family: Poppins, sans-serif;
      font-weight: 600;
    }

    .nav-tab, .btn, .action-button {
      font-family: Poppins, sans-serif;
      letter-spacing: 0.4px;
    }

  h1 {
    font-size: 24px;
    font-weight: 700;
    color: #d07095;  
    margin-top: 20px !important;
    margin-bottom: 15px;
  }
  
  h2 {
    font-size: 21px;
    font-weight: 700;
    color: #d58aa1;  
    margin-top: 20px !important;
    margin-bottom: 12px;
  }
  
  h3 {
    font-size: 18px;
    font-weight: 600;
    color: #785d84; 
    margin-top: 15px !important;
    margin-bottom: 8px;
  }

  .browse-instruction {
  background-color: #fdf6f9;
  border: 1px solid #f3d4dc;
  padding: 12px 18px;
  border-radius: 8px;
  margin-bottom: 12px;
  font-size: 14px;
  color: #4e3a5c;
  line-height: 1.5;
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
    
    position: relative;
    
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
    border-radius: 999px !important;
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
    margin: 12px 20px 6px 20px; 
    padding: 10px 16px;
    background-color: #fcf8f3;  
    border: 1px solid #f3d4dc; 
    border-radius: 12px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  }
  
  .selector-row label {
  font-weight: 500;
  color: #785d84;
  } 
  
  .selectize-control.single .selectize-input {
  font-family: Poppins, sans-serif;
  font-size: 14px;
  padding: 6px 12px;
  border: 1px solid #eecedf;
  border-radius: 6px;
  background-color: white;
  box-shadow: none;
  }
  
  table.dataTable.stripe tbody tr.odd {
    background-color: #fdf6f9;
  }
  
  .dataTables_filter input {
    border: 1px solid #eecedf;
    border-radius: 6px;
    padding: 4px 8px;
    font-family: Poppins, sans-serif;
  }
  
  .dataTable, .dataTables_wrapper {
    font-size: 15px;
    font-family: Poppins, sans-serif;
    color: #4e3a5c;
  }
  
  .dataTable th {
    font-size: 14px;
    font-weight: 600;
    color: #785d84;
  }

  table.dataTable tbody td {
  padding-top: 8px;
  padding-bottom: 8px;
  }

  .dataset-card {
    background-color: #fcf8f3;
    border: 1px solid #f3d4dc;
    padding: 16px 20px;
    border-radius: 12px;
    margin-bottom: 16px;
    height: 100%;
  }
  
  .dataset-card h3 {
    margin-top: 0;
    color: #a85888;
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
  

  
  .lo-tab .nav-tabs > li > a {
    color: #e3a5b2 !important;
    font-weight: 500;
    padding: 6px 16px;
    margin: 6px 6px 0 0;
    font-size: 18px;
    border-radius: 999px !important;
    border: none !important;
    background-color: transparent !important;
    transition: all 0.3s ease;
    box-shadow: none !important;
    
  }
  
  .lo-tab .nav-tabs > li > a:hover {
    background-color: #fce9ef !important;
    border-color: #e3a5b2 !important;
    color: #d07095 !important;
    box-shadow: 0 1px 2px rgba(227,165,178,0.25) !important;
  }
  
  .lo-tab .nav-tabs > li.active > a,
  .lo-tab .nav-tabs > li.active > a:focus,
  .lo-tab .nav-tabs > li.active > a:hover {
    background-color: #e3a5b2 !important;
    color: white !important;
    border: 1px solid #e3a5b2 !important;
    font-weight: bold;
    box-shadow: inset 0 -2px 0 #c66690;
    border-radius: 999px !important;
  }
  
  .lo-tab .nav-tabs > li > a[aria-selected='true'] {
    background-color: #e3a5b2 !important;
    color: white !important;
    font-weight: bold;
    border-radius: 999px !important;
    box-shadow: inset 0 -1px 0 #c66690 !important;
    border: 0.5px solid #e3a5b2 !important;
  }



  .sub-tab .nav-tabs > li > a {
    color: #785d84 !important;
    font-weight: 500;
    padding: 5px 14px;
    margin: 6px 6px 0 0;
    font-size: 15px; 
    border-radius: 999px !important;
    border: none !important;
    background-color: transparent !important;
    transition: all 0.3s ease;
    box-shadow: none !important;
  }
  
  .sub-tab .nav-tabs > li > a:hover {
    background-color: #d8c3e5 !important; 
    color: #4e3a5c !important;
    box-shadow: 0 1px 2px rgba(120, 93, 132, 0.3) !important;
  }
  
  .sub-tab .nav-tabs > li.active > a,
  .sub-tab .nav-tabs > li.active > a:focus,
  .sub-tab .nav-tabs > li.active > a:hover {
    background-color: #c8b0d8 !important; 
    color: white !important;
    font-weight: bold;
    border: 1px solid #785d84 !important;
    box-shadow: inset 0 -2px 0 #5e456a;
    border-radius: 999px !important;
  }
  

  .sub-tab .nav-tabs > li > a[aria-selected='true'] {
    background-color: #c8b0d8 !important;
    color: white !important;
    font-weight: bold;
    border-radius: 999px !important;
    box-shadow: inset 0 -1px 0 #5e456a !important;
    border: 0.5px solid #785d84 !important;
  }



")),
  
  div(class = "topbar-container",
      div(
        style = "flex: 1; display: flex; align-items: center; justify-content: flex-start;",
        tags$a(
          href = "https://github.com/devanshimirchandani/data3888-group07",
          target = "_blank",
          rel = "noopener noreferrer",
          title = "Click to view Preva introduction and code on Bilibili",  # tooltip
          img(
            src = "Preva.png",
            class = "logo-title",
            `data-step` = 1,
            `data-intro` = "Click the Preva logo to visit our GitHub or video page for more details, source code, and team introduction."
          )
        )
      )
      ,
      div(class = "toolbar-group", style = "flex: 5; display: flex; justify-content: center; gap: 24px;",
          actionButton("about_tab", "About Preva", class = "nav-tab",
                       `data-step` = 2,
                       `data-intro` = "Click here to learn what Preva is about."),
          actionButton("intro_tab", "Module 1", class = "nav-tab",
                       `data-step` = 3,
                       `data-intro` = "This section introduces you to breast cancer and gene data."),
          actionButton("data_tab", "Module 2", class = "nav-tab",
                       `data-step` = 4,
                       `data-intro` = "Select this tab to upload, browse, and prepare your datasets."),
          actionButton("analysis_tab", "Module 3", class = "nav-tab",
                       `data-step` = 5,
                       `data-intro` = "This is where you build machine learning models."),
          actionButton("eval_tab", "Module 4", class = "nav-tab",
                       `data-step` = 6,
                       `data-intro` = "Finally, check how well your model performed!")
      ),
      div(style = "flex: 1; display: flex; justify-content: flex-end; align-items: center;",
          actionButton("help_modal_button", label = NULL, icon = icon("question-circle"),
                       style = "background: none; border: none; color: #e3a5b2; font-size: 22px;",
                       class = "nav-tab",
                       `data-step` = 7,
                       `data-intro` = "Click here anytime to open the Help panel and revisit this tour.")
          
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
  
  quiz_submitted_m1 <- reactiveVal(FALSE)
  
  quiz_submitted_m2 <- reactiveVal(FALSE)
  
  quiz_submitted_m3 <- reactiveVal(FALSE)
  
  quiz_submitted_m4 <- reactiveVal(FALSE)
  
  
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
  
  
  observeEvent(input$module1_tabs, {
    if (input$module1_tabs == "lo2_gene") {
      updateTabsetPanel(
        session,
        "module1_lo2_subtabs",
        selected = "inner_gene"
      )
    } else if (input$module1_tabs == "check_quiz") {
      updateTabsetPanel(
        session,
        "module1_quiz_subtabs",
        selected = "quiz"
      )
    }
  })
  
  
  observeEvent(input$module2_tabs, {
    if (input$module2_tabs == "LO3: Using GEO DB/GSE Data") {
      updateTabsetPanel(session, "module2_lo3_subtabs", selected = "About Datasets")
    } else if (input$module2_tabs == "LO4: Data Pre-processing") {
      updateTabsetPanel(session, "module2_lo4_subtabs", selected = "Pre-processing Steps")
    } else if (input$module2_tabs == "Check Your Understanding") {
      updateTabsetPanel(session, "module2_quiz_subtabs", selected = "Quiz")
    }
  })
  
  observeEvent(input$module3_tabs, {
    if (input$module3_tabs == "LO5: Predictive Modelling") {
      updateTabsetPanel(session, "module3_lo5_subtabs", selected = "Model Building")
    }
    else if (input$module3_tabs == "Check Your Understanding") {
      updateTabsetPanel(session, "module3_quiz_subtabs", selected = "Quiz")
    }
  })
  
  observeEvent(input$module4_tabs, {
    if (input$module4_tabs == "LO6: Model Evaluation") {
      updateTabsetPanel(session, "module4_lo6_subtabs", selected = "Glossary")
    } else if (input$module4_tabs == "LO7: Communicating Findings") {
      updateTabsetPanel(session, "module4_lo7_subtabs", selected = "Key Findings Summary")
    }else if (input$module4_tabs == "Check Your Understanding") {
      updateTabsetPanel(session, "module4_quiz_subtabs", selected = "Quiz")
    }
  })
  
  
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
    req(input$dataset_choice, input$model_choice, input$topn_choice)
    ds    <- input$dataset_choice
    model <- input$model_choice
    gs    <- paste0("top", input$topn_choice)  # 构造 geneset 名称: "top30","top50",...
    
    if (ds == "Uploaded" && !is.null(uploaded_exprs())) {
      df <- predict_uploaded_data(uploaded_exprs(), model)
    } else {
      # 确保这个 dataset、geneset、model 都存在
      if (
        !ds %in% names(all_results_df) ||
        !gs %in% names(all_results_df[[ds]]) ||
        !model %in% names(all_results_df[[ds]][[gs]])
      ) {
        updateCheckboxGroupInput(session, "selected_genes", choices = character(0))
        return()
      }
      # 从三维列表中取出对应的 dataframe
      df <- all_results_df[[ds]][[gs]][[model]]
    }
    
    # 更新可选基因列表（如果 prediction table 有基因列的话）
    if (!is.null(df) && ncol(df) >= 6) {
      gene_cols <- colnames(df)[6:ncol(df)]
      updateCheckboxGroupInput(session, "selected_genes", choices = gene_cols)
    } else {
      updateCheckboxGroupInput(session, "selected_genes", choices = character(0))
    }
  })
  
  
  
  page <- reactiveVal("About")
  
  observeEvent(input$about_tab, {
    page("About")
    runjs("
    $('.nav-tab').removeClass('active');
    $('#about_tab').addClass('active');
  ")
  })
  
  observeEvent(input$intro_tab, {
    page("Module 1")
    runjs("
    $('.nav-tab').removeClass('active');
    $('#intro_tab').addClass('active');
  ")
    updateTabsetPanel(session, "module1_tabs", selected = "lo1_tab")
    updateTabsetPanel(session, "module1_lo1_subtabs", selected = "Breast Cancer")
  })

  observeEvent(input$data_tab, {
    page("Module 2")
    runjs("
    $('.nav-tab').removeClass('active');
    $('#data_tab').addClass('active');
  ")
    updateTabsetPanel(session, "module2_tabs", selected = "LO3: Using GEO DB/GSE Data")
    updateTabsetPanel(session, "module2_lo3_subtabs", selected = "About Datasets")
  })
                       
  observeEvent(input$analysis_tab, {
    page("Module 3")
    runjs("
    $('.nav-tab').removeClass('active');
    $('#analysis_tab').addClass('active');
  ")
    updateTabsetPanel(session, "module3_tabs", selected = "LO5: Predictive Modelling")
    updateTabsetPanel(session, "module3_lo5_subtabs", selected = "Model Building")
  })
  
  observeEvent(input$eval_tab, {
    page("Module 4")
    runjs("
    $('.nav-tab').removeClass('active');
    $('#eval_tab').addClass('active');
  ")
    updateTabsetPanel(session, "module4_tabs", selected = "LO6: Model Evaluation")
    updateTabsetPanel(session, "module4_lo6_subtabs", selected = "Glossary")
  })

  
  
  observeEvent(input$data_tab, page("Module 2"))
  observeEvent(input$analysis_tab, page("Module 3"))
  observeEvent(input$eval_tab, page("Module 4"))
  
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
  
  
  
  # Introduction Tour
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
      <p><b>Preva</b> is an educational technology (EdTech) application designed to help users explore gene expression data and understand how it relates to breast cancer progression using machine learning.</p>
      <p>It follows a structured learning design made up of <b>Modules</b> and <b>Learning Objectives (LOs)</b>:</p>
      
      <ul>
        <li><b>Module 1: Introduction</b>
          <ul>
            <li><b>LO1:</b> Understand the nature and impact of breast cancer.</li>
            <li><b>LO2:</b> Explore differentially expressed genes using expression data.</li>
          </ul>
        </li>
        
        <li><b>Module 2: Data Processing</b>
          <ul>
            <li><b>LO3:</b> Navigate and use public GEO/GSE datasets.</li>
            <li><b>LO4:</b> Pre-process gene data (cleaning, transformation, annotation).</li>
          </ul>
        </li>
        
        <li><b>Module 3: Model Building</b>
          <ul>
            <li><b>LO5:</b> Build predictive models (SVM, KNN, RF) to identify key genes by cancer stage.</li>
          </ul>
        </li>
        
        <li><b>Module 4: Model Evaluation</b>
          <ul>
            <li><b>LO6:</b> Evaluate models using F1-score, ROC-AUC, Confusion Matrix, and Cross-validation.</li>
            <li><b>LO7:</b> Summarise findings and reflect on limitations and future directions.</li>
          </ul>
        </li>
    ")
    ))
  })
  
  observe({
    if (page() == "About") {
      runjs("
      $('.nav-tab').removeClass('active');
      $('#about_tab').addClass('active');
    ")
    }
  })
  
  output$conditional_selectors <- renderUI({
    if (page() == "Module 3") {
      div(class = "selector-row",
          fluidRow(
            column(4,
                   selectInput("dataset_choice", "Dataset:",
                               choices = c("GSE15852", "GSE17907", "Combined", "Uploaded"),
                               selected = "GSE15852")
            ),
            column(4,
                   selectInput("model_choice", "Model:",
                               choices = c("SVM", "KNN", "RF"),
                               selected = "SVM")
            ),
            column(4,
                   selectInput("topn_choice", "Top genes:",
                               choices = c("30", "50", "80", "110"),
                               selected = "30")
            )
          )
      )
    } else if (page() == "Module 4") {
      div(class = "selector-row",
          fluidRow(
            column(4,
                   selectInput("dataset_choice", "Dataset:",
                               choices = c("GSE15852", "GSE17907", "Combined", "Uploaded"),
                               selected = "GSE15852")
            ),
            column(4,
                   selectInput("model_choice", "Model:",
                               choices = c("SVM", "KNN", "RF"),
                               selected = "SVM")
            ),
            column(4,
                   selectInput("topn_choice", "Top genes:",
                               choices = c("30", "50", "80", "110"),
                               selected = "30")
            )
          )
      )
    }
  })
  
  output$main_ui <- renderUI({
    if (page() %in% c("Predictive Modelling", "Evaluation")) {
      req(input$dataset_choice)
    }
    
    
    switch(page(),
           "About" = {
             tagList(
               h1("About Preva"),
               p("Breast cancer is the most common cancer among females, and one in seven women are diagnosed with it in their lifetime. This is why we created Preva. This is a tool that medical researchers can use to enhance knowledge that facilitates early breast cancer detection, and understanding the interaction between genes and breast cancer diagnoses. But also to detect it when it matters, and know what stage the cancer is at for further action."),
               p("This Shiny web application uses machine learning models to help predict and visualise breast cancer diagnoses. The application, which was developed using R and is driven by interactive data visualisation, gives users an easy-to-use way to examine patient data and prediction outcomes."),
               tags$b("Sources:"),
               tags$ul(tags$li(a("National Breast Cancer Foundation (2021)", href = "https://nbcf.org.au/about-breast-cancer/breast-cancer-stats/", target = "_blank")))
             )
           },
           
           "Module 1" = {
             div(class = "nav-tabs lo-tab",
                 tabsetPanel(id = "module1_tabs", type = "tabs",
                             
                             tabPanel("LO1: Intro to Breast Cancer", value = "lo1_tab",
                                      div(class = "nav-tabs sub-tab",
                                          tabsetPanel(id = "module1_lo1_subtabs", type = "tabs",
                                                      tabPanel("Breast Cancer", value = "Breast Cancer",
                                                               h1("Why Focus on Breast Cancer?"),
                                                               
                                                               p("Breast cancer is the most common cancer among women worldwide, with over 2.3 million new cases diagnosed in 2020. Most breast cancers originate in the ducts (ductal carcinoma) or lobules (lobular carcinoma) of the breast. Risk factors include age, family history, obesity, alcohol use, and hormonal factors."),
                                                               
                                                               p("Early detection through screening has significantly improved survival rates. Treatment may involve surgery, chemotherapy, radiation, hormone therapy, or targeted drugs, depending on the subtype and stage."),
                                                               
                                                               h3("References"),
                                                               tags$ul(
                                                                 tags$li("Loibl, S. et al. 2021. “Breast Cancer.” The Lancet 397(10286): 1750–69. ", a("[Link]", href = "https://doi.org/10.1016/S0140-6736(21)00574-1", target = "_blank")),
                                                                 tags$li("American Cancer Society. 2023. “What Is Breast Cancer?” ", a("[Link]", href = "https://www.cancer.org/cancer/types/breast-cancer/about/what-is-breast-cancer.html", target = "_blank")),
                                                                 tags$li("Akram, M. et al. 2017. “Awareness and Current Knowledge of Breast Cancer.” Biological Research 50(1): 33. ", a("[Link]", href = "https://doi.org/10.1186/s40659-017-0138-8", target = "_blank"))
                                                               )
                                                      )
                                          )
                                      )
                             ),
                             
                             tabPanel("LO2: Gene Expression",  value = "lo2_gene", 
                                      div(class = "nav-tabs sub-tab",
                                          tabsetPanel(id = "module1_lo2_subtabs", type = "tabs",
                                                        tabPanel("Gene Expression", value = "inner_gene",
                                                                 h1("Understanding Breast Cancer Gene Expression"),
                                                                 
                                                                 p("Gene expression is the process by which information from DNA is transcribed into RNA and translated into proteins that regulate cellular function. Abnormal expression patterns—due to mutations or regulatory defects—can lead to tumour development and progression."),
                                                                 
                                                                 p("In breast cancer, gene expression profiling enables the identification of molecular subtypes and supports personalised diagnosis and treatment strategies. Commonly implicated genes include TP53, BRCA1, and BRCA2, which are linked to hereditary breast cancer risk."),
                                                                 
                                                                 p("In this app, we utilise gene expression data from the GSE15852 dataset, which includes tumour samples spanning Grade I to III. These histological grades reflect how aggressive or differentiated the cancer cells are. This dataset supports the classification of tumour subtypes based on expression patterns."),
                                                                 
                                                                 p("By analysing these patterns, Preva helps visualise expression variation across grades and link gene activity to clinical outcomes."),
                                                                 
                                                                 h3("References"),
                                                                 tags$ul(
                                                                   tags$li("Alberts, B. et al. 2015. Molecular Biology of the Cell. 6th ed. Garland Science."),
                                                                   tags$li("Yarden, Y. and Sliwkowski, M.X. 2001. “Untangling the ErbB Signalling Network.” Nature Reviews Molecular Cell Biology 2(2): 127–37."),
                                                                   tags$li("Reis-Filho, J.S. and Pusztai, L. 2011. “Gene Expression Profiling in Breast Cancer.” The Lancet 378(9805): 1812–23. ", a("[Link]", href = "https://doi.org/10.1016/S0140-6736(11)61539-0", target = "_blank")),
                                                                   tags$li("Wendt, C. and Margolin, S. 2019. “Identifying Breast Cancer Susceptibility Genes.” Acta Oncologica 58(2): 135–46. ", a("[Link]", href = "https://doi.org/10.1080/0284186X.2018.1529428", target = "_blank")),
                                                                   tags$li("Lindblom, A. 1995. “Familial Breast Cancer and Genes Involved in Breast Carcinogenesis.” Breast Cancer Research and Treatment 34(1): 171–83. ", a("[Link]", href = "https://doi.org/10.1007/BF00666004", target = "_blank"))
                                                                 )
                                                        )
                                                      
                                          )
                                      )
                             ),
                             tabPanel("Check Your Understanding",
                                      div(class = "nav-tabs sub-tab",
                                          tabsetPanel(id = "module1_quiz_subtabs", type = "tabs",
                                                      tabPanel("Quiz", uiOutput("quiz_game_ui_m1")),
                                                      tabPanel("Review & Summary", uiOutput("quiz_summary_ui_m1"))
                                          )
                                      )
                             )
                             
                             
                 )
             )
           },
           
           "Module 2" = {
             div(class = "nav-tabs lo-tab",
                 tabsetPanel(id = "module2_tabs", type = "tabs",
                             tabPanel("LO3: Using GEO DB/GSE Data",
                                      div(class = "nav-tabs sub-tab",
                                          tabsetPanel(id = "module2_lo3_subtabs", type = "tabs",
                                                      tabPanel("About Datasets", 
                                                               h1("About the Datasets"),
                                                               p("This application supports several public gene expression datasets relevant to breast cancer research. Below is a summary of each dataset, including the source and brief description."),
                                                               
                                                               fluidRow(
                                                                 column(4,
                                                                        div(class = "dataset-card",
                                                                          h2("GSE15852"),
                                                                          p("This dataset contains gene expression profiles from 43 normal breast samples and 43 tumour samples spanning Grade I to III. It is frequently used to study molecular subtypes of breast cancer."),
                                                                          tags$a("View on NCBI GEO", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15852", target = "_blank"),
                                                                          tableOutput("tbl_GSE15852")
                                                                        )),
                                                                 column(4,
                                                                        div(class = "dataset-card",
                                                                          h2("GSE17907"),
                                                                          p("GSE17907 includes gene expression data from 55 breast cancer samples, enriched for high-grade cases. It is often used for studying aggressive tumour phenotypes."),
                                                                          tags$a("View on NCBI GEO", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17907", target = "_blank"),
                                                                          tableOutput("tbl_GSE17907")
                                                                        )),
                                                                 column(4,
                                                                          div(class = "dataset-card",
                                                                            h2("Combined (GSE10810 + GSE17907)"),
                                                                            p("This dataset combines GSE10810 and GSE17907, harmonised using batch correction techniques. It provides a broader spectrum of samples for better model generalisation."),
                                                                            tags$a("GSE10810", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10810", target = "_blank"), br(),
                                                                            tags$a("GSE17907", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17907", target = "_blank"),
                                                                            tableOutput("tbl_Combined")
                                                                          ))
                                                               )
                                                      ),
                                                
                                                      tabPanel("Browse Dataset", 
                                                               div(class = "browse-instruction",
                                                                   p(strong("Instructions:"), " This panel allows you to upload and preview your own dataset (.RData)."),
                                                                   tags$ol(
                                                                     tags$li("Upload your .RData file using the 'Browse' button."),
                                                                     tags$li("Select the expression matrix (genes × samples) and phenotype (metadata) object from the dropdowns."),
                                                                     tags$li("Select the patient ID column and grade column if available."),
                                                                     tags$li("Click 'Check Dataset Overview' to preview and confirm that your dataset is loaded correctly.")
                                                                   ),
                                                                   p("Once loaded, your dataset will become available for prediction and analysis throughout Preva.")
                                                               ),
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
                                                               div(class = "browse-instruction",
                                                                   h3("Preview the Selected Dataset"),
                                                                   p("This section displays a preview of the dataset after pre-processing. You can use it to verify the uploaded or selected dataset, including the corrected expression data and matched phenotype metadata."),
                                                                   p("Use the dropdown below to switch between datasets. The corresponding table will update automatically.")
                                                               ),
                                                               
                                                               selectInput("dataset_choice_data", "Select Dataset:",
                                                                           choices = c("GSE15852", "GSE10810", "GSE17907", "Combined (GSE10810 and GSE17907)"),
                                                                           selected = "GSE15852"
                                                               ),
                                                               
                                                               DTOutput("data_table")
                                                      )
                                                      
                                          )
                                      )
                             ),
                             tabPanel("LO4: Data Pre-processing",
                                      div(class = "nav-tabs sub-tab",
                                          tabsetPanel(
                                            id   = "module2_lo4_subtabs",
                                            type = "tabs",
                                            tabPanel("Pre-processing Steps",
                                                     h1("How Preva Processes Data"),
                                                     p("Before training predictive models, the input gene expression datasets are pre-processed to ensure consistency, reliability, and fairness. These steps help reduce noise, correct technical artifacts, and address class imbalance."),
                                                    tags$ul(
                                                      tags$li(tags$b("Remove missing values:"), " Exclude genes or samples with missing values (NA) to avoid bias or computation errors."),
                                                      tags$li(tags$b("Standardise grade labels:"), " Convert dataset-specific grade annotations into consistent labels (Normal, Grade 1–3)."),
                                                      tags$li(tags$b("Filter low-expression genes:"), " Remove genes with minimal expression (e.g., <1 CPM) to reduce noise and improve signal clarity."),
                                                      tags$li(tags$b("Balance classes (SMOTE):"), " Use Synthetic Minority Over-sampling Technique (SMOTE) to address imbalance in tumour grades."),
                                                      tags$li(tags$b("Batch correction (ComBat):"), " Adjust for technical variation between datasets using the ComBat algorithm.")
                                                    ),
                                                    
                                                    p("These standardised pre-processing steps are essential for ensuring that downstream machine learning models are accurate, interpretable, and generalisable across different datasets.")
                                                  ),
                                            tabPanel("PCA",
                                                     h1("PCA After Pre-processing"),
                                                     
                                                     p("This PCA visualisation is based on gene expression data after pre-processing, including filtering, batch correction (ComBat), and grade label harmonisation. It helps evaluate how well samples cluster by grade or dataset origin."),
                                                     
                                                     div(class = "selector-row",
                                                         fluidRow(
                                                           column(6,
                                                                  selectInput(
                                                                    "dataset_choice", "Dataset:",
                                                                    choices  = c("GSE15852", "GSE17907", "Combined"),
                                                                    selected = "GSE15852"
                                                                  )
                                                           ),
                                                           column(6,
                                                                  selectInput(
                                                                    "topn_choice", "Top genes:",
                                                                    choices  = c("30", "50", "80", "110"),
                                                                    selected = "30"
                                                                  )
                                                           )
                                                         )
                                                     ),
                                                     
                                                     p("You can adjust the dataset and number of top variable genes to explore how sample separation changes under different settings."),
                                                     
                                                     div(style = "text-align: center;",
                                                         plotOutput("pca_plot", width = "600px", height = "420px")
                                                     )
                                            )
                                          ))),   
                             
                             
                             tabPanel("Check Your Understanding",
                                      div(class = "nav-tabs sub-tab",
                                          tabsetPanel(id = "module2_quiz_subtabs", type = "tabs",
                                                      tabPanel("Quiz", uiOutput("quiz_game_ui_m2")),
                                                      tabPanel("Review & Summary", uiOutput("quiz_summary_ui_m2"))
                                          )
                                      )
                             )
                 )
             )
           },
           
           "Module 3" = {
             tags$div(class = "lo-tab",
                      tabsetPanel(id = "module3_tabs", type = "tabs",
                                  tabPanel("LO5: Predictive Modelling",
                                           tags$div(class = "sub-tab",
                                                    tabsetPanel(id = "module3_lo5_subtabs", type = "tabs",
                                                                tabPanel("Model Building",
                                                                         h1("Overview"),
                                                                         
                                                                         h2("Support Vector Machine (SVM)"),
                                                                         p("SVM is a binary classification model that separates samples of different categories by finding an optimal hyperplane in a high-dimensional space. It can handle problems by mapping to a higher dimension through the kernel function"),
                                                                         
                                                                         h2("K-Nearest Neighbors (KNN)"),
                                                                         p("KNN calculates the distances between the samples to be classified and all the samples in the training set, selects the K neighbors with the closest distances, and then takes the voting results of these K neighbors as the final classification."),
                                                                         
                                                                         h2("Random Forest (RF)"),
                                                                         p("RF obtains the final classification by constructing multiple decision trees and voting on their prediction results. It randomly selects a feature subset during training, thereby reducing the risk of overfitting and being insensitive to missing values and outliers."),
                                                                         
                                                                        
                                                                ),
                                                                
                                                                tabPanel("Prediction",
                                                                         h1("Prediction Table"),
                                                                         div(class = "browse-instruction",
                                                                             h3("Understanding the Prediction Table"),
                                                                             p("This table displays the predicted tumour grades for each sample based on the selected dataset, model, and number of top genes."),
                                                                             tags$ul(
                                                                               tags$li(tags$b("True Grade:"), " The actual grade label from the dataset."),
                                                                               tags$li(tags$b("Predicted Grade:"), " The predicted label assigned by the current model."),
                                                                               tags$li(tags$b("Probabilities:"), " Confidence scores for each grade category.")
                                                                             ),
                                                                             p("You can use the 'Select All' and 'Clear Selection' buttons below to explore prediction results for different subsets.")
                                                                         ),
                                                                         fluidRow(
                                                                           column(6, actionButton("select_all", "Select All")),
                                                                           column(6, actionButton("clear_all", "Clear Selection"))
                                                                         ),
                                                                         DTOutput("prediction_table")
                                                                ),
                                                                
                                                                
                                                                tabPanel("Gene Analysis",
                                                                         h1("Select Genes to Visualize"),
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
                                                                
                                                                tabPanel("Performance", value = "Performance",
                                                                         div(class = "browse-instruction",
                                                                             h3("How to Use the Performance Section"),
                                                                             p("This section lets you evaluate how well your selected model performed in predicting tumour grades."),
                                                                             tags$ul(
                                                                               tags$li("If you have selected patients from the Prediction Table, performance plots will reflect those selected patients."),
                                                                               tags$li("If no selection is made, you can check the box below to view model performance on the full dataset.")
                                                                             )
                                                                         ),
                                                                         
                                                                         checkboxInput("use_all_data", "Use full dataset (ignore selection)", value = FALSE),
                                                                         
                                                                         h1("Grade Distribution"),
                                                                         p("You can choose between two visualisation styles to explore prediction performance across grades."),
                                                                         
                                                                         tags$ul(
                                                                           tags$li(strong("Bar + Line Plot:"), " Bars show predicted grade counts, and the line represents the true grade counts, making it easy to compare them directly."),
                                                                           tags$li(strong("Confusion matrix heatmap:"), " Displays true vs predicted grades in a grid format, helping you spot correct and incorrect predictions clearly.")
                                                                         ),
                                                                         p("These plots provide an intuitive and visually distinct summary of model performance."),
                                                                         
                                                                       
                                                                         conditionalPanel(
                                                                           condition = "input.module3_lo5_subtabs == 'Performance'",
                                                                           radioButtons("bar_plot_style", "Plot Style:",
                                                                                        choices = c("Bar + Line Plot" = "facet", "Confusion Matrix Heatmap" = "heatmap"),
                                                                                        inline = TRUE),
                                                                           div(style = "display: flex; justify-content: center;",
                                                                               plotOutput("bar_plot", width = "600px", height = "350px"))
                                                                         ),
                                                                         
                                                                         
                                                                         conditionalPanel(
                                                                           condition = "input.module3_lo5_subtabs == 'Performance'",
                                                                           h1("Heatmap"),
                                                                           p("This heatmap visualises the expression levels of the top 10 most variable genes across selected patients."),
                                                                           div(style = "display: flex; justify-content: center;",
                                                                               plotOutput("heatmap_plot", width = "750px", height = "400px"))
                                                                         )
                                                                )
                                                    )
                                                    )
                                           ),
                      
                                  tabPanel("Check Your Understanding",
                                           div(class = "nav-tabs sub-tab",
                                               tabsetPanel(id = "module3_quiz_subtabs", type = "tabs",
                                                           tabPanel("Quiz", uiOutput("quiz_game_ui_m3")),
                                                           tabPanel("Review & Summary", uiOutput("quiz_summary_ui_m3"))
                                               )
                                           )
                                  )
                      ))
           },
           
           "Module 4" = {
             tags$div(class = "lo-tab",
                      tags$div(class = "lo-tab",
                               tabsetPanel(id = "module4_tabs", type = "tabs",
                                           tabPanel("LO6: Model Evaluation",
                                                    tags$div(class = "sub-tab",
                                                             tabsetPanel(id = "module4_lo6_subtabs", type = "tabs",
                                                                         tabPanel("Glossary", 
                                                                                  div(class = "browse-instruction",
                                                                                      h3("Understanding Preva Terminology"),
                                                                                      p("This glossary provides definitions for key terms, metrics, and concepts used throughout the Preva application."),
                                                                                      tags$ul(
                                                                                        tags$li(tags$b("True Grade:"), " The actual tumour grade recorded in the dataset."),
                                                                                        tags$li(tags$b("Predicted Grade:"), " The grade predicted by the selected machine learning model."),
                                                                                        tags$li(tags$b("Chance / Probability:"), " The model's confidence score for each possible grade category."),
                                                                                        tags$li(tags$b("SMOTE:"), " A method to address class imbalance by synthesising new examples of underrepresented classes."),
                                                                                        tags$li(tags$b("Batch Correction:"), " A technique used to remove technical differences between datasets (e.g., ComBat).")
                                                                                      ),
                                                                                      p("Use this reference to better interpret model results, performance metrics, and gene visualisations throughout the app.")
                                                                                  ),
                                                                                  h1("Glossary of Model Evaluation Terms"),
                                                                                  p("This glossary explains key terms used in evaluating machine learning models in Preva."),
                                                                                  DT::dataTableOutput("glossary_table")
                                                                         ),
                                                                         
                                                                         tabPanel("F score", 
                                                                                  h1("F1 Score Summary"),
                                                                                  p("The F1 score combines precision and recall for each class, and is especially useful for evaluating imbalanced datasets. ",
                                                                                    "This table shows the per-grade classification metrics for the selected model and dataset."),
                                                                                  p(strong("Precision:"), "Out of all predicted positive samples for a class, how many were actually correct."),
                                                                                  p(strong("Recall:"), "Out of all actual positive samples for a class, how many were correctly identified."),
                                                                                  p(strong("F1 Score:"), "The harmonic mean of precision and recall. A balanced measure of accuracy for each class."),
                                                                                  DT::dataTableOutput("fscore_table")
                                                                         ),
                                                                         tabPanel("ROC / AUC", 
                                                                                  h1("Receiver Operating Characteristic (ROC) Curve"),
                                                                                  
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
                                                                                      plotOutput("roc_plot", width = "600px", height = "400px")
                                                                                  )
                                                                                  
                                                                                  
                                                                         ),
                                                                         tabPanel("Confusion Matrix",
                                                                                  div(class = "browse-instruction",
                                                                                      h3("How to Read the Confusion Matrix"),
                                                                                      p("This matrix compares predicted grades (columns) to true grades (rows). Each cell shows the number of samples in that category."),
                                                                                      tags$ul(
                                                                                        tags$li(tags$b("0 = Normal")),
                                                                                        tags$li(tags$b("1 = Grade 1")),
                                                                                        tags$li(tags$b("2 = Grade 2")),
                                                                                        tags$li(tags$b("3 = Grade 3"))
                                                                                      ),
                                                                                      p("Higher values on the diagonal indicate more accurate predictions.")
                                                                                  ),
                                                                                  h1("Confusion Matrix"),
                                                                                  p("This matrix shows how often predicted grades match the actual grades. Diagonal values are correct predictions, and off-diagonal cells represent errors."),
                                                                                  
                                                                                  
                                                                                  tableOutput("confusion_matrix_table"),
                                                                                  
                                                                                  h2("False Negative Rate (FNR) Explanation"),
                                                                                  p(strong("False Negative Rate (FNR):"),
                                                                                    "Also known as the miss rate, FNR is the proportion of actual positive cases that the model failed to detect. ",
                                                                                    "It is computed as 1 − Recall (Sensitivity)."),
                                                                                  p("FNR is important in clinical applications, where missing a diagnosis may be more serious than over-predicting one.")
                                                                         ),
                                                                         tabPanel("Cross Validation",
                                                                                  div(class = "browse-instruction",
                                                                                      h3("Cross-Validation Performance"),
                                                                                      p("This section shows model performance across 5-fold cross-validation."),
                                                                                      tags$ul(
                                                                                        tags$li("Each fold uses a different split of training/testing data."),
                                                                                        tags$li("Performance metrics such as accuracy and F1-score are averaged across folds."),
                                                                                        tags$li("Helps estimate how well the model will generalise to unseen data.")
                                                                                      )
                                                                                  ),
                                                                                  h1("Cross Validation Overview"),
                                                                                  p("Cross-validation is a technique used to estimate the generalisation performance of a machine learning model. ",
                                                                                    "It works by splitting the data into multiple training/testing folds and averaging performance across folds."),
                                                                                  p("This helps prevent overfitting and provides a more robust measure of model accuracy."),
                                                                                  p("In future versions, this section will include visualisations of accuracy, F1 score, and other metrics across folds.")
                                                                         ),
                                                                         tabPanel("Metric Summary", 
                                                                                  h1("Evaluation Summary Across Gene Sets"),
                                                                                  p("This plot shows how model performance changes across different gene sets. Currently shown for demonstration only."),
                                                                                  DT::dataTableOutput("metric_table")
                                                                                  
                                                                         ),
                                                                         tabPanel("Model Comparison",
                                                                                  h1("Model Comparison"),
                                                                                  p("This section compares the predictive performance of different machine learning models across various gene selection strategies and cross-validation settings."),
                                                                                  p("You can choose the number of top genes (from the selector above) and the cross-validation scheme (5-fold or 10-fold) to explore how model performance changes."),
                                                                                  p(strong("Note:"), " We compare SVM, KNN and RF side-by-side, ignoring the global Model selector."),
                                                                                  
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           radioButtons("comp_cv", "CV Fold:",
                                                                                                        choices = c("5-fold" = "CV5"),
                                                                                                        selected = "CV5")
                                                                                    ),
                                                                                    column(8,
                                                                                           plotOutput("comparison_plot"),
                                                                                           verbatimTextOutput("mean_accuracy_text")
                                                                                    )
                                                                                  )
                                                                         )
                                                                         
                                                             ))
                                                    
                                           ),
                                           
                                           tabPanel("LO7: Communicating Findings", 
                                                    tags$div(class = "sub-tab",
                                                             tabsetPanel(id = "module4_lo7_subtabs", type = "tabs",
                                                                         tabPanel("Key Findings Summary", value = "Key Findings Summary",
                                                                                  h1("Summary of Model and Dataset Findings"),
                                                                                  
                                                                                  p("Preva evaluated predictive models across multiple breast cancer datasets. Among the models tested, Random Forest (RF) consistently showed the best performance in terms of AUC and accuracy, especially when using a smaller subset of top-ranked genes (e.g., top 30)."),
                                                                                  
                                                                                  p("K-Nearest Neighbours (KNN) ranked second in most evaluations, and Support Vector Machines (SVM) third. RF’s strong performance likely reflects its ensemble nature and ability to handle noisy data. While SVM is robust for certain tasks, it performed less consistently across all datasets."),
                                                                                  
                                                                                  h3("Top 5 Genes and Their Roles"),
                                                                                  p("The top 5 genes selected based on feature importance scores showed consistent relevance across multiple datasets and models. These genes can be further investigated in biological databases (e.g., UniProt, PubMed) to explore their functional roles in tumour progression and grade differentiation."),
                                                                                  
                                                                                  h3("Dataset Observations"),
                                                                                  tags$ul(
                                                                                    tags$li("GSE15852 has 43 tumour and 43 normal samples—balanced but relatively small."),
                                                                                    tags$li("GSE10810 includes 31 tumour and 27 control samples—also balanced but limited."),
                                                                                    tags$li("GSE17907 consists of 340 tumour samples without controls—useful for high-grade training, but lacks comparative baseline."),
                                                                                    tags$li("The Combined dataset merges GSE10810 and GSE17907 but is biased due to overrepresentation of ERBB2-positive samples in GSE17907.")
                                                                                  )
                                                                         ),
                                                                         
                                                                         tabPanel(
                                                                           "Limitations", 
                                                                           h1("Limitations"),
                                                                           p("Limitations of current study and modelling process."),
                                                                           tags$ul(
                                                                             tags$li(
                                                                               tags$strong("Uneven sample size and classification distribution"),
                                                                               tags$p("Although we oversampled a few categories through SMOTE, synthetic samples cannot completely replace real biological samples and may not be able to capture all biological variations.")
                                                                             ),
                                                                             tags$li(
                                                                               tags$strong("Limitations of feature selection methods"),
                                                                               tags$p("We only use the limma strategy based on the linear model to select the Top N genes. However, the relationship between genes and grade is not necessarily linear; it may be nonlinear.")
                                                                             ),
                                                                             tags$li(
                                                                               tags$strong("Limitations of the model"),
                                                                               tags$p("Although models such as SVM, KNN, and RF can provide better classification performance, they do not directly give the contribution of a single gene to the prediction.")
                                                                             )
                                                                           )
                                                                         ),
                                                                         
                                                                         tabPanel("Future Directions", value = "Future Directions",
                                                                                  h1("Future Improvements and Educational Use"),
                                                                                  
                                                                                  p("To further improve Preva and its impact, both technically and educationally, the following enhancements are proposed:"),
                                                                                  
                                                                                  tags$ul(
                                                                                    tags$li("Train models on more diverse and larger datasets, including more samples with matched normal controls to improve comparison."),
                                                                                    tags$li("Support user-adjustable model parameters, such as number of trees in RF or k-value in KNN, to allow exploratory tuning and hands-on learning."),
                                                                                    tags$li("Expand gene expression analysis to other cancer types and consider clinical deployment potential."),
                                                                                    tags$li("Integrate in-app educational content to explain metrics, gene biology, and model function for students and early researchers."),
                                                                                    tags$li("Collaborate with hospitals or research labs to gain access to more labelled data and validate predictions in real-world contexts.")
                                                                                  )
                                                                         )
                                                                         
                                                                         
                                                             ))),
                                           
                                           tabPanel("Check Your Understanding",
                                                    div(class = "nav-tabs sub-tab",
                                                        tabsetPanel(id = "module4_quiz_subtabs", type = "tabs",
                                                                    tabPanel("Quiz", uiOutput("quiz_game_ui_m4")),
                                                                    tabPanel("Review & Summary", uiOutput("quiz_summary_ui_m4"))
                                                        )
                                                    )
                                           )
                               )))}
           
    )
  })
  
  
  
  output$confusion_matrix_table <- renderTable({
    ds      <- input$dataset_choice
    model   <- input$model_choice
    geneset <- paste0("top", input$topn_choice)
    df      <- all_results_df[[ds]][[geneset]][[model]]
    req(df)
    
    cm <- caret::confusionMatrix(
      factor(df$"Predicted Grade", levels = c("0","1","2","3")),
      factor(df$"True Grade",      levels = c("0","1","2","3"))
    )
    
    as.data.frame.matrix(cm$table)
  }, rownames = TRUE)
  
  
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
    req(input$dataset_choice, input$model_choice, input$topn_choice)
    ds    <- input$dataset_choice
    model <- input$model_choice
    gs    <- paste0("top", input$topn_choice)
    
    if (ds == "Uploaded") {
      df <- predict_uploaded_data(uploaded_exprs(), model)
    } else {
      # 一样从 all_results_df[[ds]][[gs]][[model]] 拿
      df <- all_results_df[[ds]][[gs]][[model]]
    }
    
    datatable(df, options = list(pageLength = 10))
  })
  
  
  observeEvent(input$dataset_choice, {
    if (page() == "Module 3") {
      isolate({
        updateTabsetPanel(session, "modelling_tab", selected = input$modelling_tab)
      })
    } else if (page() == "Module 4") {
      isolate({
        updateTabsetPanel(session, "evaluation_tab", selected = input$evaluation_tab)
      })
    }
  })
  
  
  observeEvent(input$select_all, {
    req(input$dataset_choice, input$model_choice, input$topn_choice)
    ds    <- input$dataset_choice
    gs    <- paste0("top", input$topn_choice)  # geneset 名称
    model <- input$model_choice
    
    df <- all_results_df[[ds]][[gs]][[model]]
    req(!is.null(df))
    
    proxy <- dataTableProxy("prediction_table")
    selectRows(proxy, seq_len(nrow(df)))
  })
  
  observeEvent(input$clear_all, {
    proxy <- dataTableProxy("prediction_table")
    selectRows(proxy, NULL)
  })
  
  

  #Combined Comparison
  observeEvent(input$dataset_choice, {
    if (input$dataset_choice == "Combined" && is.null(comparison_results_combined())) {
      showModal(modalDialog(
        title = "Preparing model comparison for Combined dataset...",
        "This may take a moment.",
        easyClose = FALSE,
        footer = NULL
      ))
      
      result <- prepare_predictions_for("Combined")  
      comparison_results_combined(result$models)
      
      removeModal()
    }
  })
  
  output$tbl_GSE15852 <- renderTable({
    data.frame(
      Metric = c("Genes",   "Normal", "Grade1", "Grade2", "Grade3"),
      Count  = c(22282,     43,       8,        23,       12),
      stringsAsFactors = FALSE
    )
  }, colnames = TRUE, rownames = FALSE, digits = 0)
  
  output$tbl_GSE17907 <- renderTable({
    data.frame(
      Metric = c("Genes",   "Normal", "Grade1", "Grade2", "Grade3"),
      Count  = c(24577,     8,        3,        10,       34),
      stringsAsFactors = FALSE
    )
  }, colnames = TRUE, rownames = FALSE, digits = 0)
  
  output$tbl_Combined <- renderTable({
    data.frame(
      Metric = c("Genes",   "Normal", "Grade1", "Grade2", "Grade3"),
      Count  = c(16760,     35,       5,        20,       44),
      stringsAsFactors = FALSE
    )
  }, colnames = TRUE, rownames = FALSE, digits = 0)
  
  output$pca_plot <- renderPlot({
    df <- load_data(input$dataset_choice)
    if (nrow(df) == 0 || !"Grade" %in% colnames(df) || all(is.na(df$Grade))) {
      plot.new()
      text(0.5, 0.5, "No data available for PCA.", cex = 1.3)
      return()
    }
    
    expr_mat <- df[, sapply(df, is.numeric), drop = FALSE]
    rownames(expr_mat) <- df$Patient
    grade     <- df$Grade
    
    gene_set_name <- paste0("top", input$topn_choice)    # e.g. "top30"
    genesets      <- all_genesets_df[[input$dataset_choice]]
    if (!gene_set_name %in% names(genesets)) {
      plot.new()
      text(0.5, 0.5, "No gene set found for PCA.", cex = 1.3)
      return()
    }
    selected_genes <- intersect(genesets[[gene_set_name]], colnames(expr_mat))
    if (length(selected_genes) < 2) {
      plot.new()
      text(0.5, 0.5, "Not enough genes for PCA.", cex = 1.3)
      return()
    }
    
    x <- expr_mat[, selected_genes, drop = FALSE]
    
    pca    <- prcomp(x, scale. = TRUE)
    pc_df  <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Grade = grade)
    pc_var <- round(summary(pca)$importance[2,1:2] * 100, 1)
    
    library(ggplot2)
    gg <- ggplot(pc_df, aes(x = PC1, y = PC2, color = Grade)) +
      geom_point(size = 3, alpha = 0.8) +
      theme_minimal(base_size = 13) +
      labs(
        title = "PCA of Samples (by Grade)",
        x     = paste0("PC1 (", pc_var[1], "%)"),
        y     = paste0("PC2 (", pc_var[2], "% variance)")
      )
    
    library(dplyr)
    grade_counts <- pc_df %>% count(Grade)
    valid_grades <- grade_counts$Grade[grade_counts$n >= 3]
    if (length(valid_grades) > 0) {
      gg <- gg + stat_ellipse(
        data     = pc_df %>% filter(Grade %in% valid_grades),
        type     = "norm",
        linetype = "dashed"
      )
    }
    
    print(gg)
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
    ds    <- input$dataset_choice
    genes <- paste0("top", input$topn_choice)
    mod   <- input$model_choice
    
    df <- all_results_df[[ds]][[genes]][[mod]]
    req(df)
    
    req("True Grade" %in% colnames(df),
        "Predicted Grade" %in% colnames(df))
    
    selected <- input$prediction_table_rows_selected
    if (!isTRUE(input$use_all_data)) {
      if (length(selected) == 0) return(NULL)
      df <- df[selected, ]
    }
    
    df$`True Grade`      <- factor(df$`True Grade`,      levels = c("0","1","2","3"))
    df$`Predicted Grade` <- factor(df$`Predicted Grade`, levels = c("0","1","2","3"))
    
    if (input$bar_plot_style == "facet") {
      library(dplyr)
      pred_count <- df %>%
        count(`Predicted Grade`) %>%
        rename(Grade = `Predicted Grade`, PredictedCount = n)
      
      true_count <- df %>%
        count(`True Grade`) %>%
        rename(Grade = `True Grade`, TrueCount = n)
      
      plot_df <- full_join(pred_count, true_count, by = "Grade") %>%
        replace_na(list(PredictedCount = 0, TrueCount = 0)) %>%
        mutate(Grade = factor(Grade, levels = c("0","1","2","3")))
      
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
        count(`True Grade`, `Predicted Grade`)
      
      ggplot(conf_mat, aes(x = `True Grade`, y = `Predicted Grade`, fill = n)) +
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
    ds      <- input$dataset_choice
    model   <- input$model_choice
    topn    <- input$topn_choice
    geneset <- paste0("top", topn)
    
    # 1. Get the data.frame containing True/Predicted Grade
    df <- all_results_df[[ds]][[geneset]][[model]]
    req(df)
    
    # 2. Calculate the confusion matrix
    cm <- caret::confusionMatrix(
      factor(df$"Predicted Grade", levels = c("0","1","2","3")),
      factor(df$"True Grade",      levels = c("0","1","2","3"))
    )
    
    # 3. Extract the metrics from byClass and only keep Precision/Recall/F1
    stats <- as.data.frame(cm$byClass)[, c("Precision","Recall","F1")]
    
    # 4. Replace all NA (i.e., null values generated when there are no samples or no predictions) with 0
    stats[is.na(stats)] <- 0
    
    
    stats <- round(stats, 3)
    stats$Grade <- rownames(stats)
    stats <- stats[, c("Grade","Precision","Recall","F1")]
    
    DT::datatable(
      stats,
      options = list(
        paging    = FALSE,
        searching = FALSE,
        info      = FALSE,
        ordering  = FALSE,
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  })
  
  output$metric_table <- DT::renderDataTable({
    ds      <- input$dataset_choice
    model   <- input$model_choice
    geneset <- paste0("top", input$topn_choice)
    
    df <- all_results_df[[ds]][[geneset]][[model]]
    req(df)
    
    cm    <- caret::confusionMatrix(
      factor(df$`Predicted Grade`, levels = c("0","1","2","3")),
      factor(df$`True Grade`,      levels = c("0","1","2","3"))
    )
    kappa <- round(cm$overall["Kappa"], 3)
    
    tbl <- data.frame(
      Dataset = ds,
      Model   = model,
      Geneset = geneset,
      Kappa   = kappa,
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      tbl,
      options = list(
        paging    = FALSE,
        searching = FALSE,
        info      = FALSE,
        ordering  = FALSE,
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  })
  
  
  output$roc_plot <- renderPlot({
    
    ds      <- input$dataset_choice
    model   <- input$model_choice
    geneset <- paste0("top", input$topn_choice)
    df      <- all_results_df[[ds]][[geneset]][[model]]
    req(df)
    
    true <- factor(df$`True Grade`, levels = c("0","1","2","3"))
    prob_df <- df %>%
      select(ends_with("Chance")) %>%
      lapply(function(x) as.numeric(sub("%","",x))/100) %>%
      as.data.frame()
    colnames(prob_df) <- c("0","1","2","3")
    
    plot(NULL, xlim=c(0,1), ylim=c(0,1),
         xlab="False Positive Rate", ylab="True Positive Rate",
         main=paste("ROC Curves –", model, "on", ds),
         cex.lab=1.2, cex.main=1.3)
    abline(0,1, lty=2, col="gray")
    
    classes <- levels(true)
    valid  <- sapply(classes, function(cls) {
      vals <- ifelse(true==cls, 1, 0)
      length(unique(vals)) == 2
    })
    valid_classes <- classes[valid]
    colors <- c("#c8b0d8","#e3a5b2","#a2b5c4","#b1d1c5")[valid]
    
    aucs <- numeric(length(valid_classes))
    for (i in seq_along(valid_classes)) {
      cls      <- valid_classes[i]
      bin_true <- ifelse(true==cls, 1, 0)
      roc_obj  <- roc(bin_true, prob_df[[cls]], quiet=TRUE)
      aucs[i]  <- round(auc(roc_obj), 3)
      lines(roc_obj, col=colors[i], lwd=2)
    }
    
    legend_txt <- paste0("Grade ", valid_classes, " (AUC=", aucs, ")")
    legend("bottomleft", legend=legend_txt, col=colors, lwd=2, bty="n", cex=1)
    
  }, width=600, height=400)
  
  
  
  output$comparison_plot <- renderPlot({
    req(input$dataset_choice, input$topn_choice, input$comp_cv)
    ds      <- input$dataset_choice
    geneset <- paste0("top", input$topn_choice) 
    cv_n    <- as.numeric(sub("CV", "", input$comp_cv))  
    
    models <- if (ds == "Combined") {
      comparison_results_combined()[[geneset]]
    } else {
      all_model_list[[ds]][[geneset]]
    }
    
    model_names <- names(models)  # c("SVM","KNN","RF")
    
    df_list <- lapply(model_names, function(m) {
      resamp <- models[[m]]$resample
      data.frame(Model = m, Accuracy = resamp$Accuracy)
    })
    df <- do.call(rbind, df_list)
    req(nrow(df)>0)
    
    ggplot(df, aes(x = Model, y = Accuracy, fill = Model)) +
      geom_boxplot(width = 0.6, alpha = 0.8) +
      theme_minimal(base_size = 13) +
      labs(
        title = paste("Model Accuracy Comparison —", geneset, input$comp_cv),
        x = "Model", y = "Accuracy"
      ) +
      scale_fill_manual(values = c("SVM"="#c8b0d8","KNN"="#e3a5b2","RF"="#b1d1c5")) +
      theme(legend.position = "none")
  })
  
  output$mean_accuracy_text <- renderPrint({
    req(input$dataset_choice, input$topn_choice, input$comp_cv)
    ds      <- input$dataset_choice
    geneset <- paste0("top", input$topn_choice)
    cv_n    <- as.numeric(sub("CV", "", input$comp_cv))
    models <- if (ds == "Combined") {
      comparison_results_combined()[[geneset]]
    } else {
      all_model_list[[ds]][[geneset]]
    }
    
    for(m in names(models)) {
      resamp <- models[[m]]$resample
      resamp <- resamp[grepl(paste0("^Fold", cv_n), resamp$Resample, ignore.case = TRUE), , drop = FALSE]
      cat(sprintf("%s: Mean Accuracy = %.3f\n", m, mean(resamp$Accuracy)))
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
  
  # Glossary
  output$glossary_table <- DT::renderDataTable({
    glossary_df <- data.frame(
      Term = c(
        "F1 Score",
        "ROC Curve",
        "AUC (Area Under Curve)",
        "Confusion Matrix",
        "False Negative Rate (FNR)",
        "Cross-Validation",
        "Kappa Statistic",
        "Model Accuracy Comparison"
      ),
      Definition = c(
        "The harmonic mean of precision and recall, providing a balanced measure of accuracy.",
        "A plot showing the trade-off between sensitivity (TPR) and 1 − specificity (FPR).",
        "A single number summarising the ROC curve; the closer to 1, the better the model.",
        "A matrix that shows the number of true vs predicted classifications for each class.",
        "The proportion of actual positive cases that the model failed to detect: 1 − Recall.",
        "A model evaluation method that splits the data into training/testing folds to reduce bias.",
        "A statistic that measures agreement between predicted and true labels beyond chance.",
        "A boxplot showing accuracy distributions across different models (SVM, KNN, RF) using a specific gene set and cross-validation setup (e.g., 5-fold). It helps compare model performance visually."
      ),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      glossary_df,
      options = list(pageLength = 10, autoWidth = TRUE),
      rownames = FALSE
    )
  })
  
  
  
  # Quiz!
  
  quiz_list_m1 <- list(
    list(
      text = "What is the main purpose of the Preva application?",
      options = c(
        "A. To cure breast cancer",
        "B. To sell breast cancer medication to hospitals",
        "C. To detect and understand breast cancer",
        "D. To replace doctors in diagnosis"
      ),
      answer = "C. To detect and understand breast cancer"
    ),
    list(
      text = "What technology does the Preva app use to make predictions?",
      options = c(
        "A. Quantum computing",
        "B. Machine learning models",
        "C. Data Visualisation",
        "D. Blockchain"
      ),
      answer = "B. Machine learning models"
    ),
    list(
      text = "What type of data is being used for the predictions?",
      options = c(
        "A. Image data",
        "B. Gene expression data",
        "C. Specimen Data",
        "D. Stock market data"
      ),
      answer = "B. Gene expression data"
    )
  )
  
  quiz_list_m2 <- list(
    list(
      text = "Which of the following is not a variable in the gene expression Dataset?",
      options = c(
        "A. TrueGrade",
        "B. PredictedGrade",
        "C. Patient",
        "D. Grade 6 Chance"
      ),
      answer = "D. Grade 6 Chance"
    ),
    list(
      text = "Change this Question due to we delete the dataset 10810",
      options = c("A. 99", "B. 104", "C. 199", "D. 86"),
      answer = "D. 86"
    ),
    list(
      text = "Which of the following is a type of histopathology?",
      options = c(
        "A. Infiltrating ductal carcinoma",
        "B. Fibrous adenoma",
        "C. Basal squamous carcinoma",
        "D. Lobular duct syndrome"
      ),
      answer = "A. Infiltrating ductal carcinoma"
    )
  )
  
  quiz_list_m3 <- list(
    list(
      text = "Which Model was NOT used?",
      options = c("A. SVM", "B. RF", "C. KNN", "D. CNN"),
      answer = "D. CNN"
    ),
    list(
      text = "Which breast cancer stage has the highest median expression level for Gene 201650_at?",
      options = c("A. Normal", "B. 1", "C. 2", "D. 3"),
      answer = "C. 2"
    ),
    list(
      text = "What is the grade 1 chance of patient GSM398103?",
      options = c("A. 16.4%", "B. 55.9%", "C. 20.6%", "D. 8.3%"),
      answer = "A. 16.4%"
    ),
    list(
      text = "Which model had the most correct predictions using the GSE15852 dataset?",
      options = c("A. Models equal", "B. SVM", "C. RF", "D. KNN"),
      answer = "B. SVM"
    )
  )
  
  quiz_list_m4 <- list(
    list(
      text = "Which model has the highest median accuracy for top110 genes in 5-fold CV?",
      options = c("A. SVM", "B. KNN", "C. RF", "D. All equal"),
      answer = "C. RF"
    ),
    list(
      text = "Does changing from top110 to top80 improve RF median accuracy?",
      options = c("A. Improve", "B. Worsen"),
      answer = "B. Worsen"
    ),
    list(
      text = paste(
        "Please do the following:",
        "1. Select the ‘GSE15852’ dataset.",
        "2. Select the SVM model.",
        "3. Go to the ‘ROC / AUC’ tab.",
        "Which grade is falsely predicted the most?",
        sep = "\n"
      ),
      options = c("A. Grade 0", "B. Grade 1", "C. Grade 2", "D. Grade 3"),
      answer = "A. Grade 0"
    )
  )
  
  
  # Module 1
  current_question_m1 <- reactiveVal(1)
  user_answers_m1 <- reactiveValues()
  
  # Module 2
  current_question_m2 <- reactiveVal(1)
  user_answers_m2 <- reactiveValues()
  
  # Module 3
  current_question_m3 <- reactiveVal(1)
  user_answers_m3 <- reactiveValues()
  
  # Module 4
  current_question_m4 <- reactiveVal(1)
  user_answers_m4 <- reactiveValues()
  
  # M1
  output$quiz_game_ui_m1 <- renderUI({
    idx <- current_question_m1()
    q <- quiz_list_m1[[idx]]
    progress_percent <- round((idx / length(quiz_list_m1)) * 100)
    
    tagList(
      div(style = "margin-bottom: 15px;",
          div(style = "width: 100%; background-color: #f1f1f1; height: 20px; border-radius: 10px;",
              div(style = paste0("width:", progress_percent, "%; height: 100%; background-color: #c8b0d8;"))
          ),
          p(paste0(progress_percent, "% Complete"), style = "text-align: center; font-weight: bold;")
      ),
      h4(paste0(idx, ". ", q$text)),
      radioButtons("user_answer_m1", label = NULL, choices = q$options, selected = character(0)),
      div(style = "display: flex; justify-content: space-between;",
          actionButton("prev_m1", "← Back"),
          span(paste(idx, "/", length(quiz_list_m1))),
          actionButton("next_m1", "Next →")
      ),
      br(),
      div(style = "text-align: right;",
          actionButton("submit_m1", "Submit Quiz", class = "btn btn-success")
      )
    )
  })
  
  output$quiz_summary_ui_m1 <- renderUI({
    tagList(
      if (!quiz_submitted_m1()) {
        div(style = "color: orange; font-weight: bold;",
            "⚠️ You haven't submitted your answers yet. Please go to Quiz tab and click 'Submit Quiz'.")
      },
      lapply(seq_along(quiz_list_m1), function(i) {
        q <- quiz_list_m1[[i]]
        ans <- user_answers_m1[[paste0("q", i)]]
        is_correct <- identical(ans, q$answer)
        
        div(
          h4(paste(i, ".", q$text)),
          p(strong("Your Answer:"), ans %||% "(no answer)"),
          p(strong("Correct Answer:"), q$answer),
          if (is.null(ans)) {
            span("⚠️ Skipped", style = "color: orange;")
          } else if (is_correct) {
            span("✔️ Correct", style = "color: green;")
          } else {
            span("❌ Incorrect", style = "color: red;")
          }
        )
      })
    )
  })
  
  # M2
  output$quiz_game_ui_m2 <- renderUI({
    idx <- current_question_m2()
    q <- quiz_list_m2[[idx]]
    progress_percent <- round((idx / length(quiz_list_m2)) * 100)
    
    tagList(
      div(style = "margin-bottom: 15px;",
          div(style = "width: 100%; background-color: #f1f1f1; height: 20px; border-radius: 10px;",
              div(style = paste0("width:", progress_percent, "%; height: 100%; background-color: #c8b0d8;"))
          ),
          p(paste0(progress_percent, "% Complete"), style = "text-align: center; font-weight: bold;")
      ),
      h4(paste0(idx, ". ", q$text)),
      radioButtons("user_answer_m2", label = NULL, choices = q$options, selected = character(0)),
      div(style = "display: flex; justify-content: space-between;",
          actionButton("prev_m2", "← Back"),
          span(paste(idx, "/", length(quiz_list_m1))),
          actionButton("next_m2", "Next →")
      ),
      br(),
      div(style = "text-align: right;",
          actionButton("submit_m2", "Submit Quiz", class = "btn btn-success")
      )
    )
  })
  
  output$quiz_summary_ui_m2 <- renderUI({
    tagList(
      if (!quiz_submitted_m2()) {
        div(style = "color: orange; font-weight: bold;",
            "⚠️ You haven't submitted your answers yet. Please go to Quiz tab and click 'Submit Quiz'.")
      },
      lapply(seq_along(quiz_list_m2), function(i) {
        q <- quiz_list_m2[[i]]
        ans <- user_answers_m2[[paste0("q", i)]]
        is_correct <- identical(ans, q$answer)
        
        div(
          h4(paste(i, ".", q$text)),
          p(strong("Your Answer:"), ans %||% "(no answer)"),
          p(strong("Correct Answer:"), q$answer),
          if (is.null(ans)) {
            span("⚠️ Skipped", style = "color: orange;")
          } else if (is_correct) {
            span("✔️ Correct", style = "color: green;")
          } else {
            span("❌ Incorrect", style = "color: red;")
          }
        )
      })
    )
  })
  
  # M3
  output$quiz_game_ui_m3 <- renderUI({
    idx <- current_question_m3()
    q <- quiz_list_m3[[idx]]
    progress_percent <- round((idx / length(quiz_list_m3)) * 100)
    
    tagList(
      div(style = "margin-bottom: 15px;",
          div(style = "width: 100%; background-color: #f1f1f1; height: 20px; border-radius: 10px;",
              div(style = paste0("width:", progress_percent, "%; height: 100%; background-color: #c8b0d8;"))
          ),
          p(paste0(progress_percent, "% Complete"), style = "text-align: center; font-weight: bold;")
      ),
      h4(paste0(idx, ". ", q$text)),
      radioButtons("user_answer_m3", label = NULL, choices = q$options, selected = character(0)),
      div(style = "display: flex; justify-content: space-between;",
          actionButton("prev_m3", "← Back"),
          span(paste(idx, "/", length(quiz_list_m1))),
          actionButton("next_m3", "Next →")
      ),
      br(),
      div(style = "text-align: right;",
          actionButton("submit_m3", "Submit Quiz", class = "btn btn-success")
      )
    )
  })
  
  output$quiz_summary_ui_m3<- renderUI({
    tagList(
      if (!quiz_submitted_m3()) {
        div(style = "color: orange; font-weight: bold;",
            "⚠️ You haven't submitted your answers yet. Please go to Quiz tab and click 'Submit Quiz'.")
      },
      lapply(seq_along(quiz_list_m3), function(i) {
        q <- quiz_list_m3[[i]]
        ans <- user_answers_m3[[paste0("q", i)]]
        is_correct <- identical(ans, q$answer)
        
        div(
          h4(paste(i, ".", q$text)),
          p(strong("Your Answer:"), ans %||% "(no answer)"),
          p(strong("Correct Answer:"), q$answer),
          if (is.null(ans)) {
            span("⚠️ Skipped", style = "color: orange;")
          } else if (is_correct) {
            span("✔️ Correct", style = "color: green;")
          } else {
            span("❌ Incorrect", style = "color: red;")
          }
        )
      })
    )
  })
  
  # M4
  output$quiz_game_ui_m4 <- renderUI({
    idx <- current_question_m4()
    q <- quiz_list_m4[[idx]]
    progress_percent <- round((idx / length(quiz_list_m4)) * 100)
    
    tagList(
      div(style = "margin-bottom: 15px;",
          div(style = "width: 100%; background-color: #f1f1f1; height: 20px; border-radius: 10px;",
              div(style = paste0("width:", progress_percent, "%; height: 100%; background-color: #c8b0d8;"))
          ),
          p(paste0(progress_percent, "% Complete"), style = "text-align: center; font-weight: bold;")
      ),
      h4(paste0(idx, ". ", q$text)),
      radioButtons("user_answer_m4", label = NULL, choices = q$options, selected = character(0)),
      div(style = "display: flex; justify-content: space-between;",
          actionButton("prev_m4", "← Back"),
          span(paste(idx, "/", length(quiz_list_m1))),
          actionButton("next_m4", "Next →")
      ),
      br(),
      div(style = "text-align: right;",
          actionButton("submit_m4", "Submit Quiz", class = "btn btn-success")
      )
    )
  })
  
  output$quiz_summary_ui_m4<- renderUI({
    tagList(
      if (!quiz_submitted_m4()) {
        div(style = "color: orange; font-weight: bold;",
            "⚠️ You haven't submitted your answers yet. Please go to Quiz tab and click 'Submit Quiz'.")
      },
      lapply(seq_along(quiz_list_m4), function(i) {
        q <- quiz_list_m4[[i]]
        ans <- user_answers_m3[[paste0("q", i)]]
        is_correct <- identical(ans, q$answer)
        
        div(
          h4(paste(i, ".", q$text)),
          p(strong("Your Answer:"), ans %||% "(no answer)"),
          p(strong("Correct Answer:"), q$answer),
          if (is.null(ans)) {
            span("⚠️ Skipped", style = "color: orange;")
          } else if (is_correct) {
            span("✔️ Correct", style = "color: green;")
          } else {
            span("❌ Incorrect", style = "color: red;")
          }
        )
      })
    )
  })
  
  
  # M1
  observeEvent(input$next_m1, {
    if (current_question_m1() < length(quiz_list_m1)) {
      current_question_m1(current_question_m1() + 1)
    }
  })
  
  observeEvent(input$prev_m1, {
    if (current_question_m1() > 1) {
      current_question_m1(current_question_m1() - 1)
    }
  })
  
  observeEvent(input$user_answer_m1, {
    idx <- current_question_m1()
    user_answers_m1[[paste0("q", idx)]] <- input$user_answer_m1
  })
  
  # M2
  observeEvent(input$next_m2, {
    if (current_question_m2() < length(quiz_list_m2)) {
      current_question_m2(current_question_m2() + 1)
    }
  })
  
  observeEvent(input$prev_m2, {
    if (current_question_m2() > 1) {
      current_question_m2(current_question_m2() - 1)
    }
  })
  
  observeEvent(input$user_answer_m2, {
    idx <- current_question_m2()
    user_answers_m2[[paste0("q", idx)]] <- input$user_answer_m2
  })
  
  # M3
  observeEvent(input$next_m3, {
    if (current_question_m3() < length(quiz_list_m3)) {
      current_question_m3(current_question_m3() + 1)
    }
  })
  
  observeEvent(input$prev_m3, {
    if (current_question_m3() > 1) {
      current_question_m3(current_question_m3() - 1)
    }
  })
  
  observeEvent(input$user_answer_m3, {
    idx <- current_question_m3()
    user_answers_m3[[paste0("q", idx)]] <- input$user_answer_m3
  })
  
  # M4
  observeEvent(input$next_m4, {
    if (current_question_m4() < length(quiz_list_m4)) {
      current_question_m4(current_question_m4() + 1)
    }
  })
  
  observeEvent(input$prev_m4, {
    if (current_question_m4() > 1) {
      current_question_m4(current_question_m4() - 1)
    }
  })
  
  observeEvent(input$user_answer_m4, {
    idx <- current_question_m4()
    user_answers_m4[[paste0("q", idx)]] <- input$user_answer_m4
  })
  
  # M1
  observeEvent(input$submit_m1, {
    quiz_submitted_m1(TRUE)
    showNotification("Your answers have been saved!", type = "message")
  })
  
  # M2
  observeEvent(input$submit_m2, {
    quiz_submitted_m2(TRUE)
    showNotification("Your answers have been saved!", type = "message")
  })
  
  
  # M3
  observeEvent(input$submit_m3, {
    quiz_submitted_m3(TRUE)
    showNotification("Your answers have been saved!", type = "message")
  })
  
  
  # M4
  observeEvent(input$submit_m4, {
    quiz_submitted_m4(TRUE)
    showNotification("Your answers have been saved!", type = "message")
  })
  
  
  
  
}  

shinyApp(ui, server)