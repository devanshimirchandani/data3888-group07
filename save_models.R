library(caret)
library(limma)
library(smotefamily)

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
    check.names = FALSE,   
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



load("GSE15852.RData")
load("data_GSE17907.RData")
load("data_Combined.RData")



cat("Training GSE15852...\n")
tmp1 <- prepare_predictions_all_models()
saveRDS(tmp1, "GSE15852_models.rds")

cat("Training GSE17907...\n")
tmp2 <- prepare_predictions_for("GSE17907")
saveRDS(tmp2, "GSE17907_models.rds")

cat("Training Combined...\n")
tmp3 <- prepare_predictions_for("Combined")
saveRDS(tmp3, "Combined_models.rds")

cat("Preparing comparison data...\n")
res <- prepare_comparison_data()
saveRDS(res, "comparison_cache.rds")

cat("All models saved as .rds files.\n")

