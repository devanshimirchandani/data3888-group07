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

load_data <- function() {
  gse <- getGEO("GSE15852", GSEMatrix = TRUE)[[1]]
  exprs_data <- exprs(gse)
  log_exprs <- log2(exprs_data + 1)
  pdata <- pData(gse)
  
  df <- as.data.frame(t(log_exprs))
  df$Patient <- rownames(df)
  df$Grade <- gsub("grade[:]* *", "", tolower(pdata$`grade:ch1`))
  df$Grade <- ifelse(df$Grade == "normal", "Normal", gsub("grade ", "", df$Grade))
  df$Grade <- factor(df$Grade, levels = c("Normal", "1", "2", "3"))
  df$Age <- as.numeric(as.character(pdata$`age in years:ch1`))
  df$Histopathology <- pdata$`histopathological exam:ch1`
  df$Race <- pdata$`race:ch1`
  
  df <- df[, c("Patient", "Grade", "Age", "Histopathology", "Race", setdiff(names(df), c("Patient", "Grade", "Age", "Histopathology", "Race")))]
  return(df)
}

prepare_predictions <- function() {
  gse <- getGEO("GSE15852", GSEMatrix = TRUE)[[1]]
  M <- log2(exprs(gse) + 1)
  p <- pData(gse)
  grade_vec <- p$`grade:ch1`
  grade_fac <- factor(grade_vec, levels = c("normal", "grade 1", "grade 2", "grade 3"))
  
  set.seed(123)
  train_idx <- createDataPartition(grade_fac, p = 0.8, list = FALSE)
  expr_train <- M[, train_idx]
  expr_test <- M[, -train_idx]
  labels_train <- factor(grade_fac[train_idx], levels = levels(grade_fac))
  labels_test <- factor(grade_fac[-train_idx], levels = levels(grade_fac))
  
  levels(labels_train) <- make.names(levels(labels_train))
  levels(labels_test) <- make.names(levels(labels_test))
  
  grade_map <- c("normal"=0, "grade.1"=1, "grade.2"=2, "grade.3"=3)
  grade_num <- grade_map[as.character(labels_train)]
  design_trend <- model.matrix(~ grade_num)
  fit_trend <- lmFit(expr_train, design_trend)
  fit_trend <- eBayes(fit_trend)
  res_trend <- topTable(fit_trend, coef = "grade_num", number = Inf, adjust.method = "BH")
  sig_genes <- rownames(subset(res_trend, adj.P.Val < 0.05 & abs(logFC) > 1))
  
  x_train <- t(expr_train[sig_genes, , drop = FALSE])
  x_test <- t(expr_test[sig_genes, , drop = FALSE])
  
  trctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                         summaryFunction = multiClassSummary, savePredictions = "final")
  set.seed(123)
  svm_mod <- train(x = x_train, y = labels_train, method = "svmLinear",
                   trControl = trctrl, probability = TRUE)
  
  predicted <- predict(svm_mod, newdata = x_test)
  predicted_probs <- predict(svm_mod, newdata = x_test, type = "prob")
  predicted_probs_percent <- round(predicted_probs * 100, 1)
  
  grade_label_map <- c("normal"="0", "grade.1"="1", "grade.2"="2", "grade.3"="3")
  results_df <- data.frame(
    Patient = rownames(x_test),
    TrueGrade = grade_label_map[as.character(labels_test)],
    PredictedGrade = grade_label_map[as.character(predicted)],
    `Normal Chance` = paste0(predicted_probs_percent$normal, "%"),
    `Grade 1 Chance` = paste0(predicted_probs_percent$grade.1, "%"),
    `Grade 2 Chance` = paste0(predicted_probs_percent$grade.2, "%"),
    `Grade 3 Chance` = paste0(predicted_probs_percent$grade.3, "%")
  )
  return(results_df)
}