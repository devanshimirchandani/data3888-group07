cat(">>> Sourcing server.R <<<\n")

options(shiny.error         = browser)
options(shiny.sanitize.errors = FALSE)
options(shiny.trace         = TRUE)

server <- function(input, output, session) {
  comparison_results <- reactiveVal(NULL)
  results_df <- prepare_predictions()
  
  data <- reactive(load_data())
  page <- reactiveVal("About")
  
  observeEvent(input$about_tab, page("About"))
  observeEvent(input$how_tab, page("How"))
  observeEvent(input$data_tab, page("Data"))
  observeEvent(input$analysis_tab, page("Analysis"))
  observeEvent(input$eval_tab, page("Evaluation"))
  observeEvent(TRUE, {
    if (is.null(comparison_results())) {
      result <- prepare_comparison_data()
      comparison_results(result)
      saveRDS(result, "comparison_cache.rds")  # 可选
    }
  }, once = TRUE)
  
  
  
  
  output$conditional_selectors <- renderUI({
    if (page() %in% c("Analysis", "Evaluation")) {
      fluidRow(
        column(4,
               selectInput("dataset_choice", "Dataset:",
                           choices = c("GSE15852", "Dataset 1", "Dataset 2"),
                           selected = "GSE15852")
        ),
        column(4,
               conditionalPanel(
                 condition = "input.dataset_choice == 'GSE15852'",
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
    if (input$dataset_choice != "GSE15852") {
      return(tagList(
        h3("Coming soon..."),
        p("This dataset is not yet available.")
      ))
    }
    
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
           "How" = {
             tagList(
               h2("How to Use"),
               p("Use the top navigation bar to explore:"),
               tags$ul(
                 tags$li("Data – view patient demographics and input values."),
                 tags$li("Analysis – inspect predictions, gene analysis, and comparisons."),
                 tags$li("Evaluation – PCA overview of gene expression profiles.")
               )
             )
           },
           "Data" = {
             tabsetPanel(
               tabPanel("Browse Dataset",
                        fileInput("browse_data", "Upload your own dataset (.RData)", accept = ".RData"),
                        DTOutput("browse_table")),
               tabPanel("GSM15852 Data", DTOutput("data_table"))
             )
           },
           
           "Analysis" = {
             tabsetPanel(
               tabPanel("Prediction",
                        h3("Prediction Table"),
                        conditionalPanel(
                          condition = "input.model_choice == 'SVM'",
                          actionButton("select_all", "Select All"),
                          actionButton("clear_all", "Clear Selection"),
                          DTOutput("prediction_table")
                        ),
                        conditionalPanel(
                          condition = "input.model_choice != 'SVM'",
                          h4("Coming soon...")
                        )
               ),
               tabPanel("Gene Analysis",
                        h3("Top 5 Variable Genes"),
                        plotOutput("gene_box", width = "600px", height = "400px")),
               tabPanel("Comparison",
                        h3("Grade Distribution"),
                        plotOutput("bar_plot", width = "600px", height = "400px"),
                        h3("Heatmap"),
                        plotOutput("heatmap_plot", width = "600px", height = "400px"),
                        h3("Correlation Plot"),
                        p("Coming soon..."))
             )
           },
           
           "Evaluation" = {
             tabsetPanel(
               tabPanel("PCA",
                        sliderInput("topN", "Number of top variable genes to include in PCA:",
                                    min = 30, max = 110, value = 80, step = 10),
                        radioButtons("cvChoice", "Cross-validation folds:",
                                     choices = c("5-fold" = 5, "10-fold" = 10),
                                     selected = 5, inline = TRUE),
                        plotOutput("pca_plot", width = "600px", height = "400px")
               ),
               tabPanel("F Score", verbatimTextOutput("f1_score_text")),
               tabPanel("Confusion Matrix", p("Coming soon...")),
               tabPanel("Cross Validation", p("Coming soon...")),
               tabPanel("ROC Curve / AUC", verbatimTextOutput("roc_auc_text")),
               tabPanel("Feature Importance", p("Coming soon...")),
               tabPanel("Comparison",
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
           }
           
    )
  })
  
  
  output$pca_plot <- renderPlot({
    # 一定要先拿到比较结果
    req(comparison_results())
    comp_data <- comparison_results()
    df        <- data()
    
    # 构造 key：比如 top80 + CV5
    top_key <- paste0("top", input$topN)      # "top30","top50","top80" 或 "top110"
    cv_key  <- paste0("CV", input$cvChoice)   # "CV5" 或 "CV10"
    
    # 从缓存里取出基因向量
    genes <- comp_data[[top_key]][[cv_key]]$Genes
    
    # 准备 PCA 输入，只保留这批基因
    pca_df_input <- df[, genes, drop = FALSE]
    # scale PCA
    pca_res <- prcomp(pca_df_input, scale. = TRUE)
    
    # 整理绘图数据
    pca_plot_df <- data.frame(
      PC1   = pca_res$x[,1],
      PC2   = pca_res$x[,2],
      Grade = df$Grade
    )
    # 计算各主成分方差百分比
    pc_var <- round(summary(pca_res)$importance[2,1:2] * 100, 1)
    
    # 绘图
    ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Grade)) +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(type = "norm", linetype = "dashed") +
      theme_minimal() +
      labs(
        title = paste0("PCA (", input$cvChoice, "-fold CV, top ", input$topN, " genes)"),
        x     = paste0("PC1 (", pc_var[1], "%)"),
        y     = paste0("PC2 (", pc_var[2], "%)")
      )
  })
  
  
  
  
  
  
  output$data_table <- renderDT({
    df <- data()
    datatable(df[, c("Patient", "Grade", "Histopathology")],
              selection = "multiple",
              options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10),
              rownames = FALSE)
  })
  
  output$browse_table <- renderDT({
    req(input$browse_data)
    inFile <- input$browse_data
    e <- new.env()
    load(inFile$datapath, envir = e)
    df <- as.data.frame(e[[ls(e)[1]]])
    datatable(df, options = list(pageLength = 10), rownames = TRUE)
  })
  
  output$prediction_table <- renderDT({
    datatable(results_df,
              selection = "multiple",
              options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10),
              rownames = TRUE)
  })
  
  observeEvent(input$select_all, {
    proxy <- dataTableProxy("prediction_table")
    selectRows(proxy, 1:nrow(results_df))
  })
  observeEvent(input$clear_all, {
    proxy <- dataTableProxy("prediction_table")
    selectRows(proxy, NULL)
  })
  
  output$gene_box <- renderPlot({
    df <- data()
    selected <- input$prediction_table_rows_selected
    if (length(selected) == 0) return(NULL)
    df <- df[selected, ]
    gene_vars <- apply(df[, 6:ncol(df)], 2, var)
    top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:5])
    expr_long <- melt(df[, c("Grade", top_genes)])
    ggplot(expr_long, aes(x = Grade, y = value, fill = Grade)) +
      geom_boxplot(alpha = 0.6) +
      facet_wrap(~variable, scales = "free") +
      theme_minimal(base_size = 10) +  
      theme(
        plot.margin = margin(5, 5, 5, 5),
        legend.position = "bottom",
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10)
      ) +
      labs(x = "Grade", y = "Expression")
  })
  
  
  output$bar_plot <- renderPlot({
    df <- data()
    selected <- input$prediction_table_rows_selected
    if (length(selected) > 0) df <- df[selected, ]
    ggplot(df, aes(x = Grade, fill = Grade)) +
      geom_bar(width = 0.6) +
      scale_fill_manual(values = c("#fcf8f3", "#b0c4b1", "#c8b0d8", "#e3a5b2")) +
      theme_minimal(base_size = 12) +
      theme(plot.margin = margin(10, 10, 10, 10))
  })
  
  
  output$heatmap_plot <- renderPlot({
    df <- data()
    selected <- input$prediction_table_rows_selected
    if (length(selected) == 0) return(NULL)
    gene_vars <- apply(df[, 6:ncol(df)], 2, var)
    top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:10])
    expr <- as.matrix(df[selected, top_genes])
    rownames(expr) <- df$Patient[selected]
    pheatmap(t(expr), cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 10,
             color = colorRampPalette(c("#e3a5b2", "#fcf8f3", "#c8b0d8"))(50))
  })
  
  
  output$corr_plot <- renderPlot({
    plot.new()
    text(0.5, 0.5, "Correlation Plot coming soon!", cex = 1.6)
  })
  
  
  output$f1_score_text <- renderPrint({
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
    
    trctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE,
                           summaryFunction = multiClassSummary, savePredictions = "final")
    
    set.seed(123)
    svm_mod <- train(x = x_train, y = labels_train, method = "svmLinear", trControl = trctrl)
    
    pred_test <- predict(svm_mod, newdata = x_test)
    cm <- confusionMatrix(pred_test, labels_test)
    
    f1_scores <- cm$byClass[, "F1"]
    macro_f1 <- mean(f1_scores, na.rm = TRUE)
    
    cat("Macro F1 Score:", round(macro_f1, 3), "\n\n")
    print(round(f1_scores, 3))
  })
  
  output$roc_auc_text <- renderPrint({
    library(pROC)
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
    
    trctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE,
                           summaryFunction = multiClassSummary, savePredictions = "final")
    
    set.seed(123)
    svm_mod <- train(x = x_train, y = labels_train, method = "svmLinear", trControl = trctrl)
    
    predicted_probs <- predict(svm_mod, newdata = x_test, type = "prob")
    auc_list <- c()
    for (label in colnames(predicted_probs)) {
      roc_obj <- roc(response = as.numeric(labels_test == label),
                     predictor = predicted_probs[[label]])
      auc_list[label] <- auc(roc_obj)
    }
    
    cat("AUC per class:\n")
    print(round(auc_list, 3))
  })
  
  
  set.seed(123)
  output$comparison_plot <- renderPlot({
    req(comparison_results())
    comp_data <- comparison_results()
    
    top_choice <- input$comp_topn
    cv_choice <- input$comp_cv
    models <- c("SVM", "KNN", "RF")
    
    acc_data <- do.call(rbind, lapply(models, function(model) {
      acc <- comp_data[[top_choice]][[cv_choice]][[model]]
      data.frame(Model = model, Accuracy = acc)
    }))
    
    mean_acc <- aggregate(Accuracy ~ Model, acc_data, mean)
    mean_acc$Label <- paste0("Mean: ", round(mean_acc$Accuracy * 100, 1), "%")
    
    ggplot(acc_data, aes(x = Model, y = Accuracy, fill = Model)) +
      geom_boxplot(alpha = 0.7) +
      geom_text(data = mean_acc, aes(label = Label, y = 1.01), size = 4) +
      theme_minimal() +
      ylim(0.6, 1.05) +
      labs(title = paste0("Model Accuracy (", cv_choice, ", ", top_choice, ")"),
           y = "Accuracy", x = "") +
      theme(legend.position = "none")
  })
  
  
  
  
}  