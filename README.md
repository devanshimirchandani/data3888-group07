# Preva — Breast Cancer Risk Calculator from Gene Expression Data
## Project Overview
Preva is an interactive risk calculator that predicts breast cancer severity from gene expression data.
This app was developed as a group project for DATA3888 (2025 S1) at the University of Sydney. (data3888-group07)

It is designed as a modular educational tool that guides users from raw data to prediction, helping students and learners explore key bioinformatics workflows in a hands-on way.

## Key Features
Module 1: Introduction to Breast Cancer & GED

Module 2: Data Upload and Processing (batch correction, resampling, filtering)

Module 3: Predictive Model Building (SVM, KNN, RF with tunable parameters)

Module 4: Risk Evaluation & Interpretation
(Accuracy, F1, Kappa, ROC/AUC, confusion matrix, top genes)

## Repository Structure

### **Shiny/**
Contains the full Shiny application, including UI, server logic, and pre-processed data used in the app.

Key components:

* [`app.R`](./Shiny/app.R) — Main file to launch the Preva Shiny application.

* [`www/`](./Shiny/www/) — Contains image resources (`.png`) used in the app.

* Pre-processed datasets:

  * [`GSE15852.RData`](./Shiny/GSE15852.RData)
  * [`data_GSE17907.RData`](./Shiny/data_GSE17907.RData)
  * [`data_GSE10810.RData`](./Shiny/data_GSE10810.RData)
  * [`data_Combined.RData`](./Shiny/data_Combined.RData)
  These are **filtered gene expression datasets**, pre-processed and ready for use — no need to download from external databases.

* Pre-trained model files (`.rds`), generated via [`save_models.R`](./Shiny/save_models.R):

  * [`GSE15852_models.rds`](./Shiny/GSE15852_models.rds)
  * [`GSE17907_models.rds`](./Shiny/GSE17907_models.rds)
  * [`Combined_models.rds`](./Shiny/Combined_models.rds)
  * [`comparison_cache.rds`](./Shiny/comparison_cache.rds)
  These files store the results of three predictive models trained on the datasets. They are loaded directly by the app.

* Optional media:

  * [`Record_Demo.mp4`](./Shiny/Record_Demo.mp4) — A recorded demo of the app in use (local video file).

### **Report/**
  
  Contains the final report and the R Markdown code to reproduce all plots included in the report. Key files:
  - `report_plot.RData` — saved plots for report
  - report scripts and documentation

### **Pre-Project/**
  
  Contains our earlier work, including preliminary data exploration and initial Shiny app drafts. Useful for understanding our iterative workflow.

## Datasets

We used three main microarray datasets:

- **Dataset 1**: **GSE15852** — 86 samples, 43 tumors
- **Dataset 2**: **GSE10810** — 58 samples, 31 tumors
- **Dataset 3**: **GSE17907** — 55 samples, 47 tumors
