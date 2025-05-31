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

- **Shiny/**
  
  Contains the complete Shiny application, including:
  - `app.R` — the main app code
  - `www/` — images used by the Shiny app
  - data files and supporting materials

- **Report/**
  
  Contains the final report and the R Markdown code to reproduce all plots included in the report. Key files:
  - `report_plot.RData` — saved plots for report
  - report scripts and documentation

- **Pre-Project/**
  
  Contains our earlier work, including preliminary data exploration and initial Shiny app drafts. Useful for understanding our iterative workflow.

## Datasets

We used three main microarray datasets:

- **Dataset 1**: **GSE15852** — 86 samples, 43 tumors
- **Dataset 2**: **GSE10810** — 58 samples, 31 tumors
- **Dataset 3**: **GSE17907** — 55 samples, 47 tumors
