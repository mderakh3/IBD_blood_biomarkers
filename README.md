# Machine Learning Based Diagnosis of Inflammatory Bowel Disease Using Liquid Biopsies

__Table of Contents__

1. Data Preprocessing
2. Data Integration
3. Differential Gene Expression Analysis
4. Biomarker Panel Discovey
5. Real-life Cohort Panel Validation

__Step 1: Data Preprocessing__

In this step, by using the preprocessing script, the microarray expression input data was preprocessed.

__Step 2: Data Integration__

Using the batch_integration script, the datasets were integrated to create a uniform metdadata for expression analysis and batch effects removal process was followed. Additionally, genes that did not have significant signals in microarray were removed to ensure the capturability of biomarkers in experimental procedures. 

__Step 3: Differential Gene Expression Analysis__

Bulk transcriptomic profiles of cases and controls were compared using the 
