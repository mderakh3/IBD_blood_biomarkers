# Machine Learning Based Diagnosis of Inflammatory Bowel Disease Using Liquid Biopsies

__Table of Contents__

1. Data Preprocessing
2. Data Integration
3. Differential Gene Expression Analysis
4. Functional Annotation
5. Biomarker Panel Discovey
6. Real-life Cohort Panel Validation

__Step 1: Data Preprocessing__

In this step, by using the preprocessing script, the microarray expression input data was preprocessed.

__Step 2: Data Integration__

Using the batch_integration script, the datasets were integrated to create a uniform metdadata for expression analysis and batch effects removal process was followed. Additionally, genes that did not have significant signals in microarray were removed to ensure the capturability of biomarkers in experimental procedures. 

__Step 3: Differential Gene Expression Analysis__

Bulk transcriptomic profiles of cases and controls were compared using the differential_expression_analysis script. 

__Step 4: Functional Annotation__

By using the functional_annotation script, functional enrichment analysis and network analysis were followed. 

__Step 5: Biomarker Panel Discovery__

Using the ... script, IBD-specific diagnsotic were identified and 20/80 split classification was performed using support vector machine algorithm. 

__Step 6: Real-life Cohort Panel Validation__

Using the ... script, the developed diagnositc biomarker panel was evaluated to classifiy patients in the real-life cohort using qRT-PCR expression data. 
