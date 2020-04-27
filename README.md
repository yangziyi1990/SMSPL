# SMSPL: Robust multimodal approach to integrative analysis of multi-omics data

This repository provides some codes for this paper. 

The code Version 1.0.

If you find it interesting and confused about the content, please contact me.

communication E-mail: yangziyi091100@163.com

## I . Abstract

With recent advancement of technologies, it is progressively easier to produce diverse types of genome-wide data. Integrative analysis is a commonly used method to integrate these data and is expected to improve our understanding of a complex biological system. Current methods, however, have some limitations and are prone to overfitting heavy noise. High noise is one of the major challenges for multi-omics data integration, which may be susceptible to overfitting and lead to poor performance in generalization. Sample reweighting strategy is typically used to address this problem. In this paper, we propose a robust multimodal data integration method, termed as SMSPL, which can simultaneously predict subtypes of cancer and identify potentially significant multi-omics signatures. Especially, the proposed method leverages the linkages between different types of data to interactively recommend high-confidence samples, adopts a new soft weighting scheme to assign weights to the training samples of each type, and then iterates between weights recalculating and classifiers updating. Simulation and five real experiments substantiate the capability of the proposed method for classification and identification of significant multi-omics signatures in heavy noise. We expect SMSPL to take a small step in the multi-omics data integration and help researchers comprehensively understand the biological process.

If you find this code useful in your research then please cite:
```bash
@article{
  title={SMSPL: Robust Multimodal Approach to Integrative Analysis of Multi-omics Data},
  author={Yang, Zi-Yi, Wu Nai-Qi, Liang Yong, Zhang Hui, and Ren Yan-Qiong},
  journal={IEEE Transactions on Cybernetics},
  year={2020},
  publisher={IEEE}
}
```

## II. Introduce about code

### i . The repository can be divided into three parts: 

1. Codes for simulated experiments.
2. Codes for benchmark multi-omics cancer datasets. (binary classification problem)
3. Codes for breast multi-omics cancer datasets. (multiple classification problem)

### ii . The compared methods:

Current supervised multimodal data integration approaches for predicting cancer subtypes and identifying significant multi-omics signatures can be classified as concatenation-based, ensemble-based, and knowledge-driven approaches.

We applied the logistic regression model/multinomial model with Elastic Net (EN) regularization [1], Random Forest (RF) [2] and Self-paced Learning (SPL) with L1 penalty [3] in the concatenation and ensemble frameworks. The compared methods include: concatenation-based methods (Concate\_EN, Concate\_ RF, and Concate\_SPL), ensemble-based methods (Ensemble\_EN, Ensemble\_RF, and Ensemble\_SPL) and DIABLO [4]. In this paper, we proposed a novel model for multi-omics data integration, termed as SMSPL. 

### iii. The SMSPL model:

The objective function of SMSPL can be expressed as:
$$
\min_{\substack{\beta^{\left(j\right)},v^{\left(j\right)}\in\left[0,1\right],\\j=1,2,\ldots,m}} E\left(\beta^{(j)}, v^{(j)} ; \lambda^{(j)}, \gamma^{(j)}, \delta\right)=  \sum_{j=1}^{m} \sum_{i=1}^{n} v_{i}^{(j)} L\left(y_{i}, f\left(x_{i}^{(j)}, \beta^{(j)}\right)\right) +\sum_{j=1}^{m}\lambda^{(j)} \left\|\beta^{(j)}\right\|_{1} \\
-\sum_{j=1}^{m} \sum_{i=1}^{n} \gamma^{(j)} v_{i}^{(j)} + \frac{\delta}{2(m-1)} \sum_{j=1}^{m} \sum_{k\neq j}^{m}\left\|v^{(j)}-v^{(k)}\right\|_{2}^{2}
$$

## III. Codes for simulated experiments

The codes for simulated experiments contain three parts:

1. Generating simulated data.
2. Train different data sets using training data sets and get the best solution for each model.
3. Calculate the prediction and feature selection performance. (We applied five indicators to evaluate the prediction performance: accuracy, sensitivity, specificity, recall, and AUC. To evaluate the feature selection performance, $\beta$-sensitivity and $\beta$-specificity are applied.)

Running "*\_simu.R" to perform the simulated experiments. Before running the code, it is need to set the path "setwd("D:/1-simulation")". After that, running the code. For example:

```
Rscript 8-SMSPL-simu.R
```

The defined of $\beta$-sensitivity and $\beta$-specificity are as follows:
$$
True Positive (TP)= \left|\beta.\ast\hat{\beta}\right|_0,   True Negative (TN)= \left|\bar{\beta}.\ast\bar{\hat{\beta}}\right|_0 \\
False Positive (FP)= \left|\bar{\beta}.\ast\hat{\beta}\right|_0 ,   False Negative (FN)=\ \left|\beta.\ast\bar{\hat{\beta}}\right|_0 \\
\beta-sensitivity= \frac{TP}{TP+FN}\\
\beta-specificity=\ \frac{TN}{TN+FP}
$$
where the $|\cdot|_0$ represents the number of non-zero elements in a vector. The logical non operators of  $\beta$ and $\hat{\beta}$ are $\bar{\beta}$ and $\bar{\hat{\beta}}$, respectively. And $.\ast$ is the element-wise product.

***Special comments:***

Since biological samples are complex and cannot be visualized, the choice of model parameters is a challenge. In order to make the model better select the samples in each class in the process of sample selection, we added the parameters of the selected sample size in the process of self-step learning.

## IV. Codes for benchmark multi-omics cancer data sets

Four benchmark multi-omics cancer datasets (mRNA, miRNA and DNA methylation) were obtained from [5]: Glioblastoma multi-forme (GBM), Kidney renal clear cell carcinoma (KRCCC), Lung squamous cell carcinoma (LSCC), Colon adenocarcinoma (COAD). Survival times were provided for each disease cohort by [5]. By using the median survival time, we dichotomized the samples into two classes in low and high survival times. A brief description of these four benchmark datasets is summarized in Table I.

Table I. The measurements of sample sizes and the number of features in each omics for four benchmark cancer datasets.

| Benchmark datasets | Samples (high/low) | mRNA  | miRNA | methylation |
| :----------------: | :----------------: | :---: | :---: | :---------: |
|       KRCCC        |    122 (61/61)     | 17665 |  329  |    24960    |
|        LSCC        |    106 (53/53)     | 12042 |  352  |    23074    |
|        GBM         |   213 (105/108)    | 12042 |  534  |    1305     |
|        COAD        |     92 (33/59)     | 17814 |  312  |    23088    |



The codes for simulated experiments contain three parts:

1. Pre-processing benchmark cancer datasets.
2. Train different data sets using training data sets and get the best solution for each model.
3. Calculate the prediction performance. (We applied five indicators to evaluate the prediction performance: accuracy, sensitivity, specificity, recall, and AUC.)

Running "*-bench.R" to perform the benchmark experiments. Before running the code, it is need to set the path "setwd("D:/2-benchmark")". Moreover, some parameters need to adjust according to the specific problems. After that, running the code. For example:

```
Rscript 9-SMSPL-bench.R
```

Intermediate result presentation:

```
Starting the 1-th iteration.
The 3-th modality select 12 samples.
Starting the 2-th iteration.
The 3-th modality select 16 samples.
Starting the 3-th iteration.
The 3-th modality select 20 samples.
Starting the 4-th iteration.
The 3-th modality select 24 samples.
Starting the 5-th iteration.
The 3-th modality select 28 samples.
Starting the 6-th iteration.
The 3-th modality select 32 samples.
Starting the 7-th iteration.
The 3-th modality select 36 samples.
Starting the 8-th iteration.
The 3-th modality select 40 samples.
...
```

***Special comments:***

For the benchmark cancer data experiment, we evaluate the classifier for the binary classification problem.

Since biological samples are complex and cannot be visualized, the choice of model parameters is a challenge. In order to make the model better select the samples in each class in the process of sample selection, we added the parameters of the selected sample size in the process of self-step learning.

## V. Codes for breast multi-omics cancer data set

We curated breast cancer multi-omics dataset (mRNA, miRNA and methylation) from the Cancer Genome Atlas (TCGA, data version 2015 11 01 for BRCA) in order to achieve a systems characterization of breast cancer subtypes with multiple omics. 

The raw data of breast cancer multi-omics can be download from the website: http://gdac.broadinstitute.org/runs/stddata__2015_11_01/data/BRCA/20151101/.

This dataset contains four subtypes of breast cancer: Luminal A (LumA), Luminal B (LumB), Her2-enriched (Her2) and Basal-like (Basal), which have been reported the most replicated subtypes of human breast cancer. The miRNA dataset was derived from two different Illumina technologies: The Illumina Genome Analyzer and the Illumina Hiseq. The methylation data was derived from two different platforms: the Illumina Methylation 27 and the Illumina 450K.

The codes for simulated experiments contain three parts:

1. Pre-processing breast cancer multi-omics dataset.
2. Train different data sets using training data sets and get the best solution for each model.
3. Calculate the prediction performance. (We applied five indicators to evaluate the prediction performance: accuracy, sensitivity, specificity, recall, and AUC.)

Running "*\_brac.R" to perform the benchmark experiments. Before running the code, it is need to set the path "setwd("D:/3-breast cancer")". Moreover, some parameters need to adjust according to the specific problems. After that, running the code. For example:

```
Rscript 9-SMSPL_brac.R
```

A brief presentation of the results:

```
> perf.SMSPL
$coef.mRNA
       [,1]                [,2]                   
  [1,] "mrna_ABCC11"       "-0.0218210423558085"  
  [2,] "mrna_AGR3"         "-0.0568151686777767"  
  [3,] "mrna_AKR1A1"       "-0.0259073635870535"  
  [4,] "mrna_AKR7A3"       "-0.00294072007786631" 
  [5,] "mrna_ALOX15B"      "-0.00137019885619369" 
  [6,] "mrna_AMY1A"        "0.0119740585542367"   
  [7,] "mrna_ANO4"         "0.00043343797513604"  
  [8,] "mrna_ANXA8L2"      "0.0157983044611121"   
  [9,] "mrna_ARHGAP11A"    "0.00578362102698803"  
 [10,] "mrna_ARSG"         "-0.0728438250253021" 
 ...
 $coef.miRNA
       [,1]                   [,2]                   
  [1,] "mirna_hsa-let-7b"     "0.00396065817585336"  
  [2,] "mirna_hsa-let-7c"     "0.162370488567795"    
  [3,] "mirna_hsa-mir-1-2"    "0.00242943836573586"  
  [4,] "mirna_hsa-mir-100"    "0.013974246776146"    
  [5,] "mirna_hsa-mir-101-1"  "-0.0957031286884165"  
  [6,] "mirna_hsa-mir-10b"    "-0.00292397943317286" 
  [7,] "mirna_hsa-mir-1245"   "-0.00499733706999291" 
  [8,] "mirna_hsa-mir-1247"   "-0.0273164057148487"  
  [9,] "mirna_hsa-mir-125a"   "0.0290930959106511"   
 [10,] "mirna_hsa-mir-1266"   "-0.0658790022208466"  
 ...
 $coef.cpg
       [,1]                       [,2]                   
  [1,] "meth_A2ML1"               "-1.10270381663695"    
  [2,] "meth_ABCG8;ABCG5"         "0.97328628147582"     
  [3,] "meth_ABTB1"               "-0.000825076734604522"
  [4,] "meth_ACAP3"               "0.0215346116519729"   
  [5,] "meth_ACCN4"               "-0.00719535131507551" 
  [6,] "meth_ADIPOQ"              "2.16978529867462"     
  [7,] "meth_ADRBK2"              "0.0294242144417417"   
  [8,] "meth_AMBP"                "0.00545885796362712"  
  [9,] "meth_AMOT"                "-0.364242564706836"   
 [10,] "meth_AMPD1"               "0.0173046367250456" 
 ...
 $feature.num
     [,1] [,2] [,3]
[1,]  171  108  208

$Perf.Train
Confusion Matrix and Statistics

          Reference
Prediction   1   2   3   4
         1 102   0   0   0
         2   0  40   0   0
         3   0   0 342   4
         4   0   0   1 121

Overall Statistics
                                         
               Accuracy : 0.9918         
                 95% CI : (0.981, 0.9973)
    No Information Rate : 0.5623         
    P-Value [Acc > NIR] : < 2.2e-16      
                                         
                  Kappa : 0.9865         
 Mcnemar's Test P-Value : NA             

Statistics by Class:

                     Class: 1 Class: 2 Class: 3 Class: 4
Precision              1.0000  1.00000   0.9884   0.9918
Recall                 1.0000  1.00000   0.9971   0.9680
F1                     1.0000  1.00000   0.9927   0.9798
Prevalence             0.1672  0.06557   0.5623   0.2049
Detection Rate         0.1672  0.06557   0.5607   0.1984
Detection Prevalence   0.1672  0.06557   0.5672   0.2000
Balanced Accuracy      1.0000  1.00000   0.9911   0.9830

$Perf.Test
Confusion Matrix and Statistics

          Reference
Prediction   1   2   3   4
         1  75   1   0   0
         2   1  27   4   6
         3   0   0 176  12
         4   0   3  12  62

Overall Statistics
                                         
               Accuracy : 0.8971         
                 95% CI : (0.862, 0.9258)
    No Information Rate : 0.5066         
    P-Value [Acc > NIR] : < 2.2e-16      
                                         
                  Kappa : 0.8435         
 Mcnemar's Test P-Value : NA             

Statistics by Class:

                     Class: 1 Class: 2 Class: 3 Class: 4
Precision              0.9868  0.71053   0.9362   0.8052
Recall                 0.9868  0.87097   0.9167   0.7750
F1                     0.9868  0.78261   0.9263   0.7898
Prevalence             0.2005  0.08179   0.5066   0.2111
Detection Rate         0.1979  0.07124   0.4644   0.1636
Detection Prevalence   0.2005  0.10026   0.4960   0.2032
Balanced Accuracy      0.9918  0.91968   0.9262   0.8624

```

***Special comments:***

For the breast cancer multi-omics data experiment, we evaluate the classifier for the multiple classification problem.

Since biological samples are complex and cannot be visualized, the choice of model parameters is a challenge. In order to make the model better select the samples in each class in the process of sample selection, we added the parameters of the selected sample size in the process of self-step learning.



## VI. Reference:

[1] H. Zou and T. Hastie, “Regularization and variable selection via the elastic net,” Journal of the royal statistical society:  series B (statistical methodology), vol. 67, no. 2, pp. 301–320, 2005.

[2] A. Liaw, M. Wiener et al., “Classification and regression by randomforest,” R news, vol. 2, no. 3, pp. 18–22, 2002.

[3] M. P. Kumar, B. Packer, and D. Koller, “Self-paced learning for latent variable models,” in Advances in Neural Information Processing Systems, 2010, pp. 1189–1197.

[4] A. Singh, C. P. Shannon, B. Gautier, F. Rohart, M. Vacher, S. J. Tebbutt, and K.-A. Lˆe Cao, “Diablo: an integrative approach for identifying key molecular drivers from multi-omics assays,” Bioinformatics, 2019.

[5] B. Wang, A. M. Mezlini, F. Demir, M. Fiume, Z. Tu, M. Brudno, B. Haibe-Kains, and A. Goldenberg, “Similarity network fusion for aggregating data types on a genomic scale,” Nature methods, vol. 11, no. 3, p. 333, 2014.









