# Phase 2 : Inserm workshop on Machine Learning in Biology

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-workshop">About The Workshop</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#session1">Session 1</a></li>
    <li><a href="#session2">Session 2</a></li>
  </ol>
</details>


<!-- ABOUT THE WORKSHOP -->
## About The Workshop

The aim of this practical session is to learn machine and deep learning models for regulatory genomics using R. Deep learning models will be programmed with the Keras R interface of Tensorflow. Most models we will use can be used for many different problems in biology and in science in general. 
We will focus on the prediction of genomic experimental output such as transcription factor ChIP-seq or DNase-seq, using the DNA sequence as input. Once the model is built and trained, we will see how to extract and interpret the model features as DNA motifs. Lastly, we will use our model for predicting the effect of single nucleotide polymorphisms (SNPs) that were previously identified from genetic studies (GWAS, eQTL). 

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

The following R libraries must be installed:
- GenomicRanges
- Biostrings
- pROC
- tensorflow
- keras
- reticulate
- TFBSTools
- JASPAR2020
- tensorflow
- keras
- reticulate
- motifStack
- universalmotif
- motifmatchr
- kebabs
- pROC
- glmnet
- doParallel
- ranger
- e1071
- gkmSVM
library(data.table)
library(ggplot2)
library(gplots)

The following R libraries must be installed (using bioconductor):
- BSgenome.Hsapiens.UCSC.hg19
- BSgenome.Hsapiens.UCSC.hg19.masked
- GenomicRanges
- Biostrings
- TFBSTools
- JASPAR2020
- motifStack
- universalmotif
- motifmatchr
- kebabs


### Installation


<!-- SESSION 1 -->
## Session 1

<!-- SESSION 2 -->
## Session 2
