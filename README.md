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

If you have an Nvidia GPU, then you must install CUDA and cuDNN libraries. See:
https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html
https://developer.nvidia.com/cudnn
Be aware that you should check the compatibility between your graphic card and the versions of CUDA and cuDNN you want to install. This is a bit tricky and time consuming!

If you don't have an Nvidia GPU, you can use the CPU, which will be slower for deep learning computations. 

You need to install tensorflow and keras R packages. See the following link for explanations:
https://tensorflow.rstudio.com/installation/gpu/local_gpu/
Ideally, you should install libraries for GPUs if you have an Nvidia GPU. If not, you can also install tensorflow and keras for cpu. 

The following R libraries must be installed (using install.packages() function):
- pROC
- reticulate
- tensorflow
- keras
- reticulate
- glmnet
- doParallel
- ranger
- e1071
- gkmSVM
- data.table
- ggplot2
- gplots

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

<!-- SESSION 1 -->
## Day 1

During the morning of the first day, we will first preprocess ChIP-seq/DNase-seq (positive) peaks obtained from different experiments: CTCF, POL2, H3K4me3 and DNase-seq. From the peaks, we will generate random control (negative) peaks with similar length, GC content and repeat distribution. From the peaks, we will extract DNA sequences. Then, we will split the sequences in a train and a test sets for machine/deep learning model training and prediction evaluation. From the sequences, we will build features using k-mer counts and known DNA binding protein motif counts. These features will be used for building machine learning models such as logistic lasso regression, random forests and support vector machines. 




<!-- SESSION 2 -->
## Day 2
