

# conda activate /home/mourad/.local/share/r-miniconda/envs/r-reticulate
# install_tensorflow(method = 'conda', envname = '/home/raphael/.local/share/r-miniconda/envs/r-reticulate')


#### LOAD LIBRARIES

library(GenomicRanges)
library(gkmSVM)
library(Biostrings)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(kebabs)
library(pROC)
library(glmnet)
library(doParallel)
library(ranger)
library(e1071)


#### SETUP PROJECT FOLDER

#setwd("/media/mourad/diskSave/MCF_Toulouse/enseignement/atelier_INSERM/")
setwd("/media/raphael/SSD2/atelier_INSERM")


#### LOAD FUNCTIONS

source("scriptR/functions.R")

use_condaenv('r-reticulate')
tf$constant("Hello Tensorflow")


#### SOME PARAMETERS

peakSize=201
kpeaks=4000
expe="DNase"
DNAletters=c("A","C","G","T")


#### INPUT FILE NAMES

fileBedPos=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_pos.bed")
fileBedNeg=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_neg.bed")
fileFastaPos=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_pos.fa")
fileFastaNeg=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_neg.fa")


#### LOAD PROCESSED DATA

# Load ChIP-seq and control fasta sequences
peakPos.seq=readDNAStringSet(fileFastaPos)
peakNeg.seq=readDNAStringSet(fileFastaNeg)

# Bind positive and negative sequences, and make label
peakAll.seq=c(peakPos.seq,peakNeg.seq)
label=c(rep(1,length(peakPos.seq)),rep(0,length(peakNeg.seq)))

# Shuffle sequence indices
idxS=sample(1:length(peakAll.seq))
peakAllS.seq=peakAll.seq[idxS]
labelS=label[idxS]

# Split train and test indices
percTrain=0.7
idxTrain=1:(ceiling(length(labelS)*percTrain))
idxTest=(length(idxTrain)+1):length(labelS)
labelTrain=labelS[idxTrain]
labelTest=labelS[idxTest]


#### BUILD FEATURES FOR MACHINE LEARNING

# Parallel computing (for glmnet)
registerDoParallel(4)

# K-mer counts as features
specK=gappyPairKernel(k=3, m=3, normalized=F) # normalization accounts for different sequence length
kmerAllS=getExRep(peakAllS.seq,kernel=specK,sparse=T)
kmerAllS=as(kmerAllS,"Matrix")
kmerTrain=kmerAllS[idxTrain,]
kmerTest=kmerAllS[idxTest,]

# PWM motif counts as features
# Load DNA binding protein motif PWMs from JASPAR database
opts <- list(species=9606, all_versions=TRUE)
PFMatrixList <- getMatrixSet(JASPAR2020, opts)
proteinNames=name(PFMatrixList)
motifID=names(PFMatrixList)

# DNA motif counts
motif_ix=matchMotifs(PFMatrixList,peakAllS.seq,out="scores",p.cutoff=1e-4)
mcAllS=motifCounts(motif_ix)
mcTrain=mcAllS[idxTrain,]
mcTest=mcAllS[idxTest,]



#### MACHINE LEARNING METHODS

### USING KMERS

# Lasso logistic regression
# NB: glmnet is optimized for sparse X input matrix.
cvlasso_kmer=cv.glmnet(x=kmerTrain,y=labelTrain,family="binomial",parallel=TRUE)
predLasso_kmer=predict(cvlasso_kmer,newx=kmerTest,type="response")[,1]
rocLasso_kmer=pROC::roc(as.factor(labelTest),predLasso_kmer,ci=T)
plot(rocLasso_kmer,main=paste0("AUROC: ", round(pROC::auc(rocLasso_kmer),3)))


### USING DNA MOTIFS

# Lasso logistic regression
cvlasso_motif=cv.glmnet(x=mcTrain,y=labelTrain,family="binomial",parallel=TRUE)
predLasso_motif=predict(cvlasso_motif,newx=mcTest,type="response")[,1]
rocLasso_motif=pROC::roc(as.factor(labelTest),predLasso_motif,ci=T)
plot(rocLasso_motif,main=paste0("AUROC: ", round(pROC::auc(rocLasso_motif),3)))

# Random forests
dataRF_motif=data.frame(label=labelTrain,as(mcTrain,"matrix"))
RF_motif=ranger(label ~ .,data=dataRF_motif,importance="permutation")
predRF_motif=predict(RF_motif,data=data.frame(as(mcTest,"matrix")))$predictions
rocRF_motif=pROC::roc(as.factor(labelTest),predRF_motif,ci=T)
plot(rocRF_motif,main=paste0("AUROC: ", round(pROC::auc(rocRF_motif),3)))

# Variable importance
dataImportanceMotif=data.frame(motifID,proteinNames,nameID=paste0(motifID,"_",proteinNames),importance=importance(RF_motif))
dataImportanceMotif=dataImportanceMotif[order(importance(RF_motif),decreasing=T)[1:20],]
dataImportanceMotif=dataImportanceMotif[order(dataImportanceMotif[,4]),]
par(mar=rep(10,4))
barplot(dataImportanceMotif[,4],names.arg=dataImportanceMotif[,3],horiz=T,las=2,xlab="Importance")

# Support vector machine
SVM_motif=svm(label ~ ., data=dataRF_motif, kernel="linear",scale=F)
predSVM_motif=predict(SVM_motif,newdata=data.frame(as(mcTest,"matrix")))
rocSVM_motif=pROC::roc(as.factor(labelTest),predSVM_motif,ci=T)
plot(rocSVM_motif,main=paste0("AUROC: ", round(pROC::auc(rocSVM_motif),3)))




