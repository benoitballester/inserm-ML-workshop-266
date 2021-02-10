

# conda activate /home/mourad/.local/share/r-miniconda/envs/r-reticulate
# install_tensorflow(method = 'conda', envname = '/home/raphael/.local/share/r-miniconda/envs/r-reticulate')


#### LOAD LIBRARIES

library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(gkmSVM)
library(Biostrings)
library(JASPAR2018)
library(TFBSTools)
library(motifmatchr)
library(kebabs)
library(pROC)
library(glmnet)
library(doParallel)
library(ranger)
library(e1071)
library(tensorflow)
library(keras)
library(reticulate)
library(ggseqlogo)
library(motifStack)
library(data.table)


#### LOAD FUNCTIONS

source("scriptR/functions.R")

use_condaenv('r-reticulate')
tf$constant("Hello Tensorflow")


#### SETUP PROJECT FOLDER

#setwd("/media/mourad/diskSave/MCF_Toulouse/enseignement/atelier_INSERM/")
setwd("/media/raphael/SSD2/atelier_INSERM")

#### CREATE FOLDERS

dir.create("data")
dir.create("data/bed")
dir.create("data/fasta")
dir.create("data")
dir.create("results")
dir.create("results/model")


#### SOME PARAMETERS

peakSize=201
kpeaks=4000
expe="POL2"
DNAletters=c("A","C","G","T")


#### INPUT FILE NAMES

fileBedPos=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_pos.bed")
fileBedNeg=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_neg.bed")
fileFastaPos=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_pos.fa")
fileFastaNeg=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_neg.fa")



#### LOAD, PROCESS AND SAVE PROCESSED DATA

## INFO: TO RUN ONLY ONCE!

# Chromosome information for hg19
Genome=BSgenome.Hsapiens.UCSC.hg19
Chr.V=c(paste0("chr",1:22),"chrX")
SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)[Chr.V]

# Load ChIP-seq/DNase-seq peaks
file_peaks=paste0("data/bed/",expe,"_GM12878_hg19.narrowPeak")
dataPeaks=read.table(file_peaks,sep="\t",header=F)

# Resize peaks
peaks.GR=GRanges(dataPeaks[,1],IRanges(dataPeaks[,2],dataPeaks[,3]),score=dataPeaks[,7])
peaks.GR=resize(peaks.GR, width=peakSize, fix="center")

# Select k best peaks and export to bed
peaksSorted.GR=peaks.GR[order(peaks.GR$score,decreasing=T)]
dataPeaksSorted=as.data.frame(peaksSorted.GR)
dataPeaksBest=dataPeaksSorted[1:kpeaks,]
write.table(dataPeaksBest[,1:3],fileBedPos,sep='\t',col.names=F,row.names=F,quote=F)

# Generate fasta sequences from the experimental peaks
# as well as control fasta sequences
# GC content and repeat distribution similar to the experimental peaks.
# The number of control peaks drawn should be similar to experimental peaks (balanced dataset).
genNullSeqs(inputBedFN=fileBedPos,nMaxTrials=10,xfold=1.5,genomeVersion="hg19",
	outputPosFastaFN=fileFastaPos,outputBedFN=fileBedNeg,outputNegFastaFN=fileFastaNeg,length_match_tol=0)

# genNullSeqs generates sequences of length 200b instead of 201b
# We want to correct this.
dataPeaksNeg=read.table(fileBedNeg,sep="\t",header=F)
peaksNeg.GR=GRanges(dataPeaksNeg[,1],IRanges(dataPeaksNeg[,2],dataPeaksNeg[,3]+1))
write.table(as.data.frame(peaksNeg.GR)[,1:3],fileBedNeg,sep='\t',col.names=F,row.names=F,quote=F)
peaksNeg.seq=getSeq(BSgenome.Hsapiens.UCSC.hg19, names=seqnames(peaksNeg.GR),
     start=start(peaksNeg.GR), end=end(peaksNeg.GR))
writeXStringSet(peaksNeg.seq,fileFastaNeg)
peaksPos.GR=GRanges(dataPeaksBest[,1],IRanges(dataPeaksBest[,2],dataPeaksBest[,3]))
peaksPos.seq=getSeq(BSgenome.Hsapiens.UCSC.hg19, names=seqnames(peaksPos.GR),
     start=start(peaksPos.GR), end=end(peaksPos.GR))
writeXStringSet(peaksPos.seq,fileFastaPos)


#### LOAD PROCESSED DATA

# Load ChIP-seq and control fasta sequences
peakPos.seq=readDNAStringSet(fileFastaPos)
peakNeg.seq=readDNAStringSet(fileFastaNeg)

# Remove control sequences on poorly assembled chromosomes (for instance chr17_gl000204)
peakNeg.seq=peakNeg.seq[-grep("gl",names(peakNeg.seq))]

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
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
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




#### DEEP LEARNING METHODS


# One hot encoding
oneHotTrain=convertOneHot(peakAllS.seq[idxTrain,])
oneHotTest=convertOneHot(peakAllS.seq[idxTest,])

# Vocabulary size
vocab_size=4 # A, T, G, C

# CNN 1 (simple convolution + 10 dense). AUC=0.9
# Simple yet very efficient model.
kernelSize=16
model <- keras_model_sequential()
model %>% 
  layer_conv_1d(filters = 256, kernel_size = kernelSize, activation = 'relu',
                input_shape = c(peakSize,vocab_size), name="conv1d_cnn1") %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 10, activation = "relu", name="10dense_cnn1") %>%
  layer_dense(units = 1, activation = 'linear', name="2ndlastlayer_cnn1") %>%
  layer_dense(units = 1, activation = 'sigmoid', name="lastlayer_cnn1") 

# CNN 2 (parallel convolution). AUC=0.9
# Useful to model different DNA motif sizes.
if(F){
inputs <- layer_input(shape = c(peakSize-1,vocab_size)) 
conv1 <- inputs %>%
  layer_conv_1d(filters = 256, kernel_size = 8, activation = 'relu') %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0.5) 
conv2 <- inputs %>%
  layer_conv_1d(filters = 256, kernel_size = 16, activation = 'relu') %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0.5) 
conv3 <- inputs %>%
  layer_conv_1d(filters = 256, kernel_size = 24, activation = 'relu') %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0.5) 
outputs <- layer_concatenate(c(conv1,conv2,conv3)) %>%
  layer_dense(units = 10, activation = "relu") %>%
  layer_dense(units = 1, activation = 'sigmoid') 
model <- keras_model(inputs = inputs, outputs = outputs)
}

# CNN 3 (double convolution + 10 dense). AUC=0.89
# Useful for computer vision, but not for DNA sequences.
if(F){
model <- keras_model_sequential()
model %>% 
  layer_conv_1d(filters = 256, kernel_size = 8, activation = 'relu',
                input_shape = c(peakSize-1,vocab_size)) %>% 
  layer_conv_1d(filters = 256, kernel_size = 3, activation = 'relu',
                input_shape = c(peakSize-1,vocab_size)) %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 10, activation = "relu") %>%
  layer_dense(units = 1, activation = 'sigmoid') 
}

# RNN with LSTM (no CNN). AUC=0.6
# Bad results. We need to include a convolution layer before the lstm layer!
# Slow (not made for efficient parallel computing)
if(F){
model <- keras_model_sequential()
model %>% 
  layer_lstm(128,input_shape = c(peakSize-1,vocab_size)) %>%
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 10, activation = "relu") %>%
  layer_dense(1) %>%
  layer_activation("sigmoid")
}

# CNN + LTSM. AUC= . (Less good results)
# Simple yet very efficient model.
if(F){
model <- keras_model_sequential()
model %>% 
  layer_conv_1d(filters = 128, kernel_size = 24, activation = 'relu',
                input_shape = c(peakSize-1,vocab_size)) %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_lstm(12) %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 20, activation = "relu") %>%
  layer_dense(units = 1, activation = 'sigmoid') 
}

# Model summary
model %>% summary()

# Compile the model
# SGD and Adagrad = bad results
# Adam and RMSprop = good results
model %>% compile(
  optimizer = 'rmsprop', # adam or rmsprop are great here
  loss = 'binary_crossentropy',
  metrics = list('accuracy')
)

# Train the model
# For CNN1, epochs=8 (if more, overfitting)
# For CNN2, epochs=10 (if more, overfitting)
# For CNN3, epochs=10 (if more, overfitting)
# For CNN+LSTM, epochs=25
history <- model %>% fit(
  oneHotTrain,
  labelTrain,
  epochs = 8, 
  batch_size = 128,
  validation_split = 0.2,
  verbose=1
)
plot(history)

# Prediction
predCNN=predict(model,oneHotTest)
rocCNN_motif=pROC::roc(as.factor(labelTest),predCNN,ci=T)
plot(rocCNN_motif,main=paste0("AUROC: ", round(pROC::auc(rocCNN_motif),3)))

# Accuracy and loss
acc_loss <- model %>% evaluate(oneHotTest, labelTest)
acc_loss

# Save model
file_model=paste0("results/model/CNN1_model_",expe,".hdf5")
save_model_hdf5(model,file_model)


#### FEATURE EXTRACTION AND INTERPRETATION

# Extract kernels from CNN1 convolutional layer
# object "model" should be the CNN1!
conv1d_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, 'conv1d_cnn1')$output)
conv1d_output <- predict(conv1d_layer_model, oneHotTrain)
print(dim(conv1d_output)) # To see the size of the tensor

# Convert kernels to PFM matrices
activationThreshold=0.15
conv1d_posmax=apply(conv1d_output,c(1,3),function(x){if(max(x)>activationThreshold){which.max(x)}else{NA}})
print(dim(conv1d_posmax))
peakTrain.seq=peakAllS.seq[idxTrain]
kernelPFMList=list()
kernelActivatedList.seq=list()
for(i in 1:ncol(conv1d_posmax)){
 if(sum(!is.na(conv1d_posmax[,i]))>0){
  idxSeqActivatedi=which(!is.na(conv1d_posmax[,i]))
  seqActivatedi.seq=peakTrain.seq[idxSeqActivatedi]
  seqActivatedi.char=as.character(seqActivatedi.seq)
  idxStarti=conv1d_posmax[idxSeqActivatedi,i]
  idxEndi=conv1d_posmax[idxSeqActivatedi,i]+kernelSize-1
  subSeqActivatedi.char=as.character(sapply(1:length(seqActivatedi.char),function(x){substring(seqActivatedi.char[x],idxStarti[x],idxEndi[x])}))
  subSeqActivatedi.mat=do.call(rbind,strsplit(subSeqActivatedi.char,""))
  #subSeqActivatedi.df=do.call(cbind,lapply(1:ncol(subSeqActivatedi.mat),function(x){factor(subSeqActivatedi.mat[,x],levels=c("A","C","G","T"))}))
  PFMmati=apply(subSeqActivatedi.mat,2,function(x){table(factor(x,levels=c("A","C","G","T")))})
  kernelPFMList[[i]]=PFMatrix(ID=paste0("kernel",i),profileMatrix=PFMmati)
  kernelActivatedList.seq[[i]]=as.character(subSeqActivatedi.char)
  print(i)
 }
}

# Plot motif logo from PFM matrix
motif1 <- new("pcm", mat=as.matrix(kernelPFMList[[1]]), name="kernel1")
plot(motif1)

# Plot multiple motif logos from subsequences
# Plot from sequences
ggseqlogo(kernelActivatedList.seq[1:20])

# Mapping


#### DEEP LEARNING PREDICTIONS OF SNP EFFECTS

# Load CNN1 model
file_model=paste0("results/model/CNN1_model_",expe,".hdf5")
model=load_model_hdf5(file_model)


# Load SNP data (hg19)
dataSNP=as.data.frame(fread("data/SNP/Breast_Mammary_Tissue.v7.signif_variant_gene_pairs.txt.gz",header=T))
#dataSNP=as.data.frame(fread("data/SNP/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz",header=T))
dataSNP=dataSNP[abs(dataSNP$slope)>1,] # Select SNPs with high absolute effects
variantInfo=t(sapply(as.character(dataSNP[,1]),function(x){as.character(strsplit(x,'_')[[1]][1:4])}))
rownames(variantInfo)=rep("",nrow(variantInfo))
SNP.GR=GRanges(paste0("chr",variantInfo[,1]),IRanges(as.numeric(variantInfo[,2]),as.numeric(variantInfo[,2])),ref=variantInfo[,3],alt=variantInfo[,4],
	slope=dataSNP[,8],tss_distance=dataSNP$tss_distance, gene=dataSNP$gene_id) 
SNP.GR=SNP.GR[nchar(SNP.GR$ref)==1 & nchar(SNP.GR$alt)==1]


# GWAS catalog (hg38)
dataGWAS=as.data.frame(fread("data/SNP/gwas_catalog_v1.0-associations_e100_r2021-01-29.tsv.gz",header=T))
dataGWAS=dataGWAS[dataGWAS[,11]!=Inf & !is.na(dataGWAS[,11]),]
dataGWAS$slope=log(dataGWAS[,11])
dataGWAS=dataGWAS[dataGWAS$CHR_POS!="",]
dataGWAS=dataGWAS[substring(dataGWAS$SNPS,1,2)=="rs",]
riskAllele=sapply(dataGWAS[,7],function(x){strsplit(x,"-")[[1]][2]})
dataGWAS=dataGWAS[riskAllele%in%DNAletters,]
riskAllele=riskAllele[riskAllele%in%DNAletters]
SNP.GR=GRanges(paste0("chr",dataGWAS[,4]),IRanges(as.numeric(dataGWAS$CHR_POS),as.numeric(dataGWAS$CHR_POS)),alt=riskAllele,slope=dataGWAS$slope) 
SNP.GR=SNP.GR[!duplicated(SNP.GR)]
refAllele=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, SNP.GR))
SNP.GR$ref=refAllele
SNP.GR=SNP.GR[SNP.GR$ref!=SNP.GR$alt]


# Allele DB (hg19)
dataAlleleDB=read.csv("data/SNP/ASB.auto.v2.1.aug16.txt/ASB.auto.v2.1.aug16.txt",header=T,sep='\t')
alleleDB.GR=GRanges(dataAlleleDB[,1],IRanges(dataAlleleDB[,2],dataAlleleDB[,3]))
values(alleleDB.GR)=dataAlleleDB[,4:12]
alleleDB.GR=alleleDB.GR[alleleDB.GR$alt%in%c("A","T","G","C")]
alleleDB.GR$prot=sapply(as.character(alleleDB.GR$TF_indiv_ASB),function(x){strsplit(x,'_')[[1]][1]})
numAlleleTab=as.data.frame(values(alleleDB.GR))
nref=sapply(1:length(alleleDB.GR),function(x){numAlleleTab[x,substring(colnames(numAlleleTab),2)==as.character(alleleDB.GR$ref[x])]})
nalt=sapply(1:length(alleleDB.GR),function(x){numAlleleTab[x,substring(colnames(numAlleleTab),2)==as.character(alleleDB.GR$alt[x])]})
alleleDB.GR$nref=nref
alleleDB.GR$nalt=nalt
alleleDB.GR$ratio=nref/(nref+nalt)
#alleleDB.GR=alleleDB.GR[alleleDB.GR$ratio>0 & alleleDB.GR$ratio<1] # Filter out suspicious results (it seems that the sites are not heterozygous)
alleleDB.GR=alleleDB.GR[alleleDB.GR$nref>2 & alleleDB.GR$nalt>2]
SNP.GR=alleleDB.GR[alleleDB.GR$prot=="POL2"]


# Extract sequences surrounding the SNPs
region.GRr=resize(SNP.GR,peakSize,fix="center")
region.seq=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, names=seqnames(region.GRr),
     start=start(region.GRr), end=end(region.GRr)))
     
region.seq=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, names=seqnames(region.GRr),
     start=start(region.GRr), end=end(region.GRr)))


# Make sequences with reference alleles
SNPposrel=ceiling((peakSize+1)/2)
region.seqRef=region.seq
substring(region.seqRef,SNPposrel,SNPposrel)=as.character(region.GRr$ref)
names(region.seqRef)=paste0("seq",1:length(region.seqRef))
#writeXStringSet(DNAStringSet(region.seqRef),filepath=paste0("data/fasta/sequences_SNPref_hg19.fa"))

# Make sequences with alternative alleles
region.seqAlt=region.seq
substring(region.seqAlt,SNPposrel,SNPposrel)=as.character(region.GRr$alt)
names(region.seqAlt)=paste0("seq",1:length(region.seqAlt))
#writeXStringSet(DNAStringSet(region.seqAlt),filepath=paste0(pathfasta,"/sequences_SNPalt_hg19.fa"))

# Predict SNP effect on the predicted probability using the deep learning model
oneHotRef=convertOneHot(region.seqRef)
oneHotAlt=convertOneHot(region.seqAlt)
predProbRef=predict(model,oneHotRef)
predProbAlt=predict(model,oneHotAlt)
deltaProbSNP=predProbRef-predProbAlt

# Predict SNP effect on the second to last layer (just before logistic transformation) using the deep learning model
secondlast_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, '2ndlastlayer_cnn1')$output)
secondlast_outputRef <- predict(secondlast_layer_model, oneHotRef)
secondlast_outputAlt <- predict(secondlast_layer_model, oneHotAlt)
deltaSNP=secondlast_outputRef-secondlast_outputAlt

# 
boxplot(deltaSNP[SNP.GR$slope<1],deltaSNP[SNP.GR$slope>1],ylim=c(-0.2,0.2))
wilcox.test(deltaSNP[SNP.GR$slope<1],deltaSNP[SNP.GR$slope>1])


hist(deltaSNP)
hist(deltaProbSNP)
boxplot(deltaSNP[SNP.GR$ratio<0.5],deltaSNP[SNP.GR$ratio>0.5],ylim=c(-0.5,0.5))
wilcox.test(deltaSNP[SNP.GR$ratio<0.5],deltaSNP[SNP.GR$ratio>0.5])




