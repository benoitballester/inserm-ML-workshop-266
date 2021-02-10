

# conda activate /home/mourad/.local/share/r-miniconda/envs/r-reticulate
# install_tensorflow(method = 'conda', envname = '/home/raphael/.local/share/r-miniconda/envs/r-reticulate')


#### LOAD LIBRARIES

library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)
library(tensorflow)
library(keras)
library(reticulate)
library(data.table)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


#### SETUP PROJECT FOLDER

#setwd("/media/mourad/diskSave/MCF_Toulouse/enseignement/atelier_INSERM/")
setwd("/media/raphael/SSD2/atelier_INSERM")


#### LOAD FUNCTIONS

source("scriptR/functions.R")

use_condaenv('r-reticulate')
tf$constant("Hello Tensorflow")


#### CREATE FOLDERS


dir.create("results/SNP")


#### SOME PARAMETERS

peakSize=201
kpeaks=4000
expe="DNase"
DNAletters=c("A","C","G","T")

SNPanalysis="eSNP" # "eSNP" or "gwasSNP"


#### INPUT FILE NAMES

fileBedPos=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_pos.bed")
fileBedNeg=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_neg.bed")
fileFastaPos=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_pos.fa")
fileFastaNeg=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_neg.fa")




#### DEEP LEARNING PREDICTIONS OF SNP EFFECTS

# Load CNN1 model
file_model=paste0("results/model/CNN1_model_",expe,".hdf5")
model=load_model_hdf5(file_model)


# Load expression SNP (eSNP) data from GTEx project (hg19)
# https://gtexportal.org/home/datasets
if(SNPanalysis=="eSNP"){
 genomeVersion=BSgenome.Hsapiens.UCSC.hg19
 
 # Load and process breast mammary tissue eSNPs
 # variant_id: variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b37
 # slope of the linear regression: is computed as the effect of the alternative allele (ALT) relative to the reference allele (REF)
 dataSNP=as.data.frame(fread("data/SNP/Breast_Mammary_Tissue.v7.signif_variant_gene_pairs.txt.gz",header=T))
 #dataSNP=as.data.frame(fread("data/SNP/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz",header=T))
 dataSNP=dataSNP[abs(dataSNP$slope)>1,] # Select SNPs with high absolute effects
  
 # Extract chromosome, position and alleles
 variantInfo=t(sapply(as.character(dataSNP[,1]),function(x){as.character(strsplit(x,'_')[[1]][1:4])}))
 rownames(variantInfo)=rep("",nrow(variantInfo))
 
 # Create a GRanges object
 SNP.GR=GRanges(paste0("chr",variantInfo[,1]),IRanges(as.numeric(variantInfo[,2]),as.numeric(variantInfo[,2])),ref=variantInfo[,3],alt=variantInfo[,4],
	slope=dataSNP[,8],tss_distance=dataSNP$tss_distance, gene=dataSNP$gene_id) 

 # Only pick SNPs (remove small indels)
 SNP.GR=SNP.GR[nchar(SNP.GR$ref)==1 & nchar(SNP.GR$alt)==1]
 
 # Select SNPs in promoters
 #txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
 #prom.GR=promoters(genes(txdb), upstream = 1500, downstream = 500)
 #olprom=findOverlaps(SNP.GR,prom.GR)
 #SNP.GR=SNP.GR[queryHits(olprom)]
 
}


# Load GWAS SNP (gwasSNP) data from GWAS catalog (hg38)
# https://www.ebi.ac.uk/gwas/docs/file-downloads
if(SNPanalysis=="gwasSNP"){
 genomeVersion=BSgenome.Hsapiens.UCSC.hg38
 
 # Load GWAS catalog SNPs
 dataGWAS=as.data.frame(fread("data/SNP/gwas_catalog_v1.0-associations_e100_r2021-01-29.tsv.gz",header=T))
 
 # Remove SNP with OR=Inf or OR=NA
 dataGWAS=dataGWAS[dataGWAS[,11]!=Inf & !is.na(dataGWAS[,11]),]
 
 # Convert odds ratios to betas (slopes)
 dataGWAS$slope=log(dataGWAS[,11])
 
 # Remove SNPs without known positions
 dataGWAS=dataGWAS[dataGWAS$CHR_POS!="",]
 
 # Remove SNPs without an rs number
 dataGWAS=dataGWAS[substring(dataGWAS$SNPS,1,2)=="rs",]
 
 # Extract risk allele and only keep A, T, C or C alleles
 riskAllele=sapply(dataGWAS[,7],function(x){strsplit(x,"-")[[1]][2]})
 dataGWAS=dataGWAS[riskAllele%in%DNAletters,]
 riskAllele=riskAllele[riskAllele%in%DNAletters]

 # Create a GRanges object
 SNP.GR=GRanges(paste0("chr",dataGWAS[,4]),IRanges(as.numeric(dataGWAS$CHR_POS),as.numeric(dataGWAS$CHR_POS)),alt=riskAllele,slope=dataGWAS$slope) 
 SNP.GR=SNP.GR[!duplicated(SNP.GR)]
 refAllele=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, SNP.GR))
 SNP.GR$ref=refAllele
 SNP.GR=SNP.GR[SNP.GR$ref!=SNP.GR$alt]
 
}

SNP.GR=GRanges(c("chr1"),IRanges(109821511,109821511),ref="G",alt="T",slope="higher")

# rs9303277 is found from GWAS in childhood acute lymphoblastic leukemia
# Ref: https://www.nature.com/articles/s41467-017-02596-9#MOESM1
# Also an eQTL in blood (https://www.gtexportal.org/home/snp/rs9303277)
SNP.GR1=GRanges("chr17",IRanges(37976469),ref="C",alt="T",rs="rs9303277",OR_GWAS=0.86, Effet_eQTL=-0.33)
SNP.GR=SNP.GR1

# Extract sequences surrounding the SNPs
region.GRr=resize(SNP.GR,peakSize,fix="center")
region.seq=as.character(getSeq(genomeVersion, names=seqnames(region.GRr),
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
#writeXStringSet(DNAStrinregionMaxProbgSet(region.seqAlt),filepath=paste0(pathfasta,"/sequences_SNPalt_hg19.fa"))

# Predict SNP effect on the predicted probability using the deep learning model
oneHotRef=convertOneHot(region.seqRef)
oneHotAlt=convertOneHot(region.seqAlt)
predProbRef=predict(model,oneHotRef)
predProbAlt=predict(model,oneHotAlt)
deltaProbSNP=predProbAlt-predProbRef

# Predict SNP effect on the second to last layer (just before logistic transformation) using the deep learning model
secondlast_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, '2ndlastlayer_cnn1')$output)
secondlast_outputRef <- predict(secondlast_layer_model, oneHotRef)
secondlast_outputAlt <- predict(secondlast_layer_model, oneHotAlt)
delta2ndLastSNP=secondlast_outputAlt-secondlast_outputRef

# Sequences with high proba
regionMaxProb=sapply(1:length(predProbRef),function(x){max(predProbRef[x],predProbAlt[x])})

# Compare model predicted SNP effects between SNPs with slope < 1 and SNPs woth slope > 1
deltaSNPnegSlope=delta2ndLastSNP[SNP.GR$slope<1]# & regionMaxProb>0.2]
deltaSNPposSlope=delta2ndLastSNP[SNP.GR$slope>1]# & regionMaxProb>0.2]
if(SNPanalysis=="eSNP"){
 dataDeltaSNP=data.frame(slope=c(rep("Expression decreased",length(deltaSNPnegSlope)),rep("Expression increased",length(deltaSNPposSlope))),
			delta=c(deltaSNPnegSlope,deltaSNPposSlope))
}else if(SNPanalysis=="gwasSNP"){
 dataDeltaSNP=data.frame(slope=c(rep("Lower disease risk",length(deltaSNPnegSlope)),rep("Higher disease risk",length(deltaSNPposSlope))),
			delta=c(deltaSNPnegSlope,deltaSNPposSlope))
}

wt=wilcox.test(deltaSNPnegSlope,deltaSNPposSlope)
wt$p.value=format.pval(wt$p.value)
wt$p.value

file_violinplot=paste0("results/SNP/violinplot_",SNPanalysis,"_",expe,".pdf")
pdf(file_violinplot)
p <- ggplot(dataDeltaSNP, aes(x=slope, y=delta)) + 
  geom_violin() + geom_boxplot(width=0.1) + coord_cartesian(ylim=c(-0.05,0.05)) + 
  ggtitle(paste0("Pvalue=",wt$p.value)) + geom_hline(yintercept=0, linetype="dashed", color="red") 
print(p)
dev.off()

median(deltaSNPnegSlope)
median(deltaSNPposSlope)







