library(bit)
library(AnnotationDbi)
library(org.Hs.eg.db)
###### Survival package
library(lattice)

######## plot KM curves
library(reshape2)
library(data.table)
library("zoo")
library("survminer")
library("survival")

#### ROC package
library(gplots)
library(ROCR)
library(pROC)
library(glmnet)
  

setwd("D:\\Academy\\data 2\\survival")

#### TCGA Targeted Therapy
# for the first time:
# Data_all = read.csv("Data/GBMLGG.csv")
# Gene_GBM=read.csv("Data/TCGA_GBM_GeneExpression.csv")
# cc = c(colnames(Gene_GBM))
# indices <- c(1)
# result <- cc[-indices]
# print(result)
# df = Data_all[, !(colnames(Data_all) %in% result)]
# write.csv(df,"Data/TCGA_LGG_GeneExpression.csv")
# rm(df)
# rm(Data_all)

Gene_LGG = read.csv("Data/TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,2]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-c(1,2)]

Survival_LGG=read.csv("Data/TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,93)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time','Target_Therapy','_primary_disease')

#####
Gene_GBM=read.csv("Data/TCGA_GBM_GeneExpression.csv")
rownames(Gene_GBM)=Gene_GBM[,1]
colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)
Gene_GBM=Gene_GBM[,-1]

Survival_GBM=read.csv("Data/TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,109)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','Target_Therapy')

Gene_TCGA=cbind(Gene_GBM[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),],Gene_LGG[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),])
colnames(Gene_TCGA)=gsub(".", "-", colnames(Gene_TCGA), fixed = TRUE)
Survival_TCGA=rbind(Survival_GBM, Survival_LGG)

Survival_TCGA=Survival_TCGA[Survival_TCGA[,4]=='YES',]
dim(Survival_TCGA)


#Preprocessing ranking result 
temp = read.table("Data//Quantifing each node_just_in HGG.txt", header = TRUE)
temp = read.table("Data//test.txt", header = TRUE)
temp = read.table("Data//test1.txt", header = TRUE)

gene_name =c()
for (i in 1:dim(temp)[1]) {
  gene_name = c(gene_name, unlist(strsplit(temp$name[i],".",fixed=TRUE))[1])
}
temp$name = gene_name
DNG=read.table("Data//Quantifing each node_just_in HGG.txt")
DNG=read.table("Data//test.txt")
DNG=read.table("Data//test1.txt")
# Differential network genes
Score=as.matrix(DNG[,c(1,6)])  # Importance Score
# Gene=rownames(Score)
Gene= Score[2:102,1]
Score= as.integer(Score[2:102,2])
names(Score)=c(Gene)

Gene1 = c()
for (i in 1:length(Gene)) {
  Gene1 = c(Gene1, unlist(strsplit(Gene[i],".", fixed = TRUE))[1])
}
names(Score)=c(Gene1)

#####################
#####################
### 100 computations
#####################
#####################

Time=100
sample_ratio=(4:8)*0.1

AUC_Training_wLASSO=matrix(0,length(sample_ratio),Time)
AUC_Test_wLASSO=matrix(0,length(sample_ratio),Time)

sr_ind=0  # index of sample ratio

for (p in sample_ratio){
  sr_ind=sr_ind+1
  
  for (i in 1:Time)
  {
    #########  Training set - TCGA dataset 1
    ind <- sample(2, nrow(Survival_TCGA), replace = TRUE, prob = c(p, 1-p))
    
    Survival_Gene=Gene_TCGA[c(names(Score)),which(ind==1)]
    
    colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
    
    
    Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
    Survival=as.matrix(Survival)
    Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]
    
    status=as.numeric(Survival[,2])
    time=as.numeric(Survival[,3])
    Genes=rownames(na.omit(Survival_Gene))
    #Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
    Gene_marker = Survival_Gene[Genes,]
    
    time=time[!is.na(status)]
    Gene_marker=Gene_marker[,!is.na(status)]
    status=status[!is.na(status)]
    
    status[status==1&time>365*3]=0  # 3-year survival status
    
    # S=as.numeric(Score)
    # Gene_marker_w=Gene_marker#*(S-min(S))/max(S)  # Importance score-weighted expression
    
    
    ##Dynamic Network-based LASSO
    S=as.numeric(Score)
    Gene_marker_w=Gene_marker  # Importance score-weighted expression
    cvfit_wLASSO <- cv.glmnet(t(Gene_marker_w),as.double(status),family="binomial",alpha=1,y.factor=S/max(S))
    #cvfit_wLASSO <- cv.glmnet(t(Gene_marker_w),as.double(status),family="binomial",alpha=1)
    lsso_coef = cvfit_wLASSO$glmnet.fit$beta[,cvfit_wLASSO$glmnet.fit$lambd ==cvfit_wLASSO$lambda.min]
    names(which(lsso_coef != 0))
    
    
    Prediction = predict(cvfit_wLASSO,newx = t(Gene_marker_w),s="lambda.min",type="response")
    #pred <- prediction(Prediction,status)
    # perf_Training_wLASSO[i] <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
    roc0=roc(status,as.numeric(Prediction))
    AUC_Training_wLASSO[sr_ind,i]=as.numeric(auc(status,as.numeric(Prediction)))
    
    
    
    ######### Validation set - TCGA dataset 2
    Survival_Gene=Gene_TCGA[names(Score),]
    Survival_Gene=Survival_Gene[,which(ind==2)]
    
    colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
    
    
    Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
    Survival=as.matrix(Survival)
    Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]
    
    status=as.numeric(Survival[,2])
    time=as.numeric(Survival[,3])
    Genes=rownames(na.omit(Survival_Gene))
    Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
    
    time=time[!is.na(status)]
    Gene_marker=Gene_marker[,!is.na(status)]
    status=status[!is.na(status)]
    
    status[status==1&time>365*3]=0  # 3-year survival status
    
    
    ### Validation for Dynamic Network-based LASSO
    S=as.numeric(Score)
    Gene_marker_w=Gene_marker#*(S-min(S))/max(S)  # Importance score-weighted expression
    
    Prediction = predict(cvfit_wLASSO,newx = t(Gene_marker_w),s="lambda.min",type="response")
    #pred <- prediction(Prediction,status)
    # perf_Test_wLASSO[i] <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
    roc1=roc(status,as.numeric(Prediction))
    AUC_Test_wLASSO[sr_ind,i]=as.numeric(auc(status,as.numeric(Prediction)))
    
  }
}


# AUC_Training_wLASSO

AUC_Test_wLASSO
# AUC_Test_wLASSO = AUC_Test_wLASSO[,1:82]
# AUC_Training_wLASSO = AUC_Training_wLASSO[,1:82]
mean(AUC_Training_wLASSO)
# 0.9709169
mean(AUC_Test_wLASSO)
# 0.8852524
#### Barplots of AUCs_Test dataset of Dynamic Network-based LASSO
dev.new()
boxplot(t(AUC_Training_wLASSO),ylim=c(0,1), col="#185998")
# mean training = 91%

dev.new()
boxplot(t(AUC_Test_wLASSO),ylim=c(0,1),col="#126a71")

# Finished successfully

############### Differential expression profiles of 7 genes in the Normal and Tumor tissues

Gene_LGG=read.csv("Data/TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,2]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-c(1,2)]
# 
Survival_LGG=read.csv("Data/TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,87)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time','Target_Therapy')

#####
Gene_GBM=read.csv("Data/TCGA_GBM_GeneExpression.csv")
rownames(Gene_GBM)=Gene_GBM[,1]
colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)
Gene_GBM=Gene_GBM[,-1]

Survival_GBM=read.csv("Data/TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,106)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','Sample Type')


Gene_TCGA=Gene_GBM
colnames(Gene_TCGA)=gsub(".", "-", colnames(Gene_TCGA), fixed = TRUE)
Survival_TCGA=Survival_GBM[intersect(colnames(Gene_GBM),rownames(Survival_GBM)),]

Gene_Normal=Gene_TCGA[,Survival_TCGA[,4]=='Solid Tissue Normal']
Gene_Tumor=Gene_TCGA[,Survival_TCGA[,4]!='Solid Tissue Normal']
dim(Gene_Normal)
dim(Gene_Tumor)

Gene_TCGA=cbind(Gene_LGG,Gene_GBM)
Gene_Tumor=Gene_TCGA[,setdiff(colnames(Gene_TCGA),colnames(Gene_Normal))]
#rownames(Gene_Tumor)=Gene_Tumor[,1]
MarkerGenes=c("MYC"	,"EGFR", "ETS1","NOTCH3",	"YBX1","LIFR","ETS2", "GNAI2")
#Gene_Tumor = Gene_Tumor[,-1]
#MarkerGenes=c("KIF2C", "CCNA2", "NDC80", "KIF11", "KIF23", "ANLN", "CENPM")
# MarkerGenes=c("MYC", "ETS1", "FOS")
MarkerGenes= c("MYC", "EGR1", "ETS1", "NOTCH3", "EGFR", "CREB1", "CAV1", "YBX1", "PDGFRA", "FOS"
)


# ID=ensemble2symbol[intersect(rownames(ensemble2symbol),GeneName),] 

Marker_Normal=matrix(as.numeric(as.matrix(Gene_Normal[MarkerGenes,])),dim(Gene_Normal[MarkerGenes,]))
Marker_Tumor= matrix(as.numeric(as.matrix(Gene_Tumor[MarkerGenes,])),dim(Gene_Tumor[MarkerGenes,]))

pvalue=matrix(0,1,length(MarkerGenes))
for (i in 1:length(MarkerGenes))
{
  Wtest=wilcox.test(Marker_Normal[i,],Marker_Tumor[i,])
  pvalue[i]=Wtest$p.value
}

Data=as.data.frame(cbind(Marker_Normal,Marker_Tumor))
Data=t(na.omit(Data))

type=as.factor(c(rep('Normal',5),rep('Tumor',697)))
x=cbind(Data,(type))

library(vioplot)

D=x[x[,5]==1,1:4]
E=x[x[,5]==2,1:4]
dev.new()

vioplot(D[,1],D[,2],D[,3],D[,4], names=paste(as.character(MarkerGenes), 'Normal',sep='_'),col="DarkKhaki", border="gray51",las=2, lty=1, lwd=1, rectCol="gray", 
        colMed="white", pchMed=19)


vioplot(E[,1],E[,2],E[,3],E[,4], add=T, names=paste(as.character(MarkerGenes), 'Tumor',sep='_'), las=2,col="DeepSkyBlue3", border="gray51", lty=1, lwd=1, rectCol="gray", 
        colMed="white", pchMed=19)



##### New by Fereshteh:
# 
# MarkerGenes = c("PDGFA",   "PDGFRA" ,   "CREB1"  ,  
#                 "ETS1", "LAPTM4B" , "NOLC1", "PLAT" ,
#                 "PPM1D", "SLC19A1", "SOX9","FN1", 
#                 "CSNK2A1", "EGR1", "ERG", "RHOC", "TFAP2A", "TMEM158", "TNC")
# 
MarkerGenes = c("CCND2", "DST", "MYC", "CASP9", "EGR1",
                "FOS", "IGF1R",
                "EGFR", "CAV1", "GNAI2", "YBX1", "ETV4", "RUNX1", "LIFR",
                "ETS2", "ITGA2","PDGFRA", "RBL1")
MarkerGenes = c("CDK6","HDAC2","HNRNPA1","TOP1MT")
# Pathway Examination:
MarkerGenes=c("PDGFA"	,"PDGFRA", "CREB1","EGR1",	"PLAT","EGR1")
MarkerGenes=c("UBA52"	,"EGFR", "FOS",	"MYC","MYC")


Marker_LGG=matrix(as.numeric(as.matrix(Gene_LGG[MarkerGenes,])),dim(Gene_LGG[MarkerGenes,]))
Marker_HGG= matrix(as.numeric(as.matrix(Gene_GBM[MarkerGenes,])),dim(Gene_GBM[MarkerGenes,]))

pvalue=matrix(0,1,length(MarkerGenes))
for (i in 1:length(MarkerGenes)){Wtest=wilcox.test(Marker_LGG[i,],Marker_HGG[i,])
  pvalue[i]=Wtest$p.value}
p.adjust(pvalue,method = "bonferroni")

Data=as.data.frame(cbind(Marker_LGG,Marker_HGG))
Data=t(na.omit(Data))

type=as.factor(c(rep('LGG',530),rep('HGG',172)))
x=cbind(Data,(type))

library(vioplot)

D=x[x[,11]==1,1:10]
E=x[x[,11]==2,1:10]
dev.new()

vioplot(D[,1],D[,2],D[,3],D[,4],D[,5],D[,6],D[,7],D[,8],D[,9],D[,10], names=paste(as.character(MarkerGenes), 'LGG',sep='_'),col="DarkKhaki", border="gray51", lty=1, lwd=1, rectCol="gray", colMed="white", pchMed=19)


vioplot(E[,1],E[,2],E[,3],E[,4],E[,5],E[,6],E[,7],E[,8],E[,9],E[,10], names=paste(as.character(MarkerGenes), 'HGG',sep='_'),col="DeepSkyBlue3", border="gray51", lty=1, lwd=1, rectCol="gray",
        colMed="white", pchMed=19)



##################################################################################################################
##################################################################################################################

#################  K-M analysis

##################  ensemble ID  to gene symbol   #################################
ensembl2gene <- toTable(org.Hs.egENSEMBL2EG)
gene2symbol <- toTable(org.Hs.egSYMBOL)
ensemble2symbol <- merge(ensembl2gene, gene2symbol, by = 'gene_id')[2:3]
ensemble2symbol=as.matrix(ensemble2symbol)
rownames(ensemble2symbol)=ensemble2symbol[,2]

###TCGA data
Gene_LGG=read.csv("Data/TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,2]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-c(1,2)]

Survival_LGG=read.csv("Data/TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time')

#####
Gene_GBM=read.csv("Data/TCGA_GBM_GeneExpression.csv")
rownames(Gene_GBM)=Gene_GBM[,1]
colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)
Gene_GBM=Gene_GBM[,-1]

Survival_GBM=read.csv("Data/TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time')

Gene_TCGA=cbind(Gene_GBM[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),],Gene_LGG[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),])
colnames(Gene_TCGA)=gsub(".", "-", colnames(Gene_TCGA), fixed = TRUE)
Survival_TCGA=rbind(Survival_GBM,Survival_LGG)


######################## Select genes for signature   #####################################

#MNB=c("KIF2C", "CCNA2", "NDC80", "KIF11", "KIF23", "ANLN", "CENPM")
# MNB=c("MYC", "CAV1",	"EGFR",	"EGR1",	"ETS1",	"FOS",	"NOTCH3",	"TRAF2",	"YBX1",	"PDGFRA"	,"CREB1",	"ETS2"	,"RUNX1",	"MC1R",	"ITGA2")

#MNB = c("MYC"	,"EGFR", "ETS1","NOTCH3",	"YBX1","LIFR","ETS2", "GNAI2")
MNB = c("MYC", "EGR1", "ETS1", "NOTCH3", "EGFR", "CREB1", "CAV1", "YBX1", "PDGFRA", "FOS")
MNB = MarkerGenes

MNB=c("AREG"	,"EGFR", "FOS",	"MYC","MYC")

MNB=c("PDGFA"	,"PDGFRA", "CREB1","EGR1","PLAT")

Survival_Gene=Gene_TCGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

### Training dataset
status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))



time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

### COX Model
surv = Surv(as.double(time),as.double(status))
fit <- coxph(surv ~ t(Gene_marker), method="efron")
A = coef(fit)
A

# 
# surv=Surv(as.double(time),as.double(status))
# 
# cvfit<-cv.glmnet(t(Gene_marker),surv,family="cox",alpha=1)
# plot(cvfit)
# A=coef(cvfit, s = "lambda.min")
# A=coef(cvfit, s="lambda.1se")
# which(A!=0)
# sum(A!=0)
# MNB[which(A!=0)]
# 
# fit <- coxph(Surv(time, status) ~ t(Gene_marker[which(A!=0),]), method="breslow")

########  ROC for Training set
predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]
status0=as.matrix(na.omit(status))
pred <- prediction(predicted0,status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="DarkCyan") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc0=roc((status0),predicted0,smooth=F,plot=TRUE)
roc0
AUC=auc(roc0,predicted0,partial.auc=FALSE)
AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc0$sensitivities,roc0$specificities)))
## optimal cut-off point 
sort(predicted0,F)[opt]

##########  K-M survival curves for MNB in TCGA datasets
groups=matrix(0,1,length(status))
groups[predicted<=sort(predicted,F)[opt]]=1
groups[predicted>sort(predicted,F)[opt]]=2
# groups[predicted<=median(predicted)]=1
# groups[predicted>median(predicted)]=2
groups=t(groups)
groups=as.numeric(groups)


fit<- survfit(Surv(time/365, status) ~ groups, data = as.data.frame(Survival_Gene))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))




# Add risk table
# and change risk table y text colors by strata
ggsurvplot(fit, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, risk.table.y.text.col = TRUE)



