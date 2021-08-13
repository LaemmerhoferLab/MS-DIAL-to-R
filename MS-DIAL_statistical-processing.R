#Description
#open new project in folder that contains only alignment output as .txt file and no other .txt files (in case of .csv files or other formats the "sep" parameter has to be changed)
#make sure you have only 3 groups (patients, controls, QCs); Blank samples will be automatically removed


#set parameters (e.g. TRUE/FALSE)

#which statistical test would you like to do? U.Test (wilcox.test) or T.test? Please select only one as TRUE.
#Other tests like Welch test have to be set by different parameters in the function
#fold change of U.test will be calculated based on median, fold change of t.test will be based on mean
U.test<-FALSE
T.test<-TRUE

#which groups are you interested in? Enter which groups must be removed to remain 2 experimental groups
remove.groups<-"3|4|5|6"

#impute missing values via random Forest method
impute.miss=TRUE

#is data already log-transformed? (affects calculation of fold change)
log.transformed.data=FALSE
#if log.transformed.data=TRUE, to which base was it transformed (e.g. exp(1),10,2,etc...)
base=exp(1)

#log-transform data (default: transformation to base=exp(1) (natural logarithm), if transformation to other base or other scaling methods are anticipated please change this parameter)
log.transform=TRUE

#which p-value adjustment method would you like to do? Choose from "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", "SGOF"
p.adjust.method="bonferroni"

#draw boxplots/ROC curve if p-value<0.05 (FALSE) or if adjusted p-value<0.05 (TRUE)?
only.adjusted.p.values=TRUE


#load packages
suppressMessages(library(missForest))
suppressMessages(library(NormalizeMets))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(ggfortify))
suppressMessages(library(pwr))
suppressMessages(library(sgof))
suppressMessages(library(pROC))
suppressMessages(library(ggsci))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(beanplot))



#load and arrange data
data.MSDIAL <- read.delim(list.files(pattern="\\.txt")[1], sep="\t", header=FALSE, dec=".",na.strings=c("NA", ""),stringsAsFactors=FALSE)
data.feature.info<-data.MSDIAL[,!grepl("Blank",data.MSDIAL[2,])]

#remove groups (2 experimental groups + QCs allowed)
#QC.group<-data.MSDIAL[1,which(data.MSDIAL[2,]=="QC")[1]]
#my.groups<-c(QC.group,my.groups)
data.feature.info<-data.feature.info[,!grepl(remove.groups,data.feature.info[1,])]
#data.feature.info<-data.feature.info[,data.feature.info[1,] %in% my.groups]

if (length(which(data.feature.info[4,]=="Average"))>0)
{
  data.feature.info<-data.feature.info[,(which(data.feature.info[4,]=="Average")[1]*-1):-ncol(data.feature.info)]
}
data.feature.info<-data.feature.info[-1:-(which(data.feature.info[,1]!="")[1]-1),]

rownames(data.feature.info)<-data.feature.info[,1]
data.feature.info<-data.feature.info[,-1]
colnames(data.feature.info)<-data.feature.info[1,]
data.feature.info<-as.data.frame(t(data.feature.info[-1,]))
data.feature.info.trans<-as.data.frame(t(data.feature.info))

data.raw.original<-data.feature.info[-1:-(which(data.MSDIAL[1,]!="")[1]-1),]
data.raw.original[]<-as.data.frame(apply(data.raw.original,2,as.numeric))

data.sample.info<-as.data.frame(t(data.MSDIAL[1:(which(data.MSDIAL[,1]!="")[1]),-1:-(which(data.MSDIAL[1,]!="")[1]-1),]))
data.sample.info<-data.sample.info[grepl("",data.sample.info[,1]),]
data.sample.info<-data.sample.info[!grepl("Blank",data.sample.info[,2]),]
data.sample.info<-data.sample.info[!grepl(remove.groups,data.sample.info[,1]),]
rownames(data.sample.info)<-c()
colnames(data.sample.info)<-c(data.MSDIAL[1:(which(data.MSDIAL[,1]!="")[1]-1),which(data.MSDIAL[1,]!="")[1]],"Sample")
data.sample.info<-data.sample.info[-1,]
data.sample.info[]<-as.data.frame(apply(data.sample.info,2,as.factor))

Class<-data.sample.info$Class
Class.samples<-levels(factor(data.sample.info$Class[which(data.sample.info$`File type`=="Sample")]))
Filetype<-data.sample.info$`File type`
InjectionID<-as.numeric(data.sample.info$`Injection order`)


#impute missing values and replace 0 by 1e-12, think about doing missing value imputation before or after removing groups
data.raw<-data.raw.original
data.raw[data.raw.original==0]<-1e-12

if (impute.miss==TRUE) 
{data.raw<-missForest(data.raw, verbose=TRUE)$ximp}


#log-transformation
if (log.transform==TRUE) {data.raw<-log(data.raw,base=exp(1))}


#sort dataset according to class and injectionorder
data.raw.IDsorted<-data.raw[order(InjectionID),]
Class.IDsorted<-data.sample.info[order(as.numeric(data.sample.info$`Injection order`)),]$Class
data.raw.Classsorted<-data.raw[order(Class),]
Class.Classsorted<-data.sample.info[order(as.numeric(data.sample.info$Class)),]$Class


#create RLA-plots (within-group, across-group, considering injection order)
dir.create("RLA plots")
setwd("RLA plots")

tiff("RLAPlot_ag_raw_IDsorted.tiff", width=174, height=100, units="mm", res=600, compression = "lzw")
par(oma = c(1,0,0,0), mar = c(1.5,1.5,0.1,0.1),mgp=c(0,0.7,0),cex=0.5)
RlaPlots(data.raw.IDsorted, Class.IDsorted, "ag", outline=FALSE,ylim=NULL, interactiveplot = FALSE,keeporder=TRUE)
dev.off()

tiff("RLAPlot_wg_raw_IDsorted.tiff", width=174, height=100, units="mm", res=600, compression = "lzw")
par(oma = c(1,0,0,0), mar = c(1.5,1.5,0.1,0.1),mgp=c(0,0.7,0),cex=0.5)
RlaPlots(data.raw.IDsorted, Class.IDsorted, "wg", outline=FALSE,ylim=NULL, interactiveplot = FALSE,keeporder=TRUE)
dev.off()

tiff("RLAPlot_ag_raw.tiff", width=174, height=100, units="mm", res=600, compression = "lzw")
par(oma = c(1,0,0,0), mar = c(1.5,1.5,0.1,0.1),mgp=c(0,0.7,0),cex=0.5)
RlaPlots(data.raw.Classsorted, Class.Classsorted, "ag", outline=FALSE,ylim=NULL, interactiveplot = FALSE,keeporder=TRUE)
dev.off()

tiff("RLAPlot_wg_raw.tiff", width=174, height=100, units="mm", res=600, compression = "lzw")
par(oma = c(1,0,0,0), mar = c(1.5,1.5,0.1,0.1),mgp=c(0,0.7,0),cex=0.5)
RlaPlots(data.raw.Classsorted, Class.Classsorted, "wg", outline=FALSE,ylim=NULL, interactiveplot = FALSE,keeporder=TRUE)
dev.off()
setwd("..")


#precision of features in QCs
precQC<-apply(subset(data.raw,Filetype=="QC"),2,function(x) sd(x)/mean(x)*100)
quantile.precQC<-quantile(precQC, probs=c(.25,.75))
iqr.precQC<-IQR(precQC)
precQC.wo.outliers<-subset(precQC, precQC>(quantile.precQC[1]-1.5*iqr.precQC)&precQC>(quantile.precQC[2]+1.5*iqr.precQC))

dir.create("QC precision")
setwd("QC precision")

tiff(paste("QC_precision_beanplot.tiff",sep=""), width=84, height=84, units="mm", res=600, compression="lzw")
par(oma = c(0,0,0,0), mar = c(3.5,4.5,4,0.5),mgp=c(3,0.7,0),cex=0.5)
beanplot(precQC, ylab = "precision [%]", xlab="QC")
mtext("median precision [%]:", side=3, adj=0,cex=0.5)
mtext(paste(round(median(precQC),1)," ? ",round(sd(precQC),1)), side=3, adj=0.38,cex=0.5)
dev.off()

tiff(paste("QC_precision_wo-outliers_beanplot.tiff",sep=""), width=84, height=84, units="mm", res=600, compression="lzw")
par(oma = c(0,0,0,0), mar = c(3.5,4.5,4,0.5),mgp=c(3,0.7,0),cex=0.5)
beanplot(precQC.wo.outliers, ylab = "precision [%]", xlab="QC")
mtext("median precision [%]:", side=3, adj=0,cex=0.5)
mtext(paste(round(median(precQC),1)," ? ",round(sd(precQC),1)), side=3, adj=0.38,cex=0.5)
dev.off()

tiff(paste("QC_precision_boxplot.tiff",sep=""), width=84, height=84, units="mm", res=600, compression="lzw")
par(oma = c(0,0,0,0), mar = c(3.5,4.5,4,0.5),mgp=c(3,0.7,0),cex=0.5)
boxplot(precQC,notch=FALSE,boxwex=0.6, ylab = "precision [%]", xlab="QC", las = 1, family = "Arial",outline=FALSE)
stripchart(precQC, vertical = TRUE, method = "jitter", add = TRUE, pch = 1, col = 'brown',cex=0.4)
dev.off()

tiff(paste("QC_precision_wo-outliers_boxplot.tiff",sep=""), width=84, height=84, units="mm", res=600, compression="lzw")
par(oma = c(0,0,0,0), mar = c(3.5,4.5,4,0.5),mgp=c(3,0.7,0),cex=0.5)
boxplot(precQC,notch=FALSE,boxwex=0.6, ylab = "precision [%]", xlab="QC", las = 1, family = "Arial",outline=TRUE)
stripchart(precQC, vertical = TRUE, method = "jitter", add = TRUE, pch = 1, col = 'brown',cex=0.2)
dev.off()

setwd("..")


#create PCA plot
PCA_all<-autoplot(prcomp(data.raw),data=as.data.frame(Class),colour="Class",frame=TRUE, frame.type="t")
tiff("PCA_all.tiff", width=174, height=174, units="mm", res=600, compression = "lzw")
print(PCA_all)
dev.off()

Classes<-subset(Class,Filetype=="Sample")
PCA_samples<-autoplot(prcomp(subset(data.raw,Filetype=="Sample")),data=as.data.frame(Classes),colour="Classes",frame=TRUE, frame.type="t")
tiff("PCA_samples.tiff", width=174, height=174, units="mm", res=600, compression = "lzw")
print(PCA_samples)
dev.off()


#compute p-values (unpaired Wilcox-Test, change parameters if other tests are needed. exact=FALSE in case of tie values in data)
if (U.test==TRUE) {stat.test=wilcox.test}
if (T.test==TRUE) {stat.test=t.test}

stat.results<-t(apply(subset(data.raw,Filetype=="Sample"),2,function(x) 
  stat.test(x~subset(Class,Filetype=="Sample"),var.equal=FALSE,paired=FALSE,exact=TRUE, conf.int=TRUE)))
p.values<-sapply(stat.results, function(x){as.numeric(x[3])})
stat.value<-sapply(stat.results, function(x){as.numeric(x[1])})


#p-value histogram
tiff(paste("p-value histogram ",Class.samples[1], " vs ",Class.samples[2],".tiff",sep=""), width=84, height=84, units="mm", res=600, compression="lzw")
par(oma = c(0,0,0,0), mar = c(3,4,1,0.5),mgp=c(3,0.7,0),cex=0.5)
hist(p.values, main=NULL ,breaks = 20, xlab=NULL, las = 1, col="gray92")
title(xlab =expression(italic(p)~value), line=2)
dev.off()



#p-value adjustment
if (p.adjust.method=="SGOF")
{p.values.adj<-as.data.frame(cbind(p.values,seq.int(nrow(as.data.frame(p.values)))))
p.values.adj<-p.values.adj[order(p.values),]
p.values.adj[,3]<-SGoF(p.values.adj$p.values)$Adjusted.pvalues
p.values.adj<-p.values.adj[order(p.values.adj[,2]),]
p.values.adj<-as.vector(p.values.adj[,3])
} else
{
  p.values.adj<-p.adjust(p.values, method=p.adjust.method)
}


#calculate fold change
if (log.transformed.data==TRUE | log.transform==TRUE)
{
  if (U.test==TRUE)
  {foldchange<-base^(apply(data.raw[grepl(Class.samples[1],data.sample.info$Class),],2, FUN=median)-apply(data.raw[grepl(Class.samples[2],data.sample.info$Class),],2, FUN=median))
  }
  if (T.test==TRUE)
  {foldchange<-base^(apply(data.raw[grepl(Class.samples[1],data.sample.info$Class),],2, FUN=mean)-apply(data.raw[grepl(Class.samples[2],data.sample.info$Class),],2, FUN=mean))
  }
}
if (log.transformed.data==FALSE & log.transform==FALSE)
{if (U.test==TRUE)
{foldchange<-apply(data.raw[grepl(Class.samples[1],data.sample.info$Class),],2, FUN=median)/apply(data.raw[grepl(Class.samples[2],data.sample.info$Class),],2, FUN=median)
}
  if (T.test==TRUE)
  {foldchange<-apply(data.raw[grepl(Class.samples[1],data.sample.info$Class),],2, FUN=mean)/apply(data.raw[grepl(Class.samples[2],data.sample.info$Class),],2, FUN=mean)
  }
}


#power calculation
if (T.test==TRUE)
{n1=length(which(data.sample.info$Class==Class.samples[1]))
n2=length(which(data.sample.info$Class==Class.samples[2]))
Cohens.d<-(apply(data.raw[grepl(Class.samples[1],data.sample.info$Class),],2, FUN=mean)-apply(data.raw[grepl(Class.samples[2],
                                                                                                             data.sample.info$Class),],2, FUN=mean))/sqrt(((apply(data.raw[grepl(Class.samples[1],data.sample.info$Class),],2, FUN=sd)^2+
                                                                                                                                                              apply(data.raw[grepl(Class.samples[2],data.sample.info$Class),],2, FUN=sd)^2)/2))
power.t<-pwr.t2n.test(n1,n2,d=Cohens.d,sig.level=0.05,power=NULL)$power
stat.power<-power.t} else
{
  Cohens.d<-999
  power.U<-999
  stat.power<-power.U
}
#power calc. for wilcox test?
#power.U<-
stat.power<-power.U
#temporary solution


#calculate AUC
suppressMessages(AUC<-sapply(apply(subset(data.raw,Filetype!="QC"),2,function(x) roc(factor(subset(Class,Filetype!="QC")),x)), function(x){as.numeric(x[9])}))     


#Summary of statistical results
stat.summary<-data.frame(rownames(data.feature.info.trans),data.feature.info.trans$`Average Rt(min)`,as.numeric(as.character(data.feature.info.trans$`Average Mz`)),as.character(data.feature.info.trans$`Metabolite name`),
                         as.character(data.feature.info.trans$`Adduct type`),round(stat.value,2),round(p.values,5),round(p.values.adj,5),round(foldchange,4),round(Cohens.d,4),round(stat.power,4),round(AUC,4),round(precQC,1),stringsAsFactors=FALSE)
colnames(stat.summary)<-c("ID","RT [min]","m/z","Identification","Adduct type","Statistical Value","p-value",paste("adjusted p-value (",p.adjust.method,")",sep=""),paste("fold change",Class.samples[1], " vs ",Class.samples[2]),"Cohens d-value", "Power","AUC","Precision in QC [%]")
rownames(stat.summary)<-c()
write.table(stat.summary, paste("statistical results ",Class.samples[1], " vs ",Class.samples[2],".csv",sep=""), sep = ",", quote = FALSE, append = FALSE, row.names = TRUE, col.names = NA)


#boxplots of significant features
dir.create("boxplots")
setwd("boxplots")

if (only.adjusted.p.values==FALSE) {draw.plots=p.values} else {draw.plots=p.values.adj}
sigrow<-which(draw.plots<0.05)
signames<-paste("ID: ",stat.summary$ID[sigrow],"   ", " RT: ", stat.summary$`RT [min]`[sigrow],"   ", " m/z: ",stat.summary$`m/z`[sigrow], "    ",stat.summary$Identification[sigrow],sep="")

if (length(sigrow)>0)
{
  for (z in 1:length(sigrow))
  {		    
    tiff(paste(stat.summary$ID[sigrow[z]],"_", stat.summary$`RT [min]`[sigrow[z]],"_", stat.summary$`m/z`[sigrow[z]],"Boxplot.tiff",sep=""), width=84, height=84, units="mm", res=600, compression="lzw")
    par(oma = c(0,0,0,0), mar = c(5,4.5,4,0.5),mgp=c(3,0.7,0),cex=0.5)
    boxplot(data.raw[,sigrow[z]] ~ Class,notch=FALSE,boxwex=0.6, ylab = "Intensity", main = signames[z], las = 1, family = "Arial",outline=FALSE)
    stripchart(data.raw[,sigrow[z]] ~ Class, vertical = TRUE, method = "jitter", add = TRUE, pch = 1, col = 'brown',cex=0.9)
    mtext("p-value:", side=3, adj=0,cex=0.5)
    mtext(round(p.values[sigrow[z]],4), side=3, adj=0.13,cex=0.5)
    mtext("adj. p-val:", side=3, adj=0.31,cex=0.5)
    mtext(round(p.values.adj[sigrow[z]],4), side=3, adj=0.44,cex=0.5)
    mtext("fold change:", side=3, adj=0.65,cex=0.5)
    mtext(round(foldchange[sigrow[z]],2), side=3, adj=0.76,cex=0.5)
    mtext("power:", side=3, adj=0.92,cex=0.5)
    mtext(round(stat.power[sigrow[z]],2), side=3, adj=0.97,cex=0.5)
    dev.off()
  }
}
setwd("..")


#ROC curves and AUCs of significant features
dir.create("ROC_AUC")
setwd("ROC_AUC")

if (length(sigrow)>0)
{
  for (z in 1:length(signames))
  {		    
    tiff(paste(stat.summary$ID[sigrow[z]],"_", stat.summary$`RT [min]`[sigrow[z]],"_", stat.summary$`m/z`[sigrow[z]],"ROC.tiff",sep=""), width=84, height=84, units="mm", res=600, compression="lzw")
    par(oma = c(0,0,0,0), mar = c(3.5,4.5,4,0.5),mgp=c(3,0.7,0),cex=0.5)
    suppressMessages(plot.roc(factor(subset(Class,Filetype!="QC")),subset(data.raw,Filetype!="QC")[,sigrow[z]],legacy.axes=TRUE,col="darkgreen",lwd=2,print.auc=TRUE,print.auc.x=0.25,print.auc.cex=1.5))
    dev.off()
  }
}
setwd("..")


#Volcano plot
tiff(paste("Volcano plot ",Class.samples[1], " vs ",Class.samples[2],".tiff",sep=""), width=84, height=84, units="mm", res=600,compression = "lzw")
par(oma = c(0,0,0,0), mar = c(3,3.5,1,0.5),mgp=c(1.9,0.7,0),adj=0.5,cex=0.5)
VolcanoPlot(log(foldchange,2),p.values,interactiveplot = FALSE,plimit=0.05, main=NULL,coeflimit=0.5,cexcutoff=0.4,ylab=expression(paste("-log"[10],italic(" p")," value",sep="")),xlab=expression(paste("log"[2]," fold change",sep="")),labelsig=FALSE)
dev.off()


#create HeatMap
tiff("HeatMap.tiff", width=174, height=174, units="mm", res=600, compression = "lzw")
par(oma = c(1,0,0,0), mar = c(1.5,1.5,0.1,0.1),mgp=c(0,0.7,0),cex=0.5)
HeatMap(data.raw, Class, interactiveplot = FALSE, colramp = c(25,"#fbf3f3","#cc1919"), scale="row", dendrogram="both",
        distmethod="manhattan", aggmethod = "ward.D", cexRow = 0.5,cexCol = 0.5, keysize=1)
dev.off()

data.raw.sig<-data.raw[,sigrow]
#create HeatMap
tiff("HeatMap_significants.tiff", width=174, height=174, units="mm", res=600, compression = "lzw")
par(oma = c(1,0,0,0), mar = c(1.5,1.5,0.1,0.1),mgp=c(0,0.7,0),cex=0.5)
HeatMap(subset(data.raw.sig,Filetype!="QC"), subset(Class,Filetype!="QC"), interactiveplot = FALSE, colramp = c(25,"#fbf3f3","#cc1919"),
        scale="row", dendrogram="both", distmethod="manhattan", aggmethod = "ward.D", cexRow = 0.5,cexCol = 0.5, keysize=1)
dev.off()


#plots for lipid elution patterns

#create table for elution pattern plots
features.all<-stat.summary[,-6:-ncol(stat.summary)]
colnames(features.all)<-c("AlignmentID","RT","mz","ID","Adduct")
features.all<-features.all[-1:-3,]

#remove Unknowns from dataset
if (sum(grep("Unknown", features.all$ID))>0)
{
  features.ID<-features.all[-c(grep("Unknown", features.all$ID)),]
} else
{features.ID<-features.all}

#remove RIKEN features from dataset
if (sum(grep("RIKEN", features.ID$ID))>0)
{
  features.ID<-features.ID[-c(grep("RIKEN", features.ID$ID)),]
}

#remove deuterated standards (Lipidomix) from dataset
if (sum(grep("\\(d", features.ID$ID))>0)
{
  features.ID<-features.ID[-c(grep("\\(d", features.ID$ID)),]
}

#save woMS2 and rearrange
woMS2<-ifelse(grepl("w/o MS2",features.ID$ID), "*", "")
features.ID$ID<-gsub("w/o MS2:","", features.ID$ID)


#define classes
features.class<-gsub("\\d+.*", "", features.ID$ID)
features.class<-gsub("-", "", features.class)
features.class<-as.factor(features.class)

#define number of Os
features.O<-grepl(";O/",features.ID$ID,fixed=TRUE)
features.O<-gsub("TRUE","O",features.O)
features.O<-gsub("FALSE","",features.O)
features.2O<-grepl(";2O",features.ID$ID,fixed=TRUE)
features.2O<-gsub("TRUE","2O",features.2O)
features.2O<-gsub("FALSE","",features.2O)
features.3O<-grepl(";3O",features.ID$ID,fixed=TRUE)
features.3O<-gsub("TRUE","3O",features.3O)
features.3O<-gsub("FALSE","",features.3O)
features.4O<-grepl(";4O",features.ID$ID,fixed=TRUE)
features.4O<-gsub("TRUE","4O",features.4O)
features.4O<-gsub("FALSE","",features.4O)
features.5O<-grepl(";5O",features.ID$ID,fixed=TRUE)
features.5O<-gsub("TRUE","5O",features.5O)
features.5O<-gsub("FALSE","",features.5O)
features.SN1<-grepl("SN1",features.ID$ID,fixed=TRUE)
features.SN1<-gsub("TRUE","SN1",features.SN1)
features.SN1<-gsub("FALSE","",features.SN1)
features.SN2<-grepl("SN2",features.ID$ID,fixed=TRUE)
features.SN2<-gsub("TRUE","SN2",features.SN2)
features.SN2<-gsub("FALSE","",features.SN2)
features.sideinfo<-paste(features.O,features.2O,features.3O,features.4O,features.5O,features.SN1,features.SN2,sep="")

#define species
features.species<-gsub(";.*","",features.ID$ID)
features.species<-gsub("\\|.*","",features.species)
features.species<-gsub("/","_",features.species)
features.species<-gsub("SN1","",features.species)
features.species<-gsub("SN2","",features.species)
features.species<-gsub("[^_/:0-9]","",features.species)

#special form for summing (MS-DIAL inconsistency)
features.special<-str_split_fixed(features.species, "_",3)
features.special.carbons<-gsub(":.*","",features.special)
features.special.carbons<-apply(features.special.carbons,2,as.numeric)
features.special.saturations<-gsub(".*:","",features.special)
features.special.saturations<-apply(features.special.saturations,2,as.numeric)

#sum carbons & saturations
features.carbons<-rowSums(features.special.carbons,na.rm=TRUE)
features.saturations<-as.factor(rowSums(features.special.saturations,na.rm=TRUE))

#create summary table and lipid class_adduct list  
features.ID<-cbind(features.ID,woMS2,class=features.class,species=features.species,carbons=features.carbons,saturations=features.saturations,sideinfo=features.sideinfo,class_adduct=paste(features.class,features.ID$Adduct,features.sideinfo))
features.ID$RT<-as.numeric(as.character(features.ID$RT))
features.ID$mz<-as.numeric(as.character((features.ID$mz)))
lipid_list<-split(features.ID, features.ID$class_adduct)
features.ID$class_adduct<-as.factor(features.ID$class_adduct)

#create plots
dir.create("./patterns")
setwd("./patterns")
for (i in 1:length(levels(features.ID$class_adduct)))
{
  #create plot only for classes with more than 1 identified species
  if (length(rownames(lipid_list[[i]]))>1)
  {
    ggplot(lipid_list[[i]], aes(RT, mz, colour=saturations, group=saturations)) +
      geom_point()+
      geom_line()+
      geom_text(aes(label=paste(species,woMS2),hjust=0.25,vjust=1.4),size=2)+
      scale_color_manual(values=c(pal_npg()(10),pal_d3()(length(levels(lipid_list[[i]]$saturations))-10)))+
      theme_bw()+
      geom_line(aes(group=carbons),color="black",linetype="dashed",size=0.2)+
      labs(x="rt [min]", title=levels(features.ID$class_adduct)[i])+
      ggsave(paste(levels(features.ID$class_adduct)[i],".tiff",sep=""), width=174, height=100, units="mm", dpi =600, compression="lzw")
  }}

setwd("..")
