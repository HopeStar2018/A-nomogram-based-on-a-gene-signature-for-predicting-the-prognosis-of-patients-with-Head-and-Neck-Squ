rm(list = ls())
setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO")
library(survival)
library(survivalROC)
library(penalized)
library(org.Hs.eg.db)
library(xlsx)
library(pROC)

#### Read TCGA normalize Larynx data ######
#### level 3 raw data ######
# k06 <- read.table("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\rawdata\\GEOdata.txt",header = T,row.names = 1)
#### normalized data ###################
k06 <- read.table(file = 'D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\analysis\\HNSCC_mRNA_normalizeExp.txt', header = T, row.names = 1, com = '', sep = "\t",check.names = F)
####################

k06.bak <- k06 # 存入k06.ori备用？？

###### get diffmRNA matrix ######
diffmRNA_name <- read.table("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\TCGA\\step1_TCGA_diffmRNA\\diff_mRNA_names.txt")
diffmRNA_name <- unlist(diffmRNA_name)
diffmRNA_name <- as.character(diffmRNA_name)

k06 <- t(k06)
k06 <- as.data.frame(k06)

GEO_gene_name <- colnames(k06)
over_temp_name <- intersect(GEO_gene_name,diffmRNA_name)

k06 <- subset(k06,select=over_temp_name)
k06<-t(k06)
k06 <- as.data.frame(k06)
#### load TCGA Larynx clinical data #######
## No.1 load HNSCC clinical data ####
samples <- read.xlsx2(file = "D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\analysis\\HNSCC_clinical.xlsx",sheetIndex = 1,header = T,row.names = 1)
# samples<-read.table("clipboard",header=T,sep='',row.names = 1)

k06 <- t(k06)
k06 <- as.data.frame(k06)
k06$Samples <- row.names(k06)
samples$Samples <- row.names(samples)
HNSCC_merge <- merge(samples,k06,by="Samples")


#### split 80% and 20% ####
# set.seed(2018-8-22)
# ind <- sample(2,nrow(HNSCC_merge),replace = TRUE, prob = c(0.8,0.2))
# trainData <- HNSCC_merge[ind == 1,]  # 训练数据
# proData <- HNSCC_merge[ind == 2,]  # 测试数据

#### start lasso ####
#install.packages("glmnet")
#install.packages("survival")
library("glmnet")
library("survival")

#### clinical data for proData #####
mysurv <- subset(HNSCC_merge,select = c("fustatus","futime"))
rownames(mysurv) <- HNSCC_merge[,1]
mysurv$fustatus <- as.numeric(as.character(mysurv$fustatus))
mysurv$futime <- as.numeric(as.character(mysurv$futime))
#### gene exp data for proData ###
myexpr <- HNSCC_merge[,-1:-12]
rownames(myexpr) <- HNSCC_merge[,1]

myexpr <- as.matrix(myexpr)
myexpr <- t(myexpr)
# 算出lambda值
cvfit = cv.glmnet(t(myexpr), Surv(mysurv$futime,mysurv$fustatus), 
                  #10倍交叉验证，非必须限定条件，这篇文献有，其他文献大多没提
                  nfold=10,
                  family = "cox"
) 
plot(cvfit)

#两个lambda值均可采用，具体lambda选值要根据自己实验设计而定。
#此处使用`lambda min`
cvfit$lambda.min #最佳lambda值
cvfit$lambda.1se #一倍SE内的更简洁的模型


fit <- glmnet(t(myexpr), Surv(mysurv$futime,mysurv$fustatus), family = "cox") 

#用包自带的函数画图
plot(fit, label = TRUE)

#自定义颜色
mycol <- rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D",
               "#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),2)

#设置x轴最大值
xmax <- 1300

plotCoef_plus <- function (beta, norm, lambda, df, dev, label = FALSE, legend = FALSE, xvar = c("norm", 
                                                                                                "lambda", "dev"), xlab = iname, ylab = "Coefficients", ...) 
{
  which = nonzeroCoef(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  beta = as.matrix(beta[which, , drop = FALSE])
  xvar = match.arg(xvar)
  switch(xvar, norm = {
    index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
    iname = "L1 Norm"
    approx.f = 1
  }, lambda = {
    index = log(lambda)
    iname = "Log Lambda"
    approx.f = 0
  }, dev = {
    index = dev
    iname = "Fraction Deviance Explained"
    approx.f = 1
  })
  dotlist = list(...)
  type = dotlist$type
  
  if (legend){
    #在右侧留出画图例的地方
    par(xpd = T, mar = par()$mar + c(0,0,0,6))
  }
  
  #修改bty，换个更好看的边框，还可以改成，o / n / 7 / l / c / u / ]
  if (is.null(type)) 
    matplot(index, t(beta), lty = 1, lwd = 2,
            xlab = xlab, ylab = ylab, 
            xlim = c(0, xmax), #设置x轴最大值
            col = mycol,#线的颜色
            type = "l", cex.lab=1.2, cex.axis=1,
            bty="n", ...)#不画右边框
  else matplot(index, t(beta), lty = 1, lwd = 2,
               xlab = xlab, ylab = ylab, 
               xlim = c(0, xmax), 
               col = mycol,
               type = "l", cex.lab=1.2, cex.axis=1,
               bty="n", ...)
  atdf = pretty(index)
  prettydf = approx(x = index, y = df, xout = atdf, rule = 2, 
                    method = "constant", f = approx.f)$y
  axis(3, at = atdf, labels = prettydf, tcl = NA)
  
  if (label) {
    nnz = length(which)
    xpos = max(index)
    pos = 4
    if (xvar == "lambda") {
      xpos = min(index)
      pos = 2
    }
    xpos = rep(xpos, nnz)
    ypos = beta[, ncol(beta)]
    
    #原函数打印序号，修改为打印基因名
    text(xpos, ypos, paste(rownames(myexpr)[which]),
         cex = 0.8, #基因名字体大小
         #基因名的颜色跟线一样
         col = mycol,
         #如果你不想要彩色的字，就用下面这行
         #col = "black",
         pos = pos)
  }
  if (legend) {
    #画图例
    legend("topright",
           inset=c(-0.12,0),#图例画到图外面
           legend = rownames(myexpr), #图例文字
           col = mycol, #图例线的颜色，与文字对应
           lwd = 3, #图例中线的粗细
           cex = 1, #图例字体大小
           bty = "n") #不显示图例边框
  }
  par(xpd=FALSE)
}

plot.glmnet_plus <- function (x, xvar = c("norm", "lambda", "dev"), label = FALSE, legend = FALSE,
                              ...) 
{
  xvar = match.arg(xvar)
  plotCoef_plus(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio, 
                label = label, legend = legend, xvar = xvar, ...)
}


plot.glmnet_plus(fit, label = TRUE, #打印基因名
                 legend = FALSE) #不显示图例

#在图上画虚线
#你想用哪个cutoff，就在“v = ”写上相应的数字
#此处以lambda.min作为cutoff
abline(v = cvfit$lambda.min, lty = 3, #线的类型，可以改成0, 1, 2, 3, 4, 5, 6
       lwd = 2, #线的粗细
       col = "black") #线的颜色


pdf("lasso.pdf",width = 10,height = 8)
plot.glmnet_plus(fit, label = T, #不打印基因名
                 legend = F) #显示图例
abline(v = cvfit$lambda.min, lty=3, 
       lwd = 2, 
       col = "black")
dev.off()

# 输出选中的基因
#在参数“s =”后面写cutoff
#此处选用lambda.min
coef.min = coef(cvfit, s = cvfit$lambda.min) 
coef.min

#提取选中的基因名
active.min = which(coef.min != 0)
geneids <- rownames(myexpr)[active.min]
geneids

#提取选中的基因对应的coefficient
index.min = coef.min[active.min]
index.min

#输出到文件
combine <- cbind(geneids, index.min)
write.csv(combine,"gene_index.csv")

signature <- as.matrix(t(myexpr[geneids,])) %*% as.matrix(index.min)
signature <- as.data.frame(signature)
signature$Samples <- row.names(signature)
colnames(signature)[1] <- "Score"

HNSCC_clinial <- HNSCC_merge[,1:12]
HNSCC_KM <- merge(HNSCC_clinial,signature,by="Samples")


###### KM survival analysis #######
HNSCC_KM$fustatus <- as.numeric(as.character(HNSCC_KM$fustatus))
HNSCC_KM$futime <- as.numeric(as.character(HNSCC_KM$futime))
HNSCC_KM$Score <- as.numeric(as.character(HNSCC_KM$Score))

MyAUC <- function(c, x){
  ##c:death or live
  ##x: X factor
  n <- length(c)
  n1 <- length(which(c == 1))
  n0 <- n - n1
  S1 <- cbind(x[which(c == 1)], 1)
  S0 <- cbind(x[which(c == 0)], 1)
  
  U <- 0
  for( i in 1 : n1){
    for(j in 1 : n0){
      if(S1[i, 1] > S0[j, 1]) {U <- U + 1}
      if(S1[i, 1] == S0[j, 1]) {U <- U + .5}
    }
  }
  auc <- U/(n1 * n0)
  
  vpj <- vector(length=n1)
  vqk <- vector(length=n0)
  U1 <- 0
  for(i in 1 : n1){
    for(j in 1 : n0){
      if(S1[i, 1] > S0[j, 1]) {U1 <- U1 + 1}
      else {if(S1[i, 1] == S0[j, 1]) {U1 <- U1 + 0.50 }}
    }
    division <- U1 / n0
    resta <- (division - auc)^2
    vpj[i] <- resta
    U1 <- 0
  }
  U2 <- 0
  resta <- 0
  for(j in 1 : n0){
    for(i in 1 : n1){
      if(S1[i, 1] > S0[j, 1]) {U2 <- U2 + 1}
      else {if(S1[i, 1] == S0[j, 1]) {U2 <- U2 + 0.50 }}
    }
    division <- U2 / n1
    resta <- (division - auc)^2
    vqk[j] <- resta
    U2 <- 0
  }
  vpj <- vpj/(n1 * (n1 - 1))
  vqk <- vqk/(n0 * (n0 - 1))
  var <- sum(vpj) + sum(vqk)
  s <- sqrt(var)
  p <- 1 - pnorm((auc - 0.5) / s)
  return(c(AUC = auc, lower.95 = qnorm(0.025, mean = auc, sd = s), upper.95 = qnorm(0.975, mean = auc, sd = s), p = p))
}

i = 13
rt <- HNSCC_KM

outTab=data.frame()
picDir="picture"
dir.create(picDir)

cox <- coxph(Surv(futime, fustatus) ~ rt[,i], data = rt)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]



state_vector <- rt$fustatus
marker_vector <- rt[,i]
roc <- rocdemo.sca( truth=state_vector, data=marker_vector, rule=dxrule.sca)
auc_value <- AUC(roc)
result <- MyAUC(state_vector, marker_vector)
p_value = result[['p']]

if (auc_value < 0.50){
  auc_value = 1 - auc_value
  p_value = 1 - p_value
}
res1 = c(i, auc_value, p_value)
target_pred <- marker_vector
target_class <- state_vector 




pred <- prediction(target_pred, target_class)

perf = performance(pred, "acc")

perf1 = performance(pred, "tpr", "fpr")

perf2 = performance(pred, "ppv")
perf3 = performance(pred, "npv")

cutoff.list.acc <- unlist(perf@x.values[[1]])
acc.list <- unlist(perf@y.values[[1]])

optimal_index = which.max(perf@y.values[[1]])



fp_rate = performance(pred, "fpr")
tp_rate = performance(pred, "tpr")


fpr = as.numeric(fp_rate@y.values[[1]])
tpr = as.numeric(tp_rate@y.values[[1]])

dist = rep(0, length(fpr))
dist = as.numeric(dist)

dist_index = length(fpr)

for (k in 1:dist_index){
  
  dist[[k]] <- sqrt(((fpr[[k]])^2) + (((tpr[[k]])-1)^2))
  
  
}

best_cutoff_index = which.min(dist)



fpr_opt = fp_rate@y.values[[1]][best_cutoff_index]
tpr_opt = tp_rate@y.values[[1]][best_cutoff_index]

fpr_opt = format(fpr_opt, digits=2)
tpr_opt = format(tpr_opt, digits=2)

cutoff <- cutoff.list.acc[best_cutoff_index]
med <- cutoff
# med=median(rt[,i])
if(med!=0){
  a=rt[,i]>med
  rt1=rt[a,]
  b=setdiff(rownames(rt),rownames(rt1))
  rt2=rt[b,]
  n1=nrow(rt1)
  n2=nrow(rt2)
  surTab1=summary(survfit(Surv(futime, fustatus) ~ 1, data = rt1))
  surTab2=summary(survfit(Surv(futime, fustatus) ~ 1, data = rt2))
  medianTab1=surTab1$table
  medianTab2=surTab2$table
  diff=survdiff(Surv(futime, fustatus) ~a,data = rt)
  fit <- survfit(Surv(futime, fustatus) ~ a, data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=i,coxSummary$coefficients,coxSummary$conf.int,KM=pValue,
                            H_med=medianTab1["median"],H_0.95LCL=medianTab1["0.95LCL"],H_0.95UCL=medianTab1["0.95UCL"],
                            L_med=medianTab2["median"],L_0.95LCL=medianTab2["0.95LCL"],L_0.95UCL=medianTab2["0.95UCL"]))
  pval=0
  if(pValue<0.05){
    pval=signif(pValue,4)
    pval=format(pval, scientific = TRUE)
  }else{
    pval=round(pValue,3)
  }
  if(pValue<0.05){
    geneName="Score"
    tiffFile="Figure 3(A) Traning_survival.tiff"
    outTiff=paste(picDir,tiffFile,sep="\\")
    tiff(file=outTiff,width = 15,height = 15,units ="cm",compression="lzw",bg="white",res=600)
    plot(fit, col=c("blue","red"), xlab="Time (days)", ylab="Overall survival",
         main=paste(geneName,"(p=",pval, ")", sep=""),mark.time=T,ylim=c(0,1.1),xlim=c(0,1825),
         lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
    legend("topright", c(paste("Low Risk"), 
                         paste("High Risk")), 
           col=c("blue","red"), bty="n", lwd = 2, cex=0.8)
    dev.off()
  }
}

########### Primary survival analysis method II ##################################

library(survival)
library(ROC)
library(ROCR)
library(hash)
library(qvalue)
library(VennDiagram)
library(xlsx)
library(survMisc)
library(pROC)
library(ggplot2)
library(survminer)

merge_HNSCC_validation <- HNSCC_KM

###### KM survival analysis #######
merge_HNSCC_validation$fustatus <- as.numeric(as.character(merge_HNSCC_validation$fustatus))
merge_HNSCC_validation$futime <- as.numeric(as.character(merge_HNSCC_validation$futime))
merge_HNSCC_validation$Score <- as.numeric(as.character(merge_HNSCC_validation$Score))


roc1 <- roc(merge_HNSCC_validation$fustatus,merge_HNSCC_validation$Score)
# plot(roc1,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"), max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)
cutoff <- roc1$thresholds[which.max(roc1$specificities+roc1$sensitivities-1)]

merge_HNSCC_validation$Score[which(merge_HNSCC_validation$Score < cutoff)] <- 0
merge_HNSCC_validation$Score[which(merge_HNSCC_validation$Score >= cutoff)] <- 1
merge_HNSCC_validation$Score <- factor(merge_HNSCC_validation$Score,labels = c("Low-Risk","High-Risk"))


fit <- survfit(Surv(futime, fustatus)~Score, data=merge_HNSCC_validation)


outTab=data.frame()
picDir="picture"
dir.create(picDir)

tiffFile="Figure 3(A) Primary_survival.tiff"
outTiff=paste(picDir,tiffFile,sep="\\")
tiff(file=outTiff,width = 20,height = 15,units ="cm",compression="lzw",bg="white",res=600)
# ggsurvplot(fit, conf.int=F,risk.table=T, risk.table.col="strata", pval=T)
ggsurvplot(fit, conf.int=F, pval=T, risk.table=T, legend.labs = c("Low Risk","High Risk"), 
           legend.title="PGS Score", palette = c("dodgerblue2","orchid2"),
           risk.table.height = 0.3)


dev.off()





####################  save best matrix for next step #######################
write.csv(HNSCC_KM,"result_step1_HNSCC.csv")



########################################################################
################### step2 analysis #####################################
########################################################################

rm(list = ls())

library(rms)
library(Hmisc)
library(lattice)
library(survival)
library(Formula)
library(ggplot2)
library(foreign)
library(rms)
library(xlsx)
library(survMisc)

setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step2_modeldisscusion")
HNSCC_nomogram <- read.xlsx2("Training_data.xlsx",row.names = 1, header = T,sheetIndex = 1)


HNSCC_nomogram$futime <- as.numeric(as.character(HNSCC_nomogram$futime))
HNSCC_nomogram$fustatus <- as.numeric(as.character(HNSCC_nomogram$fustatus))
HNSCC_nomogram$Score <- as.numeric(as.character(HNSCC_nomogram$Score))

# HNSCC_nomogram[,c(1:7)] <- lapply(HNSCC_nomogram[,c(1:7)],as.numeric)

### as vector ####
HNSCC_nomogram$Gender <- factor(HNSCC_nomogram$Gender,labels = c("Female","male"))
HNSCC_nomogram$Smoking <- factor(HNSCC_nomogram$Smoking,labels = c("non-Smoking","Smoking"))
HNSCC_nomogram$Alcohol <- factor(HNSCC_nomogram$Alcohol,labels = c("non-Alcohol","Alcohol"))
HNSCC_nomogram$T_category <- factor(HNSCC_nomogram$T_category, labels = c("T1&T2","T3&T4"))
HNSCC_nomogram$N_category <- factor(HNSCC_nomogram$N_category, labels = c("N0&N1","N2","N3"))
HNSCC_nomogram$HPV_status <- factor(HNSCC_nomogram$HPV_status, labels = c("HPV-","HPV+"))
HNSCC_nomogram$Age <- factor(HNSCC_nomogram$Age,labels = c("Young","Middle-aged","Young-seniors","Elderly"))

############### HRs forest plot ###################################################
setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step5_cox_forestplot")


###### Training univariate ##############################
### 单因素分析 ########
library(xlsx)
library(plyr)
BaSurv <- Surv(time = HNSCC_nomogram$futime,event = HNSCC_nomogram$fustatus)

# GCox <- coxph(BaSurv ~ Gender, data = HNSCC_nomogram)
# GSum <- summary(GCox)
# HR <- round(GSum$coefficients[,2],2)
# PValue <- round(GSum$coefficients[,5],3)
# CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-")
# Unicox <- data.frame("characteristics" = "Gender",
#                      "Hazard Ratio" = HR,
#                      "95CI" = CI,
#                      "P Value" = PValue)

UniCox <- function(x){
  FML <- as.formula(paste0('BaSurv~',x))
  GCox <- coxph(FML, data = HNSCC_nomogram)
  GSum <- summary(GCox)
  HR <- round(GSum$coefficients[,2],2)
  PValue <- round(GSum$coefficients[,5],3)
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-")
  Unicox <- data.frame("characteristics" = x,
                       "Hazard Ratio" = HR,
                       "CI95" = CI,
                       "P Value" = PValue)
  return(Unicox)
}
UniCox('Gender')

VarNames <- colnames(HNSCC_nomogram)[c(1:7,10)]
UniVar <- lapply(VarNames,UniCox)
UniVar <- ldply(UniVar,data.frame)

fml <-  as.formula(paste0('BaSurv~',paste0(UniVar$characteristics[UniVar$P.Value < 0.2],collapse = '+')))
MultiCox <- coxph(fml, data = HNSCC_nomogram)
MultiSum <- summary(MultiCox)

MultiName <- as.character(UniVar$characteristics[UniVar$P.Value < 0.2])
MHR <- round(MultiSum$coefficients[,2],2)
MPV <- round(MultiSum$coefficients[,5],3)
MCIL <- round(MultiSum$conf.int[,3],2)
MCIU <- round(MultiSum$conf.int[,4],2)
MCI <- paste0(MCIL,'-',MCIU)
MulCox <- data.frame('characteristics' = MultiName,
                     'Hazard Ratio' = MHR,
                     'CI95' = MCI,
                     'P Value' = MPV)

##### merge data frame #################
Final <- merge.data.frame(UniVar,MulCox, by = "characteristics", all = T, sort = T)

write.xlsx(Final,'CoxReg.xlsx',col.names = T,row.names = F,showNA = F)

######## model discussion ###########################

temp_name <-  as.character(UniVar$characteristics[UniVar$P.Value < 0.2])
temp_name <- temp_name[-4]
fml_without_HPV <-  as.formula(paste0('BaSurv~',paste0(temp_name,collapse = '+')))
MultiCox_without_HPV <- coxph(fml_without_HPV, data = HNSCC_nomogram)
MultiSum_without_HPV <- summary(MultiCox_without_HPV)

MultiName_without_HPV <- temp_name
MHR_without_HPV <- round(MultiSum_without_HPV$coefficients[,2],2)
MPV_without_HPV <- round(MultiSum_without_HPV$coefficients[,5],3)
MCIL_without_HPV <- round(MultiSum_without_HPV$conf.int[,3],2)
MCIU_without_HPV <- round(MultiSum_without_HPV$conf.int[,4],2)
MCI_without_HPV <- paste0(MCIL_without_HPV,'-',MCIU_without_HPV)
MulCox_without_HPV<- data.frame('characteristics' = MultiName_without_HPV,
                              'Hazard Ratio' = MHR_without_HPV,
                              'CI95' = MCI_without_HPV,
                              'P Value' = MPV_without_HPV)


Model_disscusion <- merge.data.frame(MulCox_without_HPV,MulCox, by = "characteristics", all = T, sort = T)


write.xlsx(Model_disscusion,'Table 3. Model_disscusion.xlsx',col.names = T,row.names = F,showNA = F)



############ nomogram plot ########################################################
setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step2_modeldisscusion")
HNSCC_nomogram <- read.xlsx2("Training_data.xlsx",row.names = 1, header = T,sheetIndex = 1)


HNSCC_nomogram$futime <- as.numeric(as.character(HNSCC_nomogram$futime))
HNSCC_nomogram$fustatus <- as.numeric(as.character(HNSCC_nomogram$fustatus))
HNSCC_nomogram$Score <- as.numeric(as.character(HNSCC_nomogram$Score))

### as vector ####
HNSCC_nomogram$Gender <- factor(HNSCC_nomogram$Gender,labels = c("Female","male"))
HNSCC_nomogram$Smoking <- factor(HNSCC_nomogram$Smoking,labels = c("non-Smoking","Smoking"))
HNSCC_nomogram$Alcohol <- factor(HNSCC_nomogram$Alcohol,labels = c("non-Alcohol","Alcohol"))
HNSCC_nomogram$T_category <- factor(HNSCC_nomogram$T_category, labels = c("T1&T2","T3&T4"))
HNSCC_nomogram$N_category <- factor(HNSCC_nomogram$N_category, labels = c("N0&N1","N2","N3"))
HNSCC_nomogram$HPV_status <- factor(HNSCC_nomogram$HPV_status, labels = c("HPV-","HPV+"))
HNSCC_nomogram$Age <- factor(HNSCC_nomogram$Age,labels = c("Young","Middle-aged","Young-seniors","Elderly"))

bc <- na.omit(HNSCC_nomogram)  #函数na.omit()即删除数据框中所有含缺失值的观测

dd <- datadist(bc)
options(datadist="dd")

f <- cph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_stutas + Score, x=T, y=T, surv=T, data=bc)
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(365, x), function(x) surv(1095, x), function(x) surv(1825, x)), lp=F, funlabel=c("1-year survival","3-year survival", "5-year survival"), maxscale=10, fun.at=seq(0.1,0.9,0.1))   #maxscale为列线图第一行最大的分值，默认值100，这是文献中列线图普遍采用的最大分值；本例由于原文设定最大分值为10，故输入代码maxscale=10
plot(nom)

validate(f, method="boot", B=1000, dxy=T)
rcorrcens(Surv(futime, fustatus) ~ predict(f), data = bc)

f1 <- cph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_status + Score, x=T, y=T, surv=T, data=bc, time.inc=365)
cal1 <- calibrate(f1, cmethod="KM", method="boot", u=365, m=60, B=1000)
plot(cal1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlab="Nomogram-Predicted Probability of 1-Year OS in training cohort",
     ylab="Actual 1-Year OS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))

# f3 <- cph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_status + Score, x=T, y=T, surv=T, data=bc, time.inc=1095)
# cal3 <- calibrate(f3, cmethod="KM", method="boot", u=1095, m=60, B=1000)
# plot(cal3,lwd=2,lty=1,
#      errbar.col=c(rgb(0,118,192,maxColorValue=255)),
#      xlab="Nomogram-Predicted Probability of 3-Year OS",
#      ylab="Actual 3-Year OS (proportion)",
#      col=c(rgb(192,98,83,maxColorValue=255)))
# 
# 
# f5 <- cph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_status + Score, x=T, y=T, surv=T, data=bc, time.inc=1825)
# cal5 <- calibrate(f5, cmethod="KM", method="boot", u=1825, m=60, B=1000)
# plot(cal5,lwd=2,lty=1,
#      errbar.col=c(rgb(0,118,192,maxColorValue=255)),
#      xlab="Nomogram-Predicted Probability of 5-Year OS",
#      ylab="Actual 5-Year OS (proportion)",
#      col=c(rgb(192,98,83,maxColorValue=255)))

c1 <- coxph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_status + Score, data=bc)
gof(c1)


##### save result plot #######
setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step2_modeldisscusion\\Training_result")
#### plot 1 nomogram #######
pdf("Nomogram_Training.pdf",width = 12,height = 10)
plot(nom)
dev.off()

#### plot 2  1 year calibration curve #####
pdf("1 year calibration curve.pdf",width = 12,height = 10)
plot(cal1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlab="Nomogram-Predicted Probability of 1-Year OS",
     ylab="Actual 1-Year OS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))
dev.off()


# #### plot 3  3 year calibration curve #####
# pdf("3 year calibration curve.pdf",width = 12,height = 10)
# plot(cal3,lwd=2,lty=1,
#      errbar.col=c(rgb(0,118,192,maxColorValue=255)),
#      xlab="Nomogram-Predicted Probability of 3-Year OS",
#      ylab="Actual 3-Year OS (proportion)",
#      col=c(rgb(192,98,83,maxColorValue=255)))
# dev.off()
# 
# #### plot 4  5 year calibration curve #####
# pdf("5 year calibration curve.pdf",width = 12,height = 10)
# plot(cal5,lwd=2,lty=1,
#      errbar.col=c(rgb(0,118,192,maxColorValue=255)),
#      xlab="Nomogram-Predicted Probability of 5-Year OS",
#      ylab="Actual 5-Year OS (proportion)",
#      col=c(rgb(192,98,83,maxColorValue=255)))
# dev.off()

########### Training survival analysis method II ##################################
rm(list = ls())
library(survival)
library(ROC)
library(ROCR)
library(hash)
library(qvalue)
library(VennDiagram)
library(xlsx)
library(survMisc)
library(survminer)

setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step2_modeldisscusion")
HNSCC_nomogram <- read.xlsx2("Training_data.xlsx",row.names = 1, header = T,sheetIndex = 1)

merge_HNSCC_validation <- HNSCC_nomogram

###### KM survival analysis #######
merge_HNSCC_validation$fustatus <- as.numeric(as.character(merge_HNSCC_validation$fustatus))
merge_HNSCC_validation$futime <- as.numeric(as.character(merge_HNSCC_validation$futime))
merge_HNSCC_validation$Score <- as.numeric(as.character(merge_HNSCC_validation$Score))


roc1 <- roc(merge_HNSCC_validation$fustatus,merge_HNSCC_validation$Score)
# plot(roc1,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"), max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)
cutoff <- roc1$thresholds[which.max(roc1$specificities+roc1$sensitivities-1)]

merge_HNSCC_validation$Score[which(merge_HNSCC_validation$Score < cutoff)] <- 0
merge_HNSCC_validation$Score[which(merge_HNSCC_validation$Score >= cutoff)] <- 1
merge_HNSCC_validation$Score <- factor(merge_HNSCC_validation$Score,labels = c("Low-Risk","High-Risk"))


fit <- survfit(Surv(futime, fustatus)~Score, data=merge_HNSCC_validation)


outTab=data.frame()
picDir="picture"
dir.create(picDir)

tiffFile="Figure 3(B) Trainning_survival.tiff"
outTiff=paste(picDir,tiffFile,sep="\\")
tiff(file=outTiff,width = 20,height = 15,units ="cm",compression="lzw",bg="white",res=600)
# ggsurvplot(fit, conf.int=F,risk.table=T, risk.table.col="strata", pval=T)
ggsurvplot(fit, conf.int=F, pval=T, risk.table=T, legend.labs = c("Low Risk","High Risk"), 
           legend.title="PGS Score", palette = c("dodgerblue2","orchid2"),
           risk.table.height = 0.3)
dev.off()


########################### step 3 inter-validation GEO #####################################
rm(list = ls())
library(xlsx)
setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step3_inter_Validation_GEO")
validation_data <- read.xlsx2("Validation_data.xlsx",sheetIndex = 1,row.names =1)
 

validation_data$futime <- as.numeric(as.character(validation_data$futime))
validation_data$fustatus <- as.numeric(as.character(validation_data$fustatus))
validation_data$Score <- as.numeric(as.character(validation_data$Score))

#### as vector ####
validation_data$Gender <- factor(validation_data$Gender,labels = c("Female","male"))
validation_data$Smoking <- factor(validation_data$Smoking,labels = c("non-Smoking","Smoking"))
validation_data$Alcohol <- factor(validation_data$Alcohol,labels = c("non-Alcohol","Alcohol"))
validation_data$T_category <- factor(validation_data$T_category, labels = c("T1&T2","T3&T4"))
validation_data$N_category <- factor(validation_data$N_category, labels = c("N0&N1","N2","N3"))
validation_data$HPV_status <- factor(validation_data$HPV_status, labels = c("HPV-","HPV+"))
validation_data$Age <- factor(validation_data$Age,labels = c("Young","Middle-aged","Young-seniors","Elderly"))

bc <- na.omit(validation_data)  #函数na.omit()即删除数据框中所有含缺失值的观测

dd <- datadist(bc)
options(datadist="dd")

f <- cph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_status + Score, x=T, y=T, surv=T, data=bc)
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(365, x), function(x) surv(1095, x), function(x) surv(1825, x)), lp=F, funlabel=c("1-year survival","3-year survival", "5-year survival"), maxscale=10, fun.at=seq(0.1,0.9,0.1))   #maxscale为列线图第一行最大的分值，默认值100，这是文献中列线图普遍采用的最大分值；本例由于原文设定最大分值为10，故输入代码maxscale=10

validate(f, method="boot", B=1000, dxy=T)
rcorrcens(Surv(futime, fustatus) ~ predict(f), data = bc)

f1 <- cph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_status + Score, x=T, y=T, surv=T, data=bc, time.inc=365)
cal1 <- calibrate(f1, cmethod="KM", method="boot", u=365, m=20, B=1000)
plot(cal1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlab="Nomogram-Predicted Probability of 1-Year OS in internal cohort",
     ylab="Actual 1-Year OS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))


##### internal validation cohort overall survival ##########################
rm(list = ls())
library(survival)
library(ROC)
library(ROCR)
library(hash)
library(qvalue)
library(VennDiagram)
library(xlsx)
setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step3_inter_Validation_GEO")
validation_data <- read.xlsx2("Validation_data.xlsx",sheetIndex = 1,row.names =1)

###### KM survival analysis #######
validation_data$fustatus <- as.numeric(as.character(validation_data$fustatus))
validation_data$futime <- as.numeric(as.character(validation_data$futime))
validation_data$Score <- as.numeric(as.character(validation_data$Score))

MyAUC <- function(c, x){
  ##c:death or live
  ##x: X factor
  n <- length(c)
  n1 <- length(which(c == 1))
  n0 <- n - n1
  S1 <- cbind(x[which(c == 1)], 1)
  S0 <- cbind(x[which(c == 0)], 1)
  
  U <- 0
  for( i in 1 : n1){
    for(j in 1 : n0){
      if(S1[i, 1] > S0[j, 1]) {U <- U + 1}
      if(S1[i, 1] == S0[j, 1]) {U <- U + .5}
    }
  }
  auc <- U/(n1 * n0)
  
  vpj <- vector(length=n1)
  vqk <- vector(length=n0)
  U1 <- 0
  for(i in 1 : n1){
    for(j in 1 : n0){
      if(S1[i, 1] > S0[j, 1]) {U1 <- U1 + 1}
      else {if(S1[i, 1] == S0[j, 1]) {U1 <- U1 + 0.50 }}
    }
    division <- U1 / n0
    resta <- (division - auc)^2
    vpj[i] <- resta
    U1 <- 0
  }
  U2 <- 0
  resta <- 0
  for(j in 1 : n0){
    for(i in 1 : n1){
      if(S1[i, 1] > S0[j, 1]) {U2 <- U2 + 1}
      else {if(S1[i, 1] == S0[j, 1]) {U2 <- U2 + 0.50 }}
    }
    division <- U2 / n1
    resta <- (division - auc)^2
    vqk[j] <- resta
    U2 <- 0
  }
  vpj <- vpj/(n1 * (n1 - 1))
  vqk <- vqk/(n0 * (n0 - 1))
  var <- sum(vpj) + sum(vqk)
  s <- sqrt(var)
  p <- 1 - pnorm((auc - 0.5) / s)
  return(c(AUC = auc, lower.95 = qnorm(0.025, mean = auc, sd = s), upper.95 = qnorm(0.975, mean = auc, sd = s), p = p))
}

i = 10

outTab=data.frame()
picDir="picture"
dir.create(picDir)

cox <- coxph(Surv(futime, fustatus) ~ validation_data[,i], data = validation_data)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]

rt <- validation_data

state_vector <- rt$fustatus
marker_vector <- rt[,i]
roc <- rocdemo.sca( truth=state_vector, data=marker_vector, rule=dxrule.sca)
auc_value <- AUC(roc)
result <- MyAUC(state_vector, marker_vector)
p_value = result[['p']]

if (auc_value < 0.50){
  auc_value = 1 - auc_value
  p_value = 1 - p_value
}
res1 = c(i, auc_value, p_value)
target_pred <- marker_vector
target_class <- state_vector 




pred <- prediction(target_pred, target_class)

perf = performance(pred, "acc")

perf1 = performance(pred, "tpr", "fpr")

perf2 = performance(pred, "ppv")
perf3 = performance(pred, "npv")

cutoff.list.acc <- unlist(perf@x.values[[1]])
acc.list <- unlist(perf@y.values[[1]])

optimal_index = which.max(perf@y.values[[1]])



fp_rate = performance(pred, "fpr")
tp_rate = performance(pred, "tpr")


fpr = as.numeric(fp_rate@y.values[[1]])
tpr = as.numeric(tp_rate@y.values[[1]])

dist = rep(0, length(fpr))
dist = as.numeric(dist)

dist_index = length(fpr)

for (k in 1:dist_index){
  
  dist[[k]] <- sqrt(((fpr[[k]])^2) + (((tpr[[k]])-1)^2))
  
  
}

best_cutoff_index = which.min(dist)



fpr_opt = fp_rate@y.values[[1]][best_cutoff_index]
tpr_opt = tp_rate@y.values[[1]][best_cutoff_index]

fpr_opt = format(fpr_opt, digits=2)
tpr_opt = format(tpr_opt, digits=2)

cutoff <- cutoff.list.acc[best_cutoff_index]
med <- cutoff
# med=median(rt[,i])
if(med!=0){
  a=rt[,i]>med
  rt1=rt[a,]
  b=setdiff(rownames(rt),rownames(rt1))
  rt2=rt[b,]
  n1=nrow(rt1)
  n2=nrow(rt2)
  surTab1=summary(survfit(Surv(futime, fustatus) ~ 1, data = rt1))
  surTab2=summary(survfit(Surv(futime, fustatus) ~ 1, data = rt2))
  medianTab1=surTab1$table
  medianTab2=surTab2$table
  diff=survdiff(Surv(futime, fustatus) ~a,data = rt)
  fit <- survfit(Surv(futime, fustatus) ~ a, data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=i,coxSummary$coefficients,coxSummary$conf.int,KM=pValue,
                            H_med=medianTab1["median"],H_0.95LCL=medianTab1["0.95LCL"],H_0.95UCL=medianTab1["0.95UCL"],
                            L_med=medianTab2["median"],L_0.95LCL=medianTab2["0.95LCL"],L_0.95UCL=medianTab2["0.95UCL"]))
  pval=0
  if(pValue<0.05){
    pval=signif(pValue,4)
    pval=format(pval, scientific = TRUE)
  }else{
    pval=round(pValue,3)
  }
  if(pValue<0.05){
    geneName="Score"
    tiffFile="Figure 3(B) Internal_validation_survival.tiff"
    outTiff=paste(picDir,tiffFile,sep="\\")
    tiff(file=outTiff,width = 15,height = 15,units ="cm",compression="lzw",bg="white",res=600)
    plot(fit, col=c("blue","red"), xlab="Time (days)", ylab="Overall survival",
         main=paste(geneName,"(p=",pval, ")", sep=""),mark.time=T,ylim=c(0,1.1),xlim=c(0,1825),
         lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
    legend("topright", c(paste("Low Risk"), 
                         paste("High Risk")), 
           col=c("blue","red"), bty="n", lwd = 2, cex=0.8)
    dev.off()
  }
}

########### internal validation survival analysis method II ##################################
rm(list = ls())
library(survival)
library(ROC)
library(ROCR)
library(hash)
library(qvalue)
library(VennDiagram)
library(xlsx)
library(survMisc)
library(survminer)

setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step3_inter_Validation_GEO")
validation_data <- read.xlsx2("Validation_data.xlsx",sheetIndex = 1,row.names =1)

merge_HNSCC_validation <- validation_data

###### KM survival analysis #######
merge_HNSCC_validation$fustatus <- as.numeric(as.character(merge_HNSCC_validation$fustatus))
merge_HNSCC_validation$futime <- as.numeric(as.character(merge_HNSCC_validation$futime))
merge_HNSCC_validation$Score <- as.numeric(as.character(merge_HNSCC_validation$Score))


roc1 <- roc(merge_HNSCC_validation$fustatus,merge_HNSCC_validation$Score)
# plot(roc1,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"), max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)
cutoff <- roc1$thresholds[which.max(roc1$specificities+roc1$sensitivities-1)]

merge_HNSCC_validation$Score[which(merge_HNSCC_validation$Score < cutoff)] <- 0
merge_HNSCC_validation$Score[which(merge_HNSCC_validation$Score >= cutoff)] <- 1
merge_HNSCC_validation$Score <- factor(merge_HNSCC_validation$Score,labels = c("Low-Risk","High-Risk"))


fit <- survfit(Surv(futime, fustatus)~Score, data=merge_HNSCC_validation)


outTab=data.frame()
picDir="picture"
dir.create(picDir)

tiffFile="Figure 3(B) internal_validation_survival.tiff"
outTiff=paste(picDir,tiffFile,sep="\\")
tiff(file=outTiff,width = 20,height = 15,units ="cm",compression="lzw",bg="white",res=600)
# ggsurvplot(fit, conf.int=F,risk.table=T, risk.table.col="strata", pval=T)
ggsurvplot(fit, conf.int=F, pval=T, risk.table=T, legend.labs = c("Low Risk","High Risk"), 
           legend.title="PGS Score", palette = c("dodgerblue2","orchid2"),
           risk.table.height = 0.3)
dev.off()



#################### step 4 external validation cohort TCGA ##################################
rm(list = ls())
library(xlsx)
setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\TCGA\\step2_getMatrix")
HNSCC_clinical_TCGA <- read.xlsx2("External_validation_TCGA_HNSCC.xlsx",sheetIndex = 1)
HNSCC_mRNA_TCGA <- read.table("HNSCC_mRNA_normalizeExp.txt",header = T,check.names = F,row.names = 1)

#### get colnames and rownames ######
colnames(HNSCC_mRNA_TCGA) <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(HNSCC_mRNA_TCGA))

HNSCC_mRNA_TCGA <- as.data.frame(t(HNSCC_mRNA_TCGA))


for(i in 1:18108){
  name <- unlist(strsplit(colnames(HNSCC_mRNA_TCGA)[i],'[|]')) 
  colnames(HNSCC_mRNA_TCGA)[i] <- name[1]
}


##### get geneid and score ###############
gene_index <- read.csv("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step1\\gene_index.csv")
gene_name <- as.character(gene_index$geneids)


HNSCC_validation_TCGA <- subset(HNSCC_mRNA_TCGA,select = gene_name)
# HNSCC_validation_TCGA$Samples <- rownames(HNSCC_validation_TCGA)

signature <- as.matrix(HNSCC_validation_TCGA) %*% as.matrix(gene_index$index.min)
signature <- as.data.frame(signature)
signature$Samples <- row.names(signature)
colnames(signature)[1] <- "Score"
colnames(HNSCC_clinical_TCGA)[1] <- "Samples"

merge_HNSCC_validation <- merge(HNSCC_clinical_TCGA,signature,by="Samples")
write.csv(merge_HNSCC_validation,"merge_HNSCC_validation.csv")

############ analysis external validation TCGA ###########################
rm(list = ls())

library(rms)
library(Hmisc)
library(lattice)
library(survival)
library(Formula)
library(ggplot2)
library(foreign)
library(rms)
library(xlsx)
library(survMisc)

setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\TCGA\\step2_getMatrix")
merge_HNSCC_validation <- read.csv("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\TCGA\\step2_getMatrix\\merge_HNSCC_validation.csv",row.names = 1)

merge_HNSCC_validation$futime <- as.numeric(as.character(merge_HNSCC_validation$futime))
merge_HNSCC_validation$fustatus <- as.numeric(as.character(merge_HNSCC_validation$fustatus))
merge_HNSCC_validation$Score <- as.numeric(as.character(merge_HNSCC_validation$Score))

################ Socre to small ######################
# merge_HNSCC_validation$Score <- log10(merge_HNSCC_validation$Score)
#### as vector ####

merge_HNSCC_validation$T_category <- factor(merge_HNSCC_validation$T_category, labels = c("T1&T2","T3&T4"))
merge_HNSCC_validation$N_category <- factor(merge_HNSCC_validation$N_category, labels = c("N0&N1","N2","N3"))
merge_HNSCC_validation$HPV_status <- factor(merge_HNSCC_validation$HPV_status, labels = c("HPV-","HPV+"))
merge_HNSCC_validation$Age <- factor(merge_HNSCC_validation$Age,labels = c("Young","Middle-aged","Young-seniors","Elderly"))



bc <- na.omit(merge_HNSCC_validation)  #函数na.omit()即删除数据框中所有含缺失值的观测

dd <- datadist(bc)
options(datadist="dd")

f <- cph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_status + Score, x=T, y=T, surv=T, data=bc)
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(365, x), function(x) surv(1095, x)), lp=F, funlabel=c("1-year survival","3-year survival"), maxscale=10, fun.at=seq(0.1,0.9,0.1))   #maxscale为列线图第一行最大的分值，默认值100，这是文献中列线图普遍采用的最大分值；本例由于原文设定最大分值为10，故输入代码maxscale=10

validate(f, method="boot", B=1000, dxy=T)
rcorrcens(Surv(futime, fustatus) ~ predict(f), data = bc)

f1 <- cph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_status + Score, x=T, y=T, surv=T, data=bc, time.inc=365)
cal1 <- calibrate(f1, cmethod="KM", method="boot", u=365, m=10, B=1000)
plot(cal1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlab="Nomogram-Predicted Probability of 1-Year OS in External cohort",
     ylab="Actual 1-Year OS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))

# f3 <- cph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_status + Score, x=T, y=T, surv=T, data=bc, time.inc=1095)
# cal3 <- calibrate(f3, cmethod="KM", method="boot", u=1095, m=9)
# plot(cal3,lwd=2,lty=1,
#      errbar.col=c(rgb(0,118,192,maxColorValue=255)),
#      xlab="Nomogram-Predicted Probability of 3-Year OS",
#      ylab="Actual 3-Year OS (proportion)",
#      col=c(rgb(192,98,83,maxColorValue=255)))

c1 <- coxph(Surv(futime, fustatus) ~ Age + T_category + N_category + HPV_status + Score, data=bc)
gof(c1)

###### exteranl validation survival  ###############
rm(list = ls())

library(survival)
library(ROC)
library(ROCR)
library(hash)
library(qvalue)
library(VennDiagram)
library(xlsx)
library(ggplot2)
library(pROC)
library(survminer)
setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\TCGA\\step2_getMatrix")
merge_HNSCC_validation <- read.csv("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\TCGA\\step2_getMatrix\\merge_HNSCC_validation.csv",row.names = 1)

###### KM survival analysis #######
merge_HNSCC_validation$fustatus <- as.numeric(as.character(merge_HNSCC_validation$fustatus))
merge_HNSCC_validation$futime <- as.numeric(as.character(merge_HNSCC_validation$futime))
merge_HNSCC_validation$Score <- as.numeric(as.character(merge_HNSCC_validation$Score))


roc1 <- roc(merge_HNSCC_validation$fustatus,merge_HNSCC_validation$Score)
plot(roc1,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"), max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)
cutoff <- roc1$thresholds[which.max(roc1$specificities+roc1$sensitivities-1)]

merge_HNSCC_validation$Score[which(merge_HNSCC_validation$Score < cutoff)] <- 0
merge_HNSCC_validation$Score[which(merge_HNSCC_validation$Score >= cutoff)] <- 1
merge_HNSCC_validation$Score <- factor(merge_HNSCC_validation$Score,labels = c("Low-Risk","High-Risk"))


fit <- survfit(Surv(futime, fustatus)~Score, data=merge_HNSCC_validation)
outTab=data.frame()
picDir="picture"
dir.create(picDir)

tiffFile="Figure 3(C) external_validation_survival.tiff"
outTiff=paste(picDir,tiffFile,sep="\\")
tiff(file=outTiff,width = 20,height = 15,units ="cm",compression="lzw",bg="white",res=600)
ggsurvplot(fit, conf.int=F, pval=T, risk.table=T, legend.labs = c("Low Risk","High Risk"), 
           legend.title="PGS Score", palette = c("dodgerblue2","orchid2"),
           risk.table.height = 0.3)
dev.off()




################### step5 Traning DCA #####################################
rm(list = ls())

library(rms)
library(Hmisc)
library(lattice)
library(survival)
library(Formula)
library(ggplot2)
library(foreign)
library(rms)
library(xlsx)
library(survMisc)

setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step2_modeldisscusion")
HNSCC_nomogram <- read.xlsx2("Training_data.xlsx",row.names = 1, header = T,sheetIndex = 1)


HNSCC_nomogram$futime <- as.numeric(as.character(HNSCC_nomogram$futime))
HNSCC_nomogram$fustatus <- as.numeric(as.character(HNSCC_nomogram$fustatus))
HNSCC_nomogram$Score <- as.numeric(as.character(HNSCC_nomogram$Score))

#### as vector ####
HNSCC_nomogram$Gender <- factor(HNSCC_nomogram$Gender,labels = c("Female","male"))
HNSCC_nomogram$Smoking <- factor(HNSCC_nomogram$Smoking,labels = c("non-Smoking","Smoking"))
HNSCC_nomogram$Alcohol <- factor(HNSCC_nomogram$Alcohol,labels = c("non-Alcohol","Alcohol"))
HNSCC_nomogram$T_category <- factor(HNSCC_nomogram$T_category, labels = c("T1&T2","T3&T4"))
HNSCC_nomogram$N_category <- factor(HNSCC_nomogram$N_category, labels = c("N0&N1","N2","N3"))
HNSCC_nomogram$HPV_status <- factor(HNSCC_nomogram$HPV_status, labels = c("HPV-","HPV+"))
HNSCC_nomogram$Age <- factor(HNSCC_nomogram$Age,labels = c("Young","Middle-aged","Young-seniors","Elderly"))

# Basic Data Set-up
#Source file to use stdca command
#Creates a survival object with time to event variable as ttcancer and the event is 
#cancer. 
setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\TCGA\\step3_DCA")
source("stdca.R")
library(survival)
Srv = Surv(HNSCC_nomogram$futime, HNSCC_nomogram$fustatus)

#Run the cox model
coxmod = coxph(Srv ~ Age + T_category + N_category + HPV_status + Score, data=HNSCC_nomogram)
#the probability of failure is calculated by subtracting the probability of 
#survival from 1. 
HNSCC_nomogram$Gene_signature_nomogram = c(1- (summary(survfit(coxmod,newdata=HNSCC_nomogram), times=365)$surv))

#Run the decision curve analysis (with a smoother)
stdca(data=HNSCC_nomogram, outcome="fustatus", ttoutcome="futime", timepoint=365, predictors="Gene_signature_nomogram", xstop=0.5, smooth=TRUE)

#Run the cox model without HPV_status
coxmod_H = coxph(Srv ~ Age + T_category + N_category + Score, data=HNSCC_nomogram)

#the probability of failure is calculated by subtracting the probability of 
#survival from 1. 
HNSCC_nomogram$Model_disintegrating_HPV  = c(1- (summary(survfit(coxmod_H,newdata=HNSCC_nomogram), times=365)$surv))



km = stdca(data=HNSCC_nomogram, outcome="fustatus", ttoutcome="futime", timepoint=365, predictors="Gene_signature_nomogram", xstop=0.5, smooth=TRUE)

#Competing Risk Model
cr = stdca(data=HNSCC_nomogram, outcome="fustatus", ttoutcome="futime", timepoint=365, predictors="Model_disintegrating_HPV", cmprsk=T, xstop=0.5)



#Plotting the curves
plot(km$net.benefit.threshold, km$net.benefit.none, type = "l", lwd=2, xlim=c(0,.50), ylim=c(-.05, .20), xlab = "Threshold Probability", ylab = "Net Benefit")
lines(km$net.benefit$threshold, km$net.benefit$all, type="l", col=3, lwd=2)
# lines(cr$net.benefit$threshold, cr$net.benefit$all, type="l", col=8, lwd=2, lty=2)
lines(km$net.benefit$threshold, km$net.benefit$Gene_signature_nomogram, type="l", col=10 , lwd = 2)
lines(cr$net.benefit$threshold, cr$net.benefit$Model_disintegrating_HPV, type="l", col = 12, lty=2, lwd = 2)
abline(h = 0, lwd =2,col  ="black")
legend("topright", cex=1, legend=c("None", "All", "Gene signature nomogram", "Model disintegrating HPV"), col=c(17, 3, 10, 12), lwd=c(2, 2, 2, 2), lty=c(1, 1, 1,2), bty = "n")


################ NRI & IDI  不一定对  ##########################################
rm(list = ls())
setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step2_modeldisscusion")
HNSCC_nomogram <- read.xlsx2("Training_data.xlsx",row.names = 1, header = T,sheetIndex = 1)

HNSCC_nomogram$futime <- as.numeric(as.character(HNSCC_nomogram$futime))
HNSCC_nomogram$fustatus <- as.numeric(as.character(HNSCC_nomogram$fustatus))
HNSCC_nomogram$Score <- as.numeric(as.character(HNSCC_nomogram$Score))

HNSCC_nomogram[,c(1:10)] <- lapply(HNSCC_nomogram[,c(1:10)],as.numeric)

# ### as vector ####
# HNSCC_nomogram$Gender <- factor(HNSCC_nomogram$Gender,labels = c("Female","male"))
# HNSCC_nomogram$Smoking <- factor(HNSCC_nomogram$Smoking,labels = c("non-Smoking","Smoking"))
# HNSCC_nomogram$Alcohol <- factor(HNSCC_nomogram$Alcohol,labels = c("non-Alcohol","Alcohol"))
# HNSCC_nomogram$T_category <- factor(HNSCC_nomogram$T_category, labels = c("T1&T2","T3&T4"))
# HNSCC_nomogram$N_category <- factor(HNSCC_nomogram$N_category, labels = c("N0&N1","N2","N3"))
# HNSCC_nomogram$HPV_status <- factor(HNSCC_nomogram$HPV_status, labels = c("HPV-","HPV+"))
# HNSCC_nomogram$Age <- factor(HNSCC_nomogram$Age,labels = c("Young","Middle-aged","Young-seniors","Elderly"))

t0=365
indata1=HNSCC_nomogram[,c(2,5,6,7,10)]
indata0=HNSCC_nomogram[,c(2,5,6,10)]
D <- HNSCC_nomogram
n=nrow(D)
covs1<-as.matrix(indata1)
covs0<-as.matrix(indata0)
#--- inference ---
x<-IDI.INF(D[,8:9], covs0, covs1, t0, npert=200)
#--- results ---
IDI.INF.OUT(x)
#--- Graphical presentaion of the estimates ---# 
IDI.INF.GRAPH(x)







