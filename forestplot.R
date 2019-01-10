setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step5_cox_forestplot")
#install.packages("forestplot")
library("forestplot")

#��ȡ�����ļ�����������һ��������ࡢ�����������ձȼ��������䣨���޼����ޣ��ȣ���Ҫע����������ļ��Ĳ���(������������)��չʾ��������У����������á�
data <- read.csv("forestplotdata.csv", stringsAsFactors=FALSE)
np <- ifelse(!is.na(data$Count), paste(data$Count," (",data$Percent,")",sep=""), NA)#�ϲ�count��percent����һ���Ǳ�Ҫ
head(np)
#��Ҫ��ͼ��չʾ���ı�
tabletext <- cbind(c("\nSubgroup",NA,NA,data$Variable,NA),
                   c("No. of\nPatients (%)",NA,NA,np,NA),
                   c("Hazard Ratio\n(95% CI)",NA,NA,ifelse(!is.na(data$Count), paste(format(data$Point.Estimate,nsmall=2)," (",format(data$Low,nsmall = 2)," to ",format(data$High,nsmall = 2),")",sep=""), NA),NA))
head(tabletext)
forestplot(labeltext=tabletext, #ͼ�е��ı�
           mean=c(NA,NA,1,data$Point.Estimate,NA),#HR
           lower=c(NA,NA,1,data$Low,NA), #95%������������
           upper=c(NA,NA,1,data$High,NA),#95%������������
           #title="Hazard Ratio",
           graph.pos=3,#ͼ�ڱ��е���λ��
           graphwidth = unit(.4,"npc"),#ͼ�ڱ��еĿ��ȱ���
           fn.ci_norm="fpDrawDiamondCI",#box����ѡ����ʯ
           col=fpColors(box="steelblue", lines="black", zero = "black"),#box��ɫ
           boxsize=c(NA,NA,NA,data$Point.Estimate,NA)/10,#box��С��������������
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#�����������߿����ߡ���
           zero=1,#zero�ߺ�����
           lwd.zero=2,#zero�߿�
           grid = structure(c(data[1,]$Point.Estimate), gp = gpar(col = "black", lty=2,lwd=2)),#����(�ɶ���)��������ꡢ��ɫ���߿�
           xticks = c(1,5,10,15),#������̶ȸ�����Ҫ����������
           lwd.xaxis=2,#X���߿�
           xlab="     <-Favors Low Score      Favors High Score->",#X�����
           hrzl_lines=list("3" = gpar(lwd=2, col="black"),#�����ж����Ӻ��ߣ����������ֱ����λ��
                           #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#����Ӱ���������ʹ��
                           "30" = gpar(lwd=2, col="black")),#���һ�еײ��Ӻ���,""������Ϊnrow(data)+5
           txt_gp=fpTxtGp(label=gpar(cex=1.25),#���������С����
                          ticks=gpar(cex=1.25),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.25)),
           #is.summary = c(T,rep(F,27)),#����������������
           lineheight = unit(.75,"cm"),#�̶��и�
           #align=c("l","c","c"),#ÿ�����ֵĶ��뷽ʽ��ż�����õ�
           #cex=10,
           colgap = unit(0,"cm"),#�м�϶
           mar=unit(rep(1.25, times = 4), "cm"),#ͼ��ҳ�߾�
           new_page = F#�Ƿ���ҳ
)
