setwd("D:\\PHD_database\\paper_writing\\nomogram\\data\\HNSCC\\GEO\\result_step5_cox_forestplot")
#install.packages("forestplot")
library("forestplot")

#读取输入文件，输入数据一般包括分类、样本数、风险比及置信区间（上限及下限）等，需要注意的是输入文件的布局(包括文字缩进)将展示在最后结果中，所见即所得。
data <- read.csv("forestplotdata.csv", stringsAsFactors=FALSE)
np <- ifelse(!is.na(data$Count), paste(data$Count," (",data$Percent,")",sep=""), NA)#合并count和percent，这一步非必要
head(np)
#将要在图中展示的文本
tabletext <- cbind(c("\nSubgroup",NA,NA,data$Variable,NA),
                   c("No. of\nPatients (%)",NA,NA,np,NA),
                   c("Hazard Ratio\n(95% CI)",NA,NA,ifelse(!is.na(data$Count), paste(format(data$Point.Estimate,nsmall=2)," (",format(data$Low,nsmall = 2)," to ",format(data$High,nsmall = 2),")",sep=""), NA),NA))
head(tabletext)
forestplot(labeltext=tabletext, #图中的文本
           mean=c(NA,NA,1,data$Point.Estimate,NA),#HR
           lower=c(NA,NA,1,data$Low,NA), #95%置信区间下限
           upper=c(NA,NA,1,data$High,NA),#95%置信区间上限
           #title="Hazard Ratio",
           graph.pos=3,#图在表中的列位置
           graphwidth = unit(.4,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           col=fpColors(box="steelblue", lines="black", zero = "black"),#box颜色
           boxsize=c(NA,NA,NA,data$Point.Estimate,NA)/10,#box大小根据样本量设置
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           grid = structure(c(data[1,]$Point.Estimate), gp = gpar(col = "black", lty=2,lwd=2)),#虚线(可多条)及其横坐标、颜色、线宽
           xticks = c(1,5,10,15),#横坐标刻度根据需要可随意设置
           lwd.xaxis=2,#X轴线宽
           xlab="     <-Favors Low Score      Favors High Score->",#X轴标题
           hrzl_lines=list("3" = gpar(lwd=2, col="black"),#第三行顶部加黑线，引号内数字标记行位置
                           #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#加阴影，弱项不建议使用
                           "30" = gpar(lwd=2, col="black")),#最后一行底部加黑线,""中数字为nrow(data)+5
           txt_gp=fpTxtGp(label=gpar(cex=1.25),#各种字体大小设置
                          ticks=gpar(cex=1.25),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.25)),
           #is.summary = c(T,rep(F,27)),#首行字体类型设置
           lineheight = unit(.75,"cm"),#固定行高
           #align=c("l","c","c"),#每列文字的对齐方式，偶尔会用到
           #cex=10,
           colgap = unit(0,"cm"),#列间隙
           mar=unit(rep(1.25, times = 4), "cm"),#图形页边距
           new_page = F#是否新页
)

