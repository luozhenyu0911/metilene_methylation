print("################  Processing forth step: Plot the mahattan plot for DMCs sites ##################################")
library(CMplot)
args=commandArgs(T)


var <- args[1]
filef <- paste0("outfile/InGroup2.",var,".DMCs_pval.0.05.out")

data<- read.table(filef,sep = '\t')
# colnames(data)<-c('Chr', 'start', 'End', 'Diff', 
#                   'qval', 'count','Amean', 'Bmean')
head(data)

data$name <- paste0(data$V1,data$V2)

data2<-data[,c(8,1,2,3)]

CMplot(data2,plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),
       file="pdf",file.name='plot2',main = "DMCs distribution diagram",dpi=300,file.output=TRUE, verbose=TRUE)


CMplot(data2,plot.type="m",LOG10=TRUE,threshold=NULL,chr.den.col=NULL,
       file="pdf",file.name='plot3',main = "DMCs distribution diagram",dpi=300,file.output=TRUE, verbose=TRUE)

print("################  Processing forth step: Plot the mahattan plot for DMCs sites ##################################")
