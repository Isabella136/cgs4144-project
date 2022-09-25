library(readxl)
library(ggplot2)
library(dplyr)
data <- read_excel("GSE65683_Required_SRE.xlsx")
logdata <- data.frame(Element_name = data$`Element name`)
rangedata <- data.frame(Element_name = data$`Element name`,
                        range = 0)
col = colnames(data)
for (i in 6:77) {
  logdata[col[i]] <- log(data[,c(i)],2)
}
for (i in 1:528) {
  rangedata[i,2] = max(logdata[,-c(1)][i,]) - min(logdata[,-c(1)][i,])
}
rangedata %>% filter(is.finite(range)) %>% ggplot(aes(x=range)) + geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)