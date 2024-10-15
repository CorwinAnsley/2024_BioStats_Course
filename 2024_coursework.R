#install dependencies
install.packages("ggplot2")
library(ggplot2)

install.packages("ggfortify")
library(ggfortify)

#Import Data
annotations = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/Data_Report_Assessment/Annotations.csv", header=TRUE, row.names=1, sep="\t")
de_gout_vs_hc = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/Data_Report_Assessment/DE_GOUT_vs_HC.csv", header=TRUE, row.names=1, sep="\t")
de_sa_vs_hc = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/Data_Report_Assessment/DE_SA_vs_HC.csv", header=TRUE, row.names=1, sep="\t")
expression_table = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/Data_Report_Assessment/Expression_Table.csv", header=TRUE, row.names=1, sep="\t")
sample_info = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/Data_Report_Assessment/Sample_Information.csv", header=TRUE, row.names=1, sep="\t")

#Are groups well matched 
sample_info$SAMPLE_GROUP=as.factor(sample_info$SAMPLE_GROUP)
sample_info$SEX=as.factor(sample_info$SEX)

sample_info_hc = subset(sample_info, sample_info['SAMPLE_GROUP'] == 'HC')
sample_info_gout = subset(sample_info, sample_info['SAMPLE_GROUP'] == 'GOUT')
sample_info_sa = subset(sample_info, sample_info['SAMPLE_GROUP'] == 'SEPSIS')

summary(sample_info_hc)
summary(sample_info_gout)
summary(sample_info_sa)

#Plot Neutrophils for each group to compare
ggp = ggplot(sample_info, aes(x=SAMPLE_GROUP, y=NEUTROPHILS)) + geom_boxplot() + labs(x="Sample Group", y="Sample Group", title="Neutrophils in sample groups")
ggp
