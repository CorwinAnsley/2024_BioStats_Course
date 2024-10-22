#install dependencies
#install.packages("ggplot2")
library(ggplot2)

#install.packages("ggfortify")
library(ggfortify)

#Import Data
annotations = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/Data_Report_Assessment/Annotations.csv", header=TRUE, row.names=1, sep="\t")
de_gout_vs_hc = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/Data_Report_Assessment/DE_GOUT_vs_HC.csv", header=TRUE, row.names=1, sep="\t")
de_sa_vs_hc = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/Data_Report_Assessment/DE_SA_vs_HC.csv", header=TRUE, row.names=1, sep="\t")
exprn_table = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/Data_Report_Assessment/Expression_Table.csv", header=TRUE, row.names=1, sep="\t")
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
summary(sample_info)

#Plot Neutrophils for each group to compare
ggp = ggplot(sample_info, aes(x=SAMPLE_GROUP, y=NEUTROPHILS)) + geom_boxplot() + labs(x="Sample Group", y="Sample Group", title="Neutrophils in sample groups")
#ggp

#PCA test

PCA = prcomp(t(exprn_table))
pca_coordinates = data.frame(PCA$x)
print(pca_coordinates)

ggp = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = sample_info$SAMPLE_GROUP)) +
  geom_point()
#ggp

# histogram of neutrophils
ggp = ggplot(sample_info, aes(x=NEUTROPHILS)) + 
  geom_histogram() + 
  labs(x="neutrophil_score", y="count", title="Neutrophils")

#ggp

# adds annotations
get_annotated_data = function(data_frame, annotations) {
  df_annotated = merge(data_frame, annotations, by.x=0, by.y=0)
  row.names(df_annotated) = df_annotated[,'Row.names']
  df_annotated = subset(df_annotated, select= -Row.names)
  return (df_annotated)
}

# filters for values below set p threshold
filter_for_sig_genes = function(data_frame, p_max = 0.05, p_column = 'p.adj') {
  data_frame = subset(data_frame, data_frame[p_column] < p_max)
  return(data_frame)
}

de_gout_vs_hc_annotated = get_annotated_data(de_gout_vs_hc, annotations)
de_sa_vs_hc_annotated = get_annotated_data(de_sa_vs_hc, annotations)

de_gout_vs_hc_antd_sig = filter_for_sig_genes(de_gout_vs_hc_annotated)
de_sa_vs_hc_antd_sig = filter_for_sig_genes(de_sa_vs_hc_annotated)

create_plots = function(data_frame, num_plots, sort_by = 'p.adj') {
  plots = c(length(num_plots))
  data_frame = data_frame[order(-data_frame$sort_by),]
  for (i in 1:num_plots)
    plots[i] = 
    
}
