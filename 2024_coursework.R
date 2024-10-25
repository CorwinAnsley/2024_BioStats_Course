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

# Helper function to get data for specific gene
get_gene_data = function(gene, gene_frame) {
  gene_data = gene_frame[gene,]
  gene_data = t(gene_data)
  gene_data = data.frame(gene_data)
  names(gene_data)[1] = 'gene'
  return (gene_data)
}

# Gets gene data along with samplegroup
get_sample_group_gene_data = function(gene, gene_frame, sample_info) {
  gene_data = get_gene_data(gene, gene_frame)
  gene_data['sample_group'] = sample_info['SAMPLE_GROUP']
  return (gene_data)
}


de_gout_vs_hc_annotated = get_annotated_data(de_gout_vs_hc, annotations)
de_sa_vs_hc_annotated = get_annotated_data(de_sa_vs_hc, annotations)

de_gout_vs_hc_antd_sig = filter_for_sig_genes(de_gout_vs_hc_annotated)
de_sa_vs_hc_antd_sig = filter_for_sig_genes(de_sa_vs_hc_annotated)

#de_df_sorted = de_gout_vs_hc_antd_sig[order(de_gout_vs_hc_antd_sig$'p.adj'),,] 

create_plots_gene_ex = function(df, exprn_table, num_plots, sample_info) { #sort_by = "p.adj"
  #plots = c(length(num_plots))
  #de_df_sorted = de_frame[order(de_frame$'p.adj'),,]
  for (i in 1:num_plots) {
    gene_symbol = df[i,'symbol']
    print(gene_symbol)
    gene_data = get_sample_group_gene_data(rownames(df)[i], exprn_table, sample_info)
    #plot_title = gene_symbol
    ggp = ggplot(gene_data, aes(x=sample_group, y=log10(gene), colour=sample_group), title=gene_symbol) +
    geom_point(fill = "blue") + 
    labs(x="samplegroup", y="expression")
    ggp
    filename = paste('C:\\Users\\2266643A\\repos\\2024_BioStats_Course\\Coursework_plots\\plot_', gene_symbol, sep='_')
    #ggsave('C:\\Users\\2266643A\\repos\\2024_BioStats_Course\\plot.png')
    ggsave(paste(filename, '.png'))
  }
  #return(plots)
}

#sort genes by biggest difference in log2fold
get_sorted_genes = function(df1, df2) {
  #df1 <- df1 %>% 
  #  add_column(delta_log2fold = 0)
  df1[,'delta_log2fold'] = 0
  for (i in 1:nrow(df1)) { 
    row_name = rownames(df1)[i]
    df1[i,'delta_log2fold'] = abs(df1[i,'log2Fold'] - df2[row_name,'log2Fold'])
    }
  df1_sorted = df1[order(df1$'delta_log2fold'),,]
  return(df1_sorted)
}

de_sa_vs_hc_antd_sig_sorted = get_sorted_genes(de_sa_vs_hc_antd_sig, de_gout_vs_hc) 
#gout_gene_plots = 
create_plots_gene_ex(de_sa_vs_hc_antd_sig_sorted, exprn_table, 5, sample_info)
#ggp = gout_gene_plots[2]
#ggp
