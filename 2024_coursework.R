#####Install Dependencies#####

#install.packages("ggplot2")
library(ggplot2)

#install.packages("ggfortify")
#library(ggfortify)

install.packages("xtable")
library("xtable")
options(xtable.floating=FALSE)
options(xtable.timestamp="")

install.packages("ggrepel")
library(ggrepel)

install.packages("dplyr")
library(dplyr)


#####Import Data#####

#Define directories wherE data is stored and saved
#DATADIR = "C:/Users/2266643A/repos/2024_BioStats_Course/Data_Report_Assessment"
DATADIR = "C:/Users/Curry/repos/Bioinfo_Msc/2024_BioStats_Course/Data_Report_Assessment"

#Define directory/filename to save plots
#PLOTSDIR = "C:\\Users\\2266643A\\repos\\2024_BioStats_Course\\Coursework_plots\\plot_"
PLOTSDIR = "C:\\Users\\CURRY\\repos\\Bioinfo_Msc\\2024_BioStats_Course\\Coursework_plots\\"

annotations = read.table(paste(DATADIR,"/Annotations.csv",sep=""), header=TRUE, row.names=1, sep="\t")
de_gout_vs_hc = read.table(paste(DATADIR,"/DE_GOUT_vs_HC.csv",sep=""), header=TRUE, row.names=1, sep="\t")
de_sa_vs_hc = read.table(paste(DATADIR,"/DE_SA_vs_HC.csv",sep=""), header=TRUE, row.names=1, sep="\t")
exprn_table = read.table(paste(DATADIR,"/Expression_Table.csv",sep=""), header=TRUE, row.names=1, sep="\t")
sample_info = read.table(paste(DATADIR,"/Sample_Information.csv",sep=""), header=TRUE, row.names=1, sep="\t")

sample_info$SAMPLE_GROUP=as.character(sample_info$SAMPLE_GROUP)
sample_info$SAMPLE_GROUP[sample_info$SAMPLE_GROUP == 'SEPSIS'] = "SA"
sample_info$SAMPLE_GROUP=as.factor(sample_info$SAMPLE_GROUP)

sample_info$SEX=as.factor(sample_info$SEX)

sample_info_hc = subset(sample_info, sample_info['SAMPLE_GROUP'] == 'HC')
sample_info_gout = subset(sample_info, sample_info['SAMPLE_GROUP'] == 'GOUT')
sample_info_sa = subset(sample_info, sample_info['SAMPLE_GROUP'] == 'SEPSIS')

sample_info_hc_vs_gout = subset(sample_info, sample_info['SAMPLE_GROUP'] != 'SEPSIS')
sample_info_hc_vs_sa = subset(sample_info, sample_info['SAMPLE_GROUP'] != 'GOUT')
sample_info_sa_vs_gout = subset(sample_info, sample_info['SAMPLE_GROUP'] != 'HC')

summary(sample_info_hc)
summary(sample_info_gout)
summary(sample_info_sa)
summary(sample_info)


#Part 1
####Analysing  Neutrophils for each group to compare

#PLotting Neutrophils VS samplegroup
ggp = ggplot(sample_info, aes(x=SAMPLE_GROUP, y=NEUTROPHILS, colour=SAMPLE_GROUP)) + 
  geom_boxplot() + 
  labs(x="Sample Group", y="Neutrophil Count") +
  theme(legend.position="none")
ggp

ggp = ggplot(sample_info, aes(x=SAMPLE_GROUP, y=NEUTROPHILS, colour=SAMPLE_GROUP), title="Neutrophils in sample groups") +
  geom_point(fill = "blue") + 
  labs(x="samplegroup", y="neutrophils") +
  theme(legend.position="none")
ggp

#Comparing Neutrophils by sex, just to check
ggp = ggplot(sample_info, aes(x=SEX, y=NEUTROPHILS, colour=SEX)) + 
  geom_boxplot() + 
  labs(x="Sex", y="Neutrophil count", title="Neutrophils in sample groups") +
  theme(legend.position="none")
ggp

#Performing anova
model_neutrophils = lm(sample_info$NEUTROPHILS~sample_info$SAMPLE_GROUP)
xtable(anova(model_neutrophils))
anova(model_neutrophils)
summary(model_neutrophils)

#t.test(sample_info_hc_vs_gout$NEUTROPHILS~sample_info_hc_vs_gout$SAMPLE_GROUP)
#t.test(sample_info_hc_vs_sa$NEUTROPHILS~sample_info_hc_vs_sa$SAMPLE_GROUP)
#t.test(sample_info_sa_vs_gout$NEUTROPHILS~sample_info_sa_vs_gout$SAMPLE_GROUP)

#PCA test
PCA = prcomp(t(exprn_table))
pca_coordinates = data.frame(PCA$x)
print(pca_coordinates)

ggp = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = sample_info$SAMPLE_GROUP)) +
  geom_point()
#ggp


# adds annotations
get_annotated_data = function(data_frame, annotations) {
  df_annotated = merge(data_frame, annotations, by.x=0, by.y=0)
  row.names(df_annotated) = df_annotated[,'Row.names']
  df_annotated = subset(df_annotated, select= -Row.names)
  return (df_annotated)
}

# filters for values below set p threshold
filter_for_sig_genes = function(data_frame, p_max = 0.05, log2Fold_threshold = 0, p_column = 'p.adj', log2Fold_column = 'log2Fold') {
  data_frame = subset(data_frame, data_frame[p_column] < p_max)
  data_frame = subset(data_frame, abs(data_frame[log2Fold_column]) > log2Fold_threshold)
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

# Gets gene data along with info about sample
get_sample_info_gene_data = function(gene, gene_frame, sample_info) {
  gene_data = get_gene_data(gene, gene_frame)
  gene_data['sample_group'] = sample_info['SAMPLE_GROUP']
  gene_data['sex'] = sample_info['SEX']
  gene_data['neutrophils'] = sample_info['NEUTROPHILS']
  return (gene_data)
}

#sort genes by biggest difference in log2fold
get_sorted_log2fold_diff_genes = function(df1, df2) {
  #df1 <- df1 %>% 
  #  add_column(delta_log2fold = 0)
  df1[,'delta_log2fold'] = 0
  for (i in 1:nrow(df1)) { 
    row_name = rownames(df1)[i]
    df1[i,'delta_log2fold'] = abs(df1[i,'log2Fold'] - df2[row_name,'log2Fold'])
  }
  df1_sorted = df1[order(df1$'delta_log2fold',decreasing=TRUE ),,]
  return(df1_sorted)
}


#de_df_sorted = de_gout_vs_hc_antd_sig[order(de_gout_vs_hc_antd_sig$'p.adj'),,] 

#create plots of gene expression by sample group
create_plots_gene_ex = function(df, exprn_table, num_plots, sample_info, dir = '') { 
  for (i in 1:num_plots) {
    gene_symbol = df[i,'symbol']
    print(gene_symbol)
    gene_data = get_sample_info_gene_data(rownames(df)[i], exprn_table, sample_info)
    #plot_title = gene_symbol
    ggp = ggplot(gene_data, aes(x=sample_group, y=log10(gene), colour=sample_group), title=gene_symbol) +
    geom_point(fill = "blue") + 
    labs(x="samplegroup", y="expression") +
      scale_color_manual(values=c("darkcyan", "chartreuse3", "darkred")) +
    theme(legend.position="none")
    ggp
    filedir = paste(PLOTSDIR, dir, sep='\\')
    print(filedir)
    filename = paste('\\plot_', gene_symbol, sep='')
    print(filename)
    filepath = paste(filedir, filename, sep='\\')
    #ggsave('C:\\Users\\2266643A\\repos\\2024_BioStats_Course\\plot.png')
    ggsave(paste(filepath, '.pdf'))
  }
}

#create plots of gene expression vs neutrophil count
create_plots_gene_ex_neutrophils = function(df, exprn_table, num_plots, sample_info, dir = '') { 
  for (i in 1:num_plots) {
    gene_symbol = df[i,'symbol']
    print(gene_symbol)
    gene_data = get_sample_info_gene_data(rownames(df)[i], exprn_table, sample_info)
    #plot_title = gene_symbol
    ggp = ggplot(gene_data, aes(x=neutrophils, y=log10(gene), colour=sample_group), title=gene_symbol) +
      geom_point(fill = "blue") + 
      labs(x="Neutrophil Count", y="Gene Expression") +
      theme(legend.position="none")
    ggp
    filedir = paste(PLOTSDIR, dir, sep='\\')
    print(filedir)
    filename = paste('\\plot_neutrophils_', gene_symbol, sep='')
    print(filename)
    filepath = paste(filedir, filename, sep='\\')
    #ggsave('C:\\Users\\2266643A\\repos\\2024_BioStats_Course\\plot.png')
    ggsave(paste(filepath, '.pdf'))
  }
}

 

##Investigating significant genes for gout and SA

#Creates a volcano plot highlighting nad labeling the most significant genes
volcano_plot_regulated_genes = function(df, p_max = 0.05, log2Fold_threshold = 1.2, p_column = 'p.adj', log2Fold_column = 'log2Fold', symbol_labels = TRUE) {
  # adding label to de tables for up and down regulated genes
  df$diffexpr = "NO" 
  df$diffexpr[df[log2Fold_column] > log2Fold_threshold & df[p_column] < p_max] = "UP"
  df$diffexpr[df[log2Fold_column] < -log2Fold_threshold  & df[p_column] < p_max] = "DOWN"
  
  # adding name labels for the significant genes
  df$delabel = NA
  if (symbol_labels) {
    df$delabel[df$diffexpr != "NO"] = df$symbol[df$diffexpr != "NO"]
  }
  
  ggp = ggplot(data=df, aes(x=log2Fold, y=-log10(p.adj), col=diffexpr, label=delabel)) + 
    geom_point() +
    theme_minimal() +
    geom_text_repel(max.overlaps=100) +
    scale_color_manual(values=c("darkcyan", "black", "darkred")) +
    geom_vline(xintercept=c(-log2Fold_threshold , log2Fold_threshold ), col="red") +
    geom_hline(yintercept=-log10(p_max), col="red") +
    theme(legend.position="none")
  ggp
}

de_gout_vs_hc_annotated = get_annotated_data(de_gout_vs_hc, annotations)
de_sa_vs_hc_annotated = get_annotated_data(de_sa_vs_hc, annotations)

volcano_plot_regulated_genes(de_gout_vs_hc_annotated)
#volcano_plot_regulated_genes(de_sa_vs_hc_annotated, log2Fold_threshold = 8.4, symbol_labels=TRUE)

de_gout_vs_hc_antd_sig = filter_for_sig_genes(de_gout_vs_hc_annotated,log2Fold_threshold=1.2)
de_sa_vs_hc_antd_sig = filter_for_sig_genes(de_sa_vs_hc_annotated,log2Fold_threshold=8.4)

#sort the genes by p values
de_gout_vs_hc_antd_sig = de_gout_vs_hc_antd_sig[order(de_gout_vs_hc_antd_sig$p.adj,decreasing=FALSE ),,]
de_sa_vs_hc_antd_sig = de_sa_vs_hc_antd_sig[order(de_sa_vs_hc_antd_sig$p.adj,decreasing=FALSE ),,]


create_plots_gene_ex(de_gout_vs_hc_antd_sig, exprn_table, 11, sample_info, dir='GOUT_genes')
create_plots_gene_ex(de_sa_vs_hc_antd_sig, exprn_table, 11, sample_info, dir='SA_genes')

create_plots_gene_ex_neutrophils(de_gout_vs_hc_antd_sig, exprn_table, 11, sample_info, dir='GOUT_genes')
create_plots_gene_ex_neutrophils(de_sa_vs_hc_antd_sig, exprn_table, 11, sample_info, dir='SA_genes')

#de_sig_antd_sa_gout = merge(de_sa_vs_hc_antd_sig, de_gout_vs_hc_antd_sig, by.x=0, by.y=0)

#Genes to focus on:

PI3_gene_data = get_sample_info_gene_data('ENSG00000124102', exprn_table, sample_info)
SPRR2C_gene_data = get_sample_info_gene_data('ENSG00000229035', exprn_table, sample_info)
GATD3A_gene_data = get_sample_info_gene_data('ENSG00000160221', exprn_table, sample_info)

get_anova_table_neutrophil_genex = function(gene_data) {
  model1 = lm(gene_data$neutrophils~gene_data$gene)
  xtable(anova(model1))
}
#Are they affected by other factors?
get_anova_table_neutrophil_genex(PI3_gene_data)
get_anova_table_neutrophil_genex(SPRR2C_gene_data)
get_anova_table_neutrophil_genex(GATD3A_gene_data)

t.test(GATD3A_gene_data$gene~GATD3A_gene_data$sex)
t.test(PI3_gene_data$gene~PI3_gene_data$sex)
t.test(SPRR2C_gene_data$gene~SPRR2C_gene_data$sex)


##DO ALL THIS FOR DIFFERENTIAL GENES BETWEEN GOUT AND SA

# Create de table for two groups
create_de_table = function(exprn_table, sample_info, group1, group2) {
  de = data.frame() 
  samples_group1 = row.names(subset(sample_info, SAMPLE_GROUP == group1))
  samples_group2 = row.names(subset(sample_info, SAMPLE_GROUP == group2))
  
  # Task 13. Create em_young and em_old.
  em_group1 = exprn_table[samples_group1]
  em_group2 = exprn_table[samples_group2]
  for (gene in rownames(exprn_table)) {
    gene_data_group1 = as.numeric(em_group1[gene,])
    gene_data_group2 = as.numeric(em_group2[gene,])
    mean_group1 = mean(gene_data_group1)
    mean_group2 = mean(gene_data_group2)
    log2fold = log2(mean_group1) - log2(mean_group2)
    p = t.test(gene_data_group1,gene_data_group2)
    p = p$p.value
    de[gene,"log2Fold"] = log2fold
    de[gene,"p"] = p
  }
  p_adj = p.adjust (de[,"p"], method="fdr")
  de[,"p.adj"] = p_adj
  return (de)
}

exprn_table_small = exprn_table[1:10, ]
de_sa_vs_gout = create_de_table(exprn_table, sample_info, "SA", "GOUT")

de_sa_vs_gout_annotated = get_annotated_data(de_sa_vs_gout, annotations)

volcano_plot_regulated_genes(de_sa_vs_gout_annotated,log2Fold_threshold=8.4)

de_sa_vs_gout_antd_sig = filter_for_sig_genes(de_sa_vs_gout_annotated,log2Fold_threshold=8.4)

de_sa_vs_hc_antd_sig_sorted = get_sorted_genes(de_sa_vs_hc_antd_sig, de_gout_vs_hc) 
#gout_gene_plots = 
create_plots_gene_ex(de_sa_vs_hc_antd_sig_sorted, exprn_table, 20, sample_info)