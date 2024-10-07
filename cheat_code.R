

## DAY 2 ####


# load tables
setwd("C:/Users/admin/Desktop/Teaching/2021_2022/MSci_bioinformatics/R_and_Statistics/tutorials/dataset/")

ss = read.table("SS.csv", header=TRUE, row.names = 1, sep="\t")
em = read.table("EM.csv", header=TRUE, row.names = 1, sep="\t")
bg = read.table("BG.csv", header=TRUE, row.names = 1, sep="\t")
de_sex = read.table("DE_sex.csv", header=TRUE, row.names = 1, sep="\t")

# merge
em_annotated = merge(bg,em, by.x=0, by.y=0)
de_sex_annotated = merge(bg,de_sex, by.x=0, by.y=0)

# make row names gene symbols
row.names(em_annotated) = em_annotated$SYMBOL
row.names(de_sex_annotated) = de_sex_annotated$SYMBOL

# rename silly column name
names(em_annotated)[1] = "GENE_ID"
names(de_sex_annotated)[1] = "GENE_ID"

# get sig genes
de_sex_annotated_sig = subset(de_sex_annotated, p.adj < 0.05)


em_symbols_sig = em_annotated[row.names(de_sex_annotated_sig),names(em)]


# make em symbols sig
sig_symbols = row.names(de_sex_annotated_sig)
samples = names(em)
em_symbols_sig = em_annotated[sig_symbols,samples]






## DAY 3 ####

library("ggplot2")

## histograms of samples

# rep1
ggp = ggplot(em, aes(x=log10(Rep_1))) + 
  geom_histogram() + 
  labs(x="expression (log10)", y="number of genes", title="Sample")
ggp

# rep201
ggp = ggplot(em, aes(x=log10(Rep_201))) + 
  geom_histogram() + 
  labs(x="expression (log10)", y="number of genes")
ggp


#--                                         --#
#-- NOW CHOOSE YOUR OWN SAMPLE AND TRY THEM --#
#--                                         --#



## histograms / density of genes

# FGR
gene_data = data.frame(t(em_symbols_sig["FGR",]))
names(gene_data) = "gene"

ggp = ggplot(gene_data, aes(x=gene)) + 
  geom_dotplot(fill="blue") + 
  labs(x="expression (log10)", y="number of samples")
ggp

ggp = ggplot(gene_data, aes(x=gene)) + 
  geom_histogram(fill="blue") + 
  labs(x="expression (log10)", y="number of samples")
ggp

ggp = ggplot(gene_data, aes(x=gene)) + 
  geom_density(fill="blue") + 
  labs(x="expression (log10)", y="density of samples")
ggp

# ZFY
gene_data = data.frame(t(em_symbols_sig["ZFY",]))
names(gene_data) = "gene"

ggp = ggplot(gene_data, aes(x=gene)) + 
  geom_dotplot(fill="blue") + 
  labs(x="expression (log10)", y="number of samples")
ggp

ggp = ggplot(gene_data, aes(x=gene)) + 
  geom_histogram(fill="blue") + 
  labs(x="expression (log10)", y="number of samples")
ggp

ggp = ggplot(gene_data, aes(x=gene)) + 
  geom_density(fill="blue") + 
  labs(x="expression (log10)", y="density of samples")
ggp



#--                                        --#
#-- NOW CHOOSE YOUR OWN GENES AND TRY THEM --#
#--                                        --#


## scatterplots of two genes 

# Correlation
genes = c("MPO", "CEACAM8")
gene_data = data.frame(t(em_symbols_sig[genes,]))
names(gene_data) = c("gene1", "gene2")

ggp = ggplot(gene_data, aes(x=gene1, y=gene2)) + 
  geom_point() + 
  labs(x="expression MPO", y="expression CEACAM8")
ggp

ggp = ggplot(gene_data, aes(x=gene1, y=gene2)) + 
  geom_density_2d() + 
  labs(x="expression TM9SF2", y="expression EVL")
ggp

# Anti-correlation
genes = c("TM9SF2", "EVL")
gene_data = data.frame(t(em_symbols_sig[genes,]))
names(gene_data) = c("gene1", "gene2")

ggp = ggplot(gene_data, aes(x=gene1, y=gene2)) + 
  geom_point() + 
  labs(x="expression TM9SF2", y="expression EVL")
ggp

ggp = ggplot(gene_data, aes(x=gene1, y=gene2)) + 
  geom_density_2d() + 
  labs(x="expression TM9SF2", y="expression EVL")
ggp


#--                                        --#
#-- NOW CHOOSE YOUR OWN GENES AND TRY THEM --#
#--                                        --#


## gene vs Age

# get ss in same sample order as em
ss = ss[names(em),]

# get gene data
gene = "ENSG00000163520"
gene_data = data.frame(t(em[gene,]))
names(gene_data) = "gene"

# add samples column
gene_data$age = ss$AGE

ggp = ggplot(gene_data, aes(x=gene, y=age)) + 
  geom_point() + 
  labs(x="ENSG00000163520", y="Age")
ggp



## DAY 4 ####

## make a gene density split

# MAKE SURE SS is in same order as EM
ss = ss[names(em),]

# FGR
gene_data = data.frame(t(em_symbols_sig["FGR",]))
names(gene_data) = "gene"
gene_data$sex = ss$SEX

ggp = ggplot(gene_data, aes(x=gene, fill=sex)) + 
  geom_density(alpha = 0.5) + 
  labs(x="expression", y="density of samples")
ggp


# ZFY - HIST
gene_data = data.frame(t(em_symbols_sig["ZFY",]))
names(gene_data) = "gene"
gene_data$sex = ss$SEX

ggp = ggplot(gene_data, aes(x=gene, fill=sex)) + 
  geom_histogram(alpha = 0.5) + 
  labs(x="expression", y="density of samples")
ggp


## make a boxplot split

# FGR
gene_data = data.frame(t(em_symbols_sig["FGR",]))
names(gene_data) = "gene"
gene_data$sex = ss$SEX

ggp = ggplot(gene_data, aes(x=sex,y=gene, fill=sex)) + 
  geom_boxplot() + 
  labs(x="sex", y="expression")
ggp

# ZFY
gene_data = data.frame(t(em_symbols_sig["ZFY",]))
names(gene_data) = "gene"
gene_data$sex = ss$SEX

ggp = ggplot(gene_data, aes(x=sex,y=gene, fill=sex)) + 
  geom_boxplot() + 
  labs(x="sex", y="expression")
ggp


#--                                        --#
#-- NOW CHOOSE YOUR OWN GENES AND TRY THEM --#
#--                                        --#


## Get DE age for one gene


# split age into 2 groups
halves = cut(ss$AGE, 2, labels=c("youngest","oldest"))
ss$AGE_halves = halves

# get the younger and older expression tables
samples_youngest = row.names(subset(ss, AGE_halves == "youngest"))
samples_oldest = row.names(subset(ss, AGE_halves == "oldest"))
em_youngest = em[,samples_youngest]
em_oldest = em[,samples_oldest]

# Test a Gene
gene_data_youngest = as.numeric(em_youngest["ENSG00000000003",])
gene_data_oldest = as.numeric(em_oldest["ENSG00000000003",])

mean_youngest = mean(gene_data_youngest)
mean_oldest = mean(gene_data_oldest)
log2fold = log2(mean_oldest) - log2(mean_youngest)

p = t.test(gene_data_oldest,gene_data_youngest)
p = p$p.value


## Make DE age for all genes

# create a new empty table to store all the results in
de_age = as.data.frame(matrix(0, ncol = 2, nrow = nrow(em)))
names(de_age) = c("Log2fold","P")
row.names(de_age) = row.names(em)

# start the loop - nore this is simply the code from "test a gene" with minor modifications.
for (row in 1:nrow(em))
{
  gene_data_youngest = as.numeric(em_youngest[row,])
  gene_data_oldest = as.numeric(em_oldest[row,])
  mean_youngest = mean(gene_data_youngest)
  mean_oldest = mean(gene_data_oldest)
  log2fold = log2(mean_oldest) - log2(mean_youngest)
  p = t.test(gene_data_oldest,gene_data_youngest)
  p = p$p.value
  
  de_age[row,"Log2fold"] = log2fold
  de_age[row,"P"] = p
}

#What is the top gene by P?
#ENSG00000163520

#--                                                          --#
#-- NOW CHOOSE YOUR OWN GENES FROM DE AGE AND MAKE A BOXPLOT --#
#--                                                          --#

#ENSG00000163520
gene_data = data.frame(t(em["ENSG00000163520",]))
names(gene_data) = "gene"
gene_data$AGE_quartile = ss$AGE_quartile

ggp = ggplot(gene_data, aes(x=AGE_quartile,y=gene, fill=AGE_quartile)) + 
  geom_boxplot() + 
  labs(x="Age Group", y="expression")
ggp

