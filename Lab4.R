
#install dependencies
#install.packages("ggplot2") # <- UNCOMMENT IF NOT INSTALLED
library(ggplot2)

ss = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/tutorial\ data-20240923/SS.csv", header=TRUE, row.names=1, sep="\t")
bg = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/tutorial\ data-20240923/BG.csv", header=TRUE, row.names=1, sep="\t")
em = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/tutorial\ data-20240923/EM.csv", header=TRUE, row.names=1, sep="\t")
de_sex = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/tutorial\ data-20240923/DE_sex.csv", header=TRUE, row.names=1, sep="\t")

em[1,]
em[2,]
em[3,]

em[,1]
em[,2]
em[,3]

em[1,1]
em[1,2]
em[1,3]

# The entire 10th column
em[,10]

# The entire 15000th row
em[15000,]

# The value of the 8,123rd row and 2nd column.
em[8123,2]


first_row = em[1,]
first_column = em[,1]
first_cell = em[1,1]

# Now try and select some rows and columns from one of the de_sex by name (not index).
# 1) The row for ENSG00000165246
de_sex['ENSG00000165246',]

#2) The row for ENSG00000072501
de_sex['ENSG00000072501',]

#3) The entire log2fold column
de_sex[,'log2fold']

#4) The entire p column
de_sex[,'p']

#5) The p value for ENSG00000268750
de_sex['ENSG00000268750','p']

# Merge tables
em_annotated = merge(em, bg, by.x=0, by.y=0)
de_sex_annotated = merge(de_sex, bg, by.x=0, by.y=0)

#Rename the rows
row.names(em_annotated) = em_annotated[,'SYMBOL']
row.names(de_sex_annotated) = de_sex_annotated[,'SYMBOL']

#Rename Columns
names(em_annotated)[which(colnames(em_annotated)=='Row.names')] = "GENE_ID"
names(de_sex_annotated)[which(colnames(de_sex_annotated)=='Row.names')] = "GENE_ID"


#Please create a new table called de_annotated_sig which contains rows (genes) in de_annotated
#with a p.adj less than 0.05.

de_sex_annotated_sig = subset(de_sex_annotated, de_sex_annotated['p.adj'] < 0.01)

#Count significant genes
nrow(de_sex_annotated_sig)

#Task 15. Create an Expression Table of Only Significant Genes
#This task is a little trickier. We need to combine some of the things we have learnt today. Our aim is 
#to create a new table called em_symbols_sig. This table should:
#  • contain ONLY genes that are significant
#• contain ONLY columns for samples and NOT things like GENE_ID or CHROMOSOME
#• the bold row names must be the gene symbols (not gene IDs).
#You might need to put these three pieces of code together:
#  row.names(de_annotated_sig) # gives a vector of just the symbols of the significant genes
#names(em) # gives a vector of just the names of the samples
# em_symbols_sig = em_annotated[???,???] # select genes (rows) and samples (columns) needed.

row.names(de_sex_annotated_sig)
names(em)
em_symbols_sig = em_annotated[row.names(de_sex_annotated_sig), names(em)]


ggp = ggplot(em, aes(x=log10(Rep_1))) + geom_histogram(colour="red",  fill="blue", size=2, alpha=0.5) + labs(x="Gene", y="Gene Expression", title="Gene Expression")

ggp = ggplot(em, aes(x=log10(Rep_201))) + geom_density(colour="red",  fill="blue", size=2, alpha=0.5) + labs(x="Gene", y="Gene Expression", title="Gene Expression")

get_gene_data = function(gene, gene_frame) {
  gene_data = gene_frame[gene,]
  gene_data = t(gene_data)
  gene_data = data.frame(gene_data)
  names(gene_data)[1] = 'gene'
  return (gene_data)
}

MPO_gene_data = get_gene_data("MPO", em_symbols_sig)
CD99_gene_data = get_gene_data("CD99", em_symbols_sig)
FIGN_gene_data = get_gene_data("FIGN", em_symbols_sig)
ZFY_gene_data = get_gene_data("ZFY", em_symbols_sig)

ggp = ggplot(ZFY_gene_data, aes(x=log10(gene))) + geom_density(colour="blue", size=2, alpha=0.5) + labs(x="log10(Gene Expression)", y="Number of Samples", title="ZFY gene expression")

get_2gene_data = function(gene1, gene2, gene_frame) {
  gene_data1 = get_gene_data(gene1, gene_frame)
  gene_data2 = get_gene_data(gene2, gene_frame)
  gene_data = merge(gene_data1, gene_data2, by.x=0, by.y=0)
  row.names(gene_data) = gene_data[,'Row.names']
  gene_data = gene_data[-c(1)]
  names(gene_data) = c("gene1", "gene2")
  return (gene_data)
}

MPO_CEACAM8_gene_data = get_2gene_data("MPO","CEACAM8",em_symbols_sig)

ggp = ggplot(MPO_CEACAM8_gene_data, aes(x=log10(gene1), y= log10(gene2))) + geom_point() + labs(x="log10(MPO Gene Expression)", y="log10(CEACAM8 Gene Expression)", title="MPO vs CEACAM8")

ZFY_EVL_gene_data = get_2gene_data("ZFY","EVL",em_symbols_sig)

ggp = ggplot(ZFY_EVL_gene_data, aes(x=log10(gene1), y= log10(gene2))) + geom_point() + labs(x="log10(ZFY Gene Expression)", y="log10(EVL Gene Expression)", title="ZFY vs EVL")

ss = ss[names(em),]

age_gene_data = em["ENSG00000163520",]
age_gene_data = t(age_gene_data)
age_gene_data = merge(age_gene_data, ss, by.x=0, by.y=0)
#age_gene_data$AGE = ss$AGE

ggp = ggplot(age_gene_data, aes(x=AGE, y=log10(ENSG00000163520))) + geom_point()

get_sex_gene_data = function(gene, gene_frame, sample_frame) {
  gene_data = get_gene_data(gene, gene_frame)
  gene_data['sex'] = sample_frame['SEX']
  return (gene_data)
  }

# Lab 4 task 3+4
MPO_sex_gene_data = get_sex_gene_data("MPO", em_symbols_sig, ss)

ggp = ggplot(MPO_sex_gene_data, aes(x=log10(gene), fill=sex)) + geom_density(colour="blue", size=2, alpha=0.5) + labs(x="log10(Gene Expression)", y="Number of Samples", title="MPO gene expression")

# Lab 4 task 5
ZFY_sex_gene_data = get_sex_gene_data("ZFY", em_symbols_sig, ss)

ggp = ggplot(ZFY_sex_gene_data, aes(x=log10(gene), fill=sex)) + geom_density(colour="blue", size=2, alpha=0.5) + labs(x="log10(Gene Expression)", y="Number of Samples", title="ZFY gene expression")

# Lab 4 task 6

ggp = ggplot(MPO_sex_gene_data, aes(x=sex, y=log10(gene))) + geom_boxplot() + labs(x="Sex", y="log10(Gene Expression)", title="MPO gene expression")
ggp = ggplot(ZFY_sex_gene_data, aes(x=sex, y=log10(gene))) + geom_boxplot() + labs(x="Sex", y="log10(Gene Expression)", title="ZFY gene expression")

# task 8
CD99_sex_gene_data = get_sex_gene_data("CD99", em_symbols_sig, ss)

ggp = ggplot(CD99_sex_gene_data, aes(x=sex, y=log10(gene))) + geom_boxplot() + labs(x="Sex", y="log10(Gene Expression)", title="CD99 gene expression")

EVL_sex_gene_data = get_sex_gene_data("EVL", em_symbols_sig, ss)

ggp = ggplot(EVL_sex_gene_data, aes(x=sex, y=log10(gene))) + geom_boxplot() + labs(x="Sex", y="log10(Gene Expression)", title="EVL gene expression")

# task 10
 
halves = cut(ss$AGE, 2, labels=c('young','old'))

#task 11
ss$AGE_halves = halves

#task 12
samples_young = row.names(subset(ss, AGE_halves == 'young'))
samples_old = row.names(subset(ss, AGE_halves == 'old'))

# Task 13. Create em_young and em_old.
em_young = em[samples_young]
em_old = em[samples_old]

# Task 14. Get the Log2Fold and P for the first gene
gene_compare_young_old = function(gene_row, em_young, em_old) {
  gene_data_young = as.numeric(em_young[gene_row,]) # note we select the first row (1) 
  gene_data_old = as.numeric(em_old[gene_row,]) # note we select the first row (1) 
  mean_young = mean(gene_data_young)
  mean_oldest = mean(gene_data_old)
  log2fold = log2(mean_young) - log2(mean_oldest) # mean_old - mean_young
  p = t.test(gene_data_young,gene_data_old) # gene_data_young vs old…
  p = p$p.value
  return (c(log2fold, p))
}

print(gene_compare_young_old(1, em_young, em_old))

#Task 15. Run for the 212th gene
print(gene_compare_young_old(212, em_young, em_old))

#Task 16. Create a table to store the results.

de_age = as.data.frame(matrix(0, ncol = 2, nrow = nrow(em)))
names(de_age) = c('Log2fold','P')
row.names(de_age) = row.names(em)

#Task 17. Add the values for gene 212 to the DE age table
age_stats = gene_compare_young_old(212, em_young, em_old)
age_log2fold = age_stats[1]
age_p = age_stats[2]

print(row.names(em)[212])
de_age[row.names(em)[212],"Log2fold"] = age_log2fold # add to row 212, col 1, of de_age the log2fold value
de_age[row.names(em)[212],"P"] = age_p # add to row 212, col 2, of de_age the p value

#Task 18. Get the values for every gene, using a loop.
for (row in 1:nrow(em))
{
  age_stats = gene_compare_young_old(row, em_young, em_old)
  log2fold = age_stats[1]
  p = age_stats[2]
  
  de_age[row.names(em)[row],"Log2fold"] = log2fold
  de_age[row.names(em)[row],"P"] = p
}

#task 19

get_age_gene_data = function(gene, gene_frame, sample_frame) {
  gene_data = get_gene_data(gene, gene_frame)
  gene_data['age'] = sample_frame['AGE_halves']
  return (gene_data)
}

lowp_age_gene_data = get_age_gene_data("ENSG00000163520", em, ss)

ggp = ggplot(lowp_age_gene_data, aes(x=age, y=log10(gene))) + geom_boxplot() + labs(x="Age group", y="log10(Gene Expression)", title="ENSG00000163520 gene expression")

ggp
