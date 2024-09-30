
#install dependencies
install.packages("ggplot2")
library(ggplot2)

ss = read.table("M:/R Lab/tutorial\ data-20240930/SS.csv", header=FALSE, sep="\t")

ss = read.table("M:/R Lab/tutorial\ data-20240930/SS.csv", header=TRUE, row.names=1, sep="\t")
bg = read.table("M:/R Lab/tutorial\ data-20240930/BG.csv", header=TRUE, row.names=1, sep="\t")
em = read.table("M:/R Lab/tutorial\ data-20240930/EM.csv", header=TRUE, row.names=1, sep="\t")
de_sex = read.table("M:/R Lab/tutorial\ data-20240930/DE_sex.csv", header=TRUE, row.names=1, sep="\t")

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

ggp
