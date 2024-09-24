ss = read.table("C:/Users/2266643A/repos/2024_BioStats_Course/tutorial\ data-20240923/SS.csv", header=FALSE, sep="\t")

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


