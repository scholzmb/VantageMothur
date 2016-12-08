#export biomfile in mothur:
# make.biom(shared=example.shared, label=0.03, constaxonomy=example.cons.taxonomy, matrixtype=dense, reftaxonomy=silva.bacteria.silva.tax)
# copy appropriate tre file

source("https://bioconductor.org/biocLite.R")

library("ggplot2")
library("phyloseq")
theme_set(theme_bw())

#files (example.x), to be renamed later
tree_file="MGSAT.0.03.thetayc.tre"
list_file="MGSAT.list"
group_file="MGSAT.groups"
shared_file="MGSAT.shared"
constaxonomy_file="MGSAT.cons.taxonomy"
tax_table="Rank4.txt"
biom_file="MGSAT.0.03.biom"
#rarefaction = "example.rarefaction"
nmds = "MGSAT.nmds.axes.nmds"


#import tree
tree<-import_mothur(mothur_tree_file=tree_file)

#importMothur <- import_mothur(mothur_tree_file = tree_file, mothur_list_file = list_file, mothur_group_file = group_file, mothur_shared_file = shared_file, mothur_constaxonomy_file = constaxonomy_file, cutoff=0.03)

#NMDS
c<-read.table(nmds, header=T)

#import biom
biom_otu_tax <- import_biom(biom_file)

pdf("16SCharts.pdf")

#let's have fune with the biomfile
GP <- prune_taxa(taxa_sums(biom_otu_tax) > 0, biom_otu_tax)

GP2= sample_data(c)
GP2 = merge_phyloseq(GP2, GP)
GP3 <- GP2
#NGP= transform_sample_counts(GP3, threshrankfun(100))
NGP <- GP3
plot_richness(GP3, measures=c("Observed","InvSimpson", "Shannon"))

#plot_bar(NGP, fill="Rank4")
#p <- plot_richness(GP, measures=alpha_meas)
#p

#rarefaction plot (not implemented)

#stacked bar
#mothur biom uses RankX instead of phylogeny labels in biom file...
plot_bar(biom_otu_tax, fill="Rank2")

#family=c("Enterobacteriales")
for(level in c("4")){
  tax_table = paste("Rank", level, ".txt", sep="")
  taxa <- read.table(tax_table, header=F, skip=1)
  family <- as.vector(taxa[,1])
  NGP <- GP3
  for (i in family){
    MBS <- paste(i,sep="")
    GP4 <- subset_taxa(NGP, Rank4==MBS) 
    print(plot_bar(GP4, fill="Treatment", title=MBS) +aes(group= Treatment ))
    #print(plot_heatmap(GP4, fill="OTU", title=MBS, taxa.label="Rank6"))
  }
}

for(level in c("5")){
  tax_table = paste("Rank", level, ".txt", sep="")
  taxa <- read.table(tax_table, header=F, skip=1)
  family <- as.vector(taxa[,1])
  NGP <- GP3
  for (i in family){
    MBS <- paste(i,sep="")
    GP4 <- subset_taxa(NGP, Rank5==MBS) 
    print(plot_bar(GP4, fill="Treatment", title=MBS) +aes(group= Treatment ))
    #print(plot_heatmap(GP4, fill="OTU", title=MBS, taxa.label="Rank6"))
  }
}


#heatmaps
NGP <- GP3
plot_heatmap(NGP, sample.order="Treatment", sample.label="Treatment", taxa.label="Rank4", title="All Taxa")
plot_tree(tree, label="taxa_names", ladderize="left", plot.margin = 0.5)

NGP = prune_taxa(names(sort(taxa_sums(GP3),TRUE)[1:50]), GP3)
plot_heatmap(NGP, taxa.label="Rank4", sample.label="Treatment", sample.order = "Treatment", title = "50 Most abundant Taxa")

ordinations= c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")
for (i in ordinations){
  ord <- paste(i,sep="")
  label <- paste(c(i, "Ordination Plot", sep=" "))
  PCA = ordinate(NGP, method=ord, distance = "bray")
  j<- plot_ordination(NGP, PCA, color = "Treatment", title=label)+geom_point(size=8)
  print(j)
  
}

dev.off()