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
biom_file="MGSAT.0.03.biom"
#rarefaction = "example.rarefaction"
nmds = "MGSAT.nmds.axes.nmds"


#import tree
tree<-import_mothur(mothur_tree_file=tree_file)
#importMothur <- import_mothur(mothur_tree_file = tree_file, mothur_list_file = list_file, mothur_group_file = group_file, mothur_shared_file = shared_file, mothur_constaxonomy_file = constaxonomy_file, cutoff=0.03)


#import biom
biom_otu_tax <- import_biom(biom_file)


#let's have fune with the biomfile
GP <- prune_taxa(taxa_sums(biom_otu_tax) > 0, biom_otu_tax)
GP1 <- subset_taxa(GP, Rank4=="Lactobacillales")
plot_bar(GP1, fill="Rank6")
p <- plot_richness(GP, measures=alpha_meas)
p

#rarefaction plot (not implemented)

#stacked bar
#mothur biom uses RankX instead of phylogeny labels in biom file...
plot_bar(biom_otu_tax, fill="Rank4")


#heatmap
plot_heatmap(biom_otu_tax)

plot_tree(tree, label="taxa_names", ladderize="left", plot.margin = 0.5)

#NMDS
c<-read.table(nmds, header=T)

#fix to use experimental.designc[
high_salt<-c[c$Treatment=="high_salt", c(2,3)]
normal_salt<-c[c$Treatment=="normal_salt", c(2,3)]

plot(high_salt, xlim=c(-.9,.9), ylim=c(-.9, .9), pch=16, cex=2, col="black", xlab="NMDS Axis 1", ylab="NMDS Axis 2")
points(normal_salt, pch=16, col="grey", cex=2)

legend(-.99, .99,inset=.5, bty="n", cex=.9, c("high_salt", "normal_salt"), col=c("black", "grey"), pch=c(16,16,16))

