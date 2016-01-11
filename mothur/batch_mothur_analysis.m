#######################
# written by Matthew Scholz
# take outputs from batch_mothur_shared.m and use them for analysis


#first read in variables from last logfile
load.logfile(name=mothur.shared)
set.logfile(name=mothur.analysis)


##charts for alpha diversity accuracy
#alpha = richness/diversity
#richness = # OTUs
#diversity = evenness of distribution of OTUs


#generate plotable stats for Chao1 and invsimpson for each sample, 100 otu intervals

#chao2 is measure of richness
#invsimpson = diversity
collect.single(shared=current,calc=chao-invsimpson,freq=100)

#generate rarefaction data using sobs (observations) 100 otu intervals
rarefaction.single(stability=current,calc=sobs,freq=100)

#summary single ##double check notes for which calcs are useful
summary.single(shared=stability.an.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=T)

##End alpha diversity

##Beta Diversity
#Beta diversity at some point requires table of samples and conditions (example from MiSeq SOP data)

#format: SampleName\tCondition
#for additional groupings,create separate files
#fancy bashy option, create spreadsheet, saved as csv with all metadata 1/column, for each analysis:
#for i in {2..(last column of data)}; do cut ....; done

#pretty pictures
#determine num OTUs to use before running start with 50:
heatmap.bin(shared=current,scale=log2,numotu=50)

#distribution measurements (using thetayc)
#subsample = subsample to the lowest number of OTUs in all samples
#now using Jaccard
dist.shared(shared=current, calc=jclass, subsample=T)
heatmap.sim(phylip=current)

dist.shared(shared=current, calc=thetayc, subsample=T)
heatmap.sim(phylip=current)


#thetayc tree is used for following analysis


#summary.shared() makes an analysis based on declared groups: not really useful for batch analysis
#summary.shared(calc=sharedchao,groups=XXX-XYY-XYX)
#venn(groups=XXX-XYY-XYX)

#probably requires groups to work.
tree.shared(phylyp=current)

#assume group file = "experimental.design"

#let's look at the likelihood that parsimony is similar/diff between groups:
parsimony(tree=current,group=experimental.design,groups=all)

#what about using unifrac weighting to tree distance ?
unifrac.weighted(tree=current,group=experimental.design,groups=all)

#let's look at PCOA for the thetayc output
#we probably won't like it compared to nmds:

pcoa(phylip=current)


#Nmds using 2-4 dimensions:
nmds(phylip=current, groups=experimental.design, mindim=2, maxdim=4)

#what about the ACTUAL beta diversity of groups?  
#amova tests difference between groups vs. variance of groups using "centroid" of each groups' OTU differences
amova(phylip=current, design=experimental.design)

#what about the variance itself? Does it change between groups?
homova(phylip=current, design=experimental.design)


