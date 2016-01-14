##############
#
# batch_mothur_shared.m
#
#  Original Script/comments written by by Dr. Tracy Teal

# This is a batch file for running through the clustering and dist.seqs steps of mothur
# resulting in a shared file
#
# --------------------------------------------------------------------------------- 
#
# We will make use of the 'current' settings, so we don't have to type out the whole 
# file names in each step. This also reduces typo errors and ensures you're using
# the file that was created from the previous step.
#
# This example uses the data and steps in the Schloss MiSeq SOP
# http://www.mothur.org/wiki/MiSeq_SOP
# It does not include the summary.seqs commands and assumes that you have looked
# at your output to ensure that you have made the appropriate decisions for your
# quality filtering parameters. If you haven't gone through this process at least
# once and examined your output, do that first!
# 
##############

load.logfile(logfile=mothur.quality)
set.logfile(name=mothur.shared)

# Set the current files to the ones we generated from the quality sequencing
# At the end we renamed our really long files to something shorter, so we
# can use those names here. Put in your own names if you did something different.

# Calculate the distance matrix
dist.seqs(fasta=current, cutoff=0.2)

# Cluster the sequences
# We're using cluster.split so it can finish more quickly (i.e. in our lifetime)
cluster.split(column=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.15)

# Make the shared file for 0.03 distance
make.shared(list=current, count=current, label=0.03)

#make rarefaction curve
rarefaction.single(shared=current)

#reclassify
classify.otu(list=current,count=current,taxonomy=current,label=0.03)

#generate phylotype information
#phylotype(taxonomy=current)

#make genus level table for phylogeny
#make.shared(list=current, count=current,label=1)

#classify them (names for OTUs)
#classify.otu(list=current,count=current,taxonomy=current,label=1)

#######################
# written by Matthew Scholz
# take outputs from batch_mothur_shared.m and use them for analysis
#first read in variables from last logfile

#phylogenetic diversity assignments:
dist.seqs(fasta=current,output=lt)
clearcut(phylip=current)

count.groups(shared=current)
sub.sample(shared=current)

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


