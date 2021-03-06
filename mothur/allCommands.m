#Mothur first steps assembled 2/8/2016 by matthew scholz
set.logfile(name=testingSetCurrent.log)

#assemble/qc/qa/count
#make.contigs(file=stability.files, processors=8)
set.current(fasta=stability.trim.contigs.good.fasta, group=stability.contigs.good.groups)
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275)
unique.seqs(fasta=current)
count.seqs(name=current, group=current)
summary.seqs(count=current)
quit()


get.current()
set.current(fasta=stability.trim.contigs.good.unique.fasta, group=stability.contigs.good.groups, name=stability.trim.contigs.good.names, qfile=stability.trim.contigs.qual, count=stability.trim.contigs.good.count_table, processors=8,summary=stability.trim.contigs.good.unique.summary)

#remove all lineages present in blank 
summary.seqs(fasta=current, count=current)
#get.groups(fasta=current, count=current, groups=Blank_1-Blank_2-Blank_3-Blank_4)
#list.seqs(fasta=current)
#remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, accnos=current, name=current)
#remove.seqs(count=stability.trim.contigs.good.unique.good.filter.count_table, accnos=current)
#list.seqs(count=current)
#summary.seqs(fasta=current, count=current)
#get.current()

#Defaults for refDB
#To use silva/greengenes
#(reference=silva.nr_v123.align, taxonomy=silva.nr_v123.tax) 
#to use RDP (recommended? for fecal microbiome?)
#(reference=trainset14_032015.rdp.fasta, taxonomy=trainset14_032015.rdp.tax)

align.seqs(fasta=current, template=silva.nr_v123.align)
summary.seqs(fasta=current, count=current)

#for Silva/greengenes
screen.seqs(fasta=current, count=current, summary=current, start=13862, end=23444, maxhomop=8, processors=8)
#for RDP
#screen.seqs(fasta=current, count=current, summary=current, start=13862, end=23444, maxhomop=8, processors=8)

summary.seqs(fasta=current, count=current)
filter.seqs(fasta=current, vertical=T, trump=.)
unique.seqs(fasta=current, count=current)
get.current()

#classify/remove useless
align.seqs(fasta=current, reference=silva.nr_v123.align)
classify.seqs(fasta=current, count=current, reference=silva.nr_v123.align, taxonomy=silva.nr_v123.tax, cutoff=80)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
get.current()

#cluster
pre.cluster(fasta=current, count=current, diffs=2)
chimera.uchime(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current)
summary.seqs(fasta=current, count=current)


#OTU Analysis
classify.seqs(fasta=current, count=current, reference=silva.nr_v123.align, taxonomy=silva.nr_v123.tax, cutoff=80)
cluster.split(fasta=current, taxonomy=current, count=current, method=dgc)
make.shared(list=current, count=current, label=0.03) 
classify.otu(list=current, count=current, taxonomy=current, label=0.03)

get.current()

#OTU/Phylotype
phylotype(taxonomy=current)
#shared, name, group, processors, large, groups, seed
make.shared(list=current, count=current, label=1)
classify.otu(list=current, count=current, taxonomy=current, label=1)

get.current()

#OTU/Phylogenetic
dist.seqs(fasta=current, output=lt)
clearcut(phylip=current)

#OTU analysis
count.groups(shared=current)
sub.sample(shared=current)

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

set.logfile(name=mothur.shared)

# OTU based analysis (no phylogeny involved)
#
###################################

##charts for alpha diversity accuracy
#alpha = richness/diversity
#richness = # OTUs
#diversity = evenness of distribution of OTUs

#generate plotable stats for Chao1 and invsimpson for each sample, 100 otu intervals
#chao2 is measure of richness
#invsimpson = diversity
collect.single(shared=current,calc=chao-invsimpson,freq=100)

#generate rarefaction data using sobs (observations) 100 otu intervals
rarefaction.single(shared=current,calc=sobs,freq=100)

#summary single ##double check notes for which calcs are useful
summary.single(shared=current, calc=nseqs-coverage-sobs-invsimpson, subsample=T)

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

get.current()

#thetayc tree is used for following analysis
#summary.shared() makes an analysis based on declared groups: not really useful for batch analysis
#summary.shared(calc=sharedchao,groups=XXX-XYY-XYX)
#venn(groups=XXX-XYY-XYX)
#probably requires groups to work.

tree.shared(phylip=current)

#assume group file = "experimental.design"
#let's look at the likelihood that parsimony is similar/diff between groups:
parsimony(tree=current,group=experimental.design,groups=all)

#what about using unifrac weighting to tree distance ?
unifrac.weighted(tree=current,group=experimental.design,groups=all)

#let's look at PCOA for the thetayc output
#we probably won't like it compared to nmds:

pcoa(phylip=current)

#Nmds using 2-4 dimensions:
nmds(phylip=current, mindim=2, maxdim=4)

#what about the ACTUAL beta diversity of groups?  
#amova tests difference between groups vs. variance of groups using "centroid" of each groups' OTU differences
amova(phylip=current, design=experimental.design)

#what about the variance itself? Does it change between groups?
homova(phylip=current, design=experimental.design)

#let's find axes responsible for variation for OTUs responsible for most variance along X axes
#can change numaxes depending on structure

#if metatadata, include metadata=metadata.txt
corr.axes(axes=stability.trim.contigs.good.unique.good.filter.unique.pick.good.filter.unique.pick.precluster.pick.nr_v123.wang.tx.1.subsample.thetayc.1.lt.ave.nmds.axes, shared=current, method=spearman, numaxes=3)

#let's predict how many communities the data indicate, without groups:
get.communitytype(shared=current)
#analyze the relabund file to see if it matches our expectaions


#metastats (non-parametric t-test)
metastats(shared=current,design=experimental.design)

#lefse
lefse(shared=current,design=experimental.design)

#indicator analysis 
indicator(shared=current, design=experimental.design)

#classify.rf = random forest to find discriminatory OTUs
classify.rf(shared=current, design=experimental.design)


######################################
#TODO: Redo all analyses using phylotypes
#
######################################

###################
#Begin Phylogeny level analysis
#
##################
