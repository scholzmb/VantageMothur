#######################
# written by Matthew Scholz
# take outputs from batch_mothur_shared.m and use them for analysis


#first read in variables from last logfile
load.logfile(name=mothur.shared)
set.logfile(name=mothur.analysis)


##charts for alpha diversity accuracy

#generate plotable stats for Chao1 and invsimpson for each sample, 100 otu intervals
collect.single(shared=current,calc=chao-invsimpson,freq=100)

#generate rarefaction data using sobs (observations) 100 otu intervals
rarefaction.single(stability=current,calc=sobs,freq=100)

#summary single ##double check notes for which calcs are useful
summary.single(shared=stability.an.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=T)

##End alpha diversity

##Beta Diversity
#Beta diversity requires table of samples and conditions (example from MiSeq SOP data)

#format: SampleName\tCondition
#for additional groupings,create separate files
#fancy bashy option, create spreadsheet, saved as csv with all metadata 1/column, for each analysis:
#for i in {2..(last column of data}; do cut ....; done

#pretty pictures
#determine num OTUs to use before running
heatmap.bin(shared=current,scale=log2,numotu=50)

#distribution measurements (using thetayc)
dist.shared(shared=current, calc=thetayc, subsample=T)
heatmap.sim(phylip=current)
dist.shared(shared=current, calc=jclass, subsample=T)
heatmap.sim(phylip=current)


