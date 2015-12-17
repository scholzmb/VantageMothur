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
phylotype(taxonomy=current)

#make genus level table for phylogeny
make.shared(list=current, count=current,label=1)

#classify them (names for OTUs)
classify.otu(list=current,count=current,taxonomy=current,label=1)

get.current()
