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

# Set the current files to the ones we generated from the quality sequencing
# At the end we renamed our really long files to something shorter, so we
# can use those names here. Put in your own names if you did something different.
set.current(fasta=quality_sequences.fasta, count=quality_sequences.count_table, taxonomy=quality_sequences.taxonomy, processors=8, group=quality_sequences.groups)

# Calculate the distance matrix
dist.seqs(fasta=current, cutoff=0.2)

# Cluster the sequences
# We're using cluster.split so it can finish more quickly (i.e. in our lifetime)
cluster.split(column=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.15)

# Make the shared file for 0.03 distance
make.shared(list=current, count=current, label=0.03)

rarefaction.single(shared=current)

classify.otu(list=current,count=current,taxonomy=current,label=0.03)

phylotype(taxonomy=current)

make.shared(list=current, count=current,label=1)

classify.otus(list=current,count=current,taxonomy=current,label=1)
get.current()
