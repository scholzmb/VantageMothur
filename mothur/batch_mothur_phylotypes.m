##############
#
# batch_mothur_phylotypes.m
#
# Running Phylotype/Phylogenetic analysis from generic current inputs (results of system(mv) from batch_mothur_shared.m

#  Original Script/comments written by by Dr. Tracy Teal

#
# --------------------------------------------------------------------------------- 
#
# We will make use of the 'current' settings, so we don't have to type out the whole 
# file names in each step. This also reduces typo errors and ensures you're using
# the file that was created from the previous step.
#
# It does not include the summary.seqs commands and assumes that you have looked
# at your output to ensure that you have made the appropriate decisions for your
# quality filtering parameters. If you haven't gone through this process at least
# once and examined your output, do that first!
# 
##############

# Set the current files to the ones we generated from the quality sequencing
# At the end we renamed our really long files to something shorter, so we
# can use those names here. Put in your own names if you did something different.
#set.current(fasta=quality_sequences.fasta, count=quality_sequences.count_table, taxonomy=quality_sequences.taxonomy, processors=8, group=quality_sequences.groups)

load.logfile(mothur.1447253583.logfile)



