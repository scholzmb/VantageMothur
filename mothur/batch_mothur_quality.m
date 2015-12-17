##############
#
# This is a batch file for running through the quality filtering steps of mothur
#
#  Original Script/comments written by by Dr. Tracy Teal

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
set.logfile(name=mothur.quality)


# Make contigs
# Join the paired ends
make.contigs(file=stability.files, processors=8)

# Screen seqs
# Get rid of sequences that are the wrong length or have ambiguous base pairs.
# We're limiting the length to 275 because we expect our sequences to be about 251 bp
# and this gives some flexibility
screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=275, processors=8)

# Get only the unique sequences
# This reduces our data set size for further processing
# No sequences are lost during this step. We're just generating a FASTA file of the unique
# sequences and keeping track of their numbers in each sample in a different file.
unique.seqs(fasta=current)

# Simplify the names and groups files 
#Generate new count file using uniques from previous command
count.seqs(name=current, group=current)

# We're assuming that we already made the pcr file that's the appropriate size
# for our amplicons.
# Here it's been created, is in this directory and is called silva.v4.fasta
# Align our sequences to that reference
# This is one of the most computationally intensive steps in this file
# reference can be changed.
align.seqs(fasta=current, reference=silva.bacteria.fasta, flip=t)


# Summarize the alignment information
# We'll need the output of this in the next step
summary.seqs(fasta=current, count=current)

# Get rid of sequences that don't align well
# ASSUMES V4 PCR
screen.seqs(fasta=current, count=current, summary=current, start=13862, end=23444, maxhomop=8)

# To reduce dataset size, get rid of columns with no information
# Get rid of columns in the file that only contain gaps '-' or overhang '.'
# MAKES FASTA uniform, deletes un-useful data, if screen is done wrong, will result in size 0 files
filter.seqs(fasta=current, vertical=T, trump=.)

# Now that we've aligned the sequences, we might see more that are identical
# Pare the dataset down to just the unique aligned ones
# NEED TO COUNT again
unique.seqs(fasta=current, count=current)
count.seqs(name=current,group=current)


# Pre-cluster to reduce dataset size
pre.cluster(fasta=current, count=current, diffs=2)

# Check for chimeras and remove them
# We're doing the chimera checking using the abundant reads as the reference
chimera.uchime(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current)

# Classify the sequences and remove ones that aren't bacterial
# We're using the Silva training sets as references: trainset9_032013
classify.seqs(fasta=current, count=current, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax, cutoff=80)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

# Take out the Mock sample from the dataset
# remove.groups(fasta=current, count=current, taxonomy=current, groups=Mock)

# Now we have a good set of quality sequences!

# Print out the last things we used so we know the last set of files generated.
get.current()

# We'll copy these to a shorter filename and use them in the next steps
# We will need the *.fasta, *.taxonomy, *.count_table files
# e.g
#accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos
#fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta
#group=stability.contigs.good.groups
#name=stability.trim.contigs.good.names
#taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy
#count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table
#processors=8
#summary=stability.trim.contigs.good.unique.summary
system(mv stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta quality_sequences.fasta)
system(mv stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy quality_sequences.taxonomy)
system(mv stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table quality_sequences.count_table)
system(mv stability.contigs.good.groups quality_sequences.groups)

quit()
