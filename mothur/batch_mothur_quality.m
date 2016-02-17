#Mothur first steps assembled 2/8/2016 by matthew scholz

#assemble/qc/qa/count
make.contigs(file=stability.files, processors=8)
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275)
unique.seqs(fasta=current)
count.seqs(name=current, group=current)
summary.seqs(count=current)

#skip this because we should already have it
#pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F)
#system(mv silva.bacteria.pcr.fasta silva.v4.fasta)

#align to silva
align.seqs(fasta=current, reference=silva.v4.fasta)
summary.seqs(fasta=current, count=current)
screen.seqs(fasta=current, count=current, summary=current, start=1968, end=11550, maxhomop=8)
summary.seqs(fasta=current, count=current)
filter.seqs(fasta=current, vertical=T, trump=.)
unique.seqs(fasta=current, count=current)

#cluster
pre.cluster(fasta=current, count=current, diffs=2)
chimera.uchime(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current)
summary.seqs(fasta=current, count=current)

#classify/remove useless
classify.seqs(fasta=current, count=current, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax, cutoff=80)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#OTU Analysis
cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.15)
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
dist.seqs(fasta=current, output=lt, processors=8)
clearcut(phylip=current)

#OTU analysis
count.groups(shared=current)
sub.sample(shared=current)

