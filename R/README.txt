#information about R scripts in this directory.


R scripts are written to work in a directory with the following files:

driver.R (to be run first)
MGSAT.0.03.biom
MGSAT.cons.taxonomy
MGSAT.cons.tax.summary
metadata.txt

phyloseq.R
MGSAT.0.03.biom
MGSAT.nmds.axes.nmds (generated from nmds output from mothur, followed by parsing using perl scrip in scripts dir)
MGSAT.0.03.thetayc.tre
Rank4.txt (can be a copy of the Rank4 analysis from MGSAT, or a list of Families to examine, one entry per line)
Rank5.txt (same as above, at genus level)

Each of these files is generated by mothur (biom file from make.biom.m), or comes from MGSAT.

Running these two scripts should generate standard report from MGSAT, followed by 16sCharts.pdf containing a number of useful charts for analysis from phyloseq outputs.

