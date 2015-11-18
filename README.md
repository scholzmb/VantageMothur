

Sample Mothur workflows

git init
git remote add origin -t <branchname> git@github.com:scholzmb/VantageMothur.git
git pull
git checkout <branchname>


mothur/ directory holds batch.m files for varying steps.

order:
mothur_batch_quality.m
mothur_batch_shared.m
mothur_batch_phylotypes.m

sample sheets:
stability.files (sample name\tread files)
sampleTable.txt (sample name\tSampleType)


For alignment, you will need to symlink/copy Silva.v4. and silva.bacteria to mothur/
