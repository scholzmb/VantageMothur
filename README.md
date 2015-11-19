

Sample Mothur workflows

MAKE SURE YOUR GLOBAL CONFIG IS SET TO ONLY UPDATE CURRENT BRANCHES DURING PUSHES
git config --global push.default=current


git init
git remote add origin master git@github.com:scholzmb/VantageMothur.git
git pull
git branch <branchname>
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
