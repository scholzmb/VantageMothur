

#Sample Mothur workflows

##MAKE SURE YOUR GLOBAL CONFIG IS SET TO ONLY UPDATE CURRENT BRANCHES DURING PUSHES

```bash
git config --global push.default=current
git init
git remote add origin master git@github.com:scholzmb/VantageMothur.git
git pull
git branch <branchname>
git checkout <branchname>
```

mothur/ directory holds batch.m files for varying steps.
run, in order 
1. mothur_batch_quality.m
2. mothur_batch_shared.m
3. mothur_batch_phylotypes.m


data/ contains sample sheets (also link.copy data here):
stability.files (sample name\tread files)
sampleTable.txt (sample name\tSampleType)

For alignment, you will need to symlink/copy Silva.v4. and silva.bacteria data/

run mothur batch files from data/
