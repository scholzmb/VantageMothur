## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)

load.meta.rsv <- function(file.name,batch=NULL,aggr.var=NULL) {
  
  meta =read.delim(file.name, header=TRUE,stringsAsFactors=T, sep="\t")
  make.global(meta)
  meta = meta[!duplicated(meta$SampleID),]
  
  

    
  meta$quantitative.quant = quantcut.ordered(meta$quantitative)
  
  
  row.names(meta) = meta$SampleID
  
  return (meta)
}


## This function should carry out analysis specific to metadata fields by themselves, without
## relation to the abundance profiles. You can write it to do nothing (empty body).

summary.meta.rsv <- function(m_a) {
  
  report$add.header("Summary of metadata variables")
  
  meta = m_a$attr
  
  report$add(summary(meta),caption="Summary of metadata variables")
  
  xtabs.formulas = list("~quantitative.quant+Treatment")
  for(xtabs.formula in xtabs.formulas) {
    fact.xtabs = xtabs(as.formula(xtabs.formula),data=meta,drop.unused.levels=T)
    report$add(fact.xtabs,caption=paste("Sample cross tabulation",xtabs.formula))
    report$add.printed(summary(fact.xtabs))
  }
    
}

## This function must generate a lits with analysis tasks

gen.tasks.rsv <- function() {
  
  task0 = within( mgsat.16s.task.template, {
  
  taxa.levels = c(2,3,4,5,6,"otu")
  
  descr = "All Samples"
  
  main.meta.var = "quantitative.quant"
  main.meta.var.cont = "quantitative"
  
  read.data.task = within(read.data.task, {
    taxa.summary.file = "MGSAT.cons.tax.summary"
    otu.shared.file="MGSAT.shared"
    cons.taxonomy.file="MGSAT.cons.taxonomy"
    meta.file="metadata.txt"
    load.meta.method=load.meta.rsv
    load.meta.options=list()    
  })
  
  summary.meta.method=summary.meta.rsv
  
  summary.meta.task = within(summary.meta.task, {
    meta.x.vars = c()
    group.vars = c(main.meta.var)
  })  
  
  })
  
task1 = within( task0, {
  
  do.summary.meta = T
  
  do.tests = T

  descr = "All samples, association with quantitative"
  
  test.counts.task = within(test.counts.task, {
    
    norm.count.task = within(norm.count.task, {
      method="norm.ihs.prop"
    })
    
    do.deseq2 = T
    do.adonis = T
    do.genesel = T
    do.stabsel = T
    do.glmer = F
    #do.divrich = c()
    
    do.plot.profiles.abund=T
    do.heatmap.abund=T
    
    divrich.task = within(divrich.task,{
      #n.rar.rep=4
      is.raw.count.data=T
      group.attr = main.meta.var
      counts.glm.task = within(counts.glm.task,{
        formula.rhs = main.meta.var.cont
      })      
    })
    
    deseq2.task = within(deseq2.task, {
      formula.rhs = main.meta.var.cont
    })
    
    genesel.task = within(genesel.task, {
      group.attr = main.meta.var
    })
    
    stabsel.task = within(stabsel.task, {
      resp.attr=main.meta.var.cont
      args.fitfun = within(args.fitfun, {
        family="gaussian"
      })
    })
    
    adonis.task = within(adonis.task, {
      
      tasks = list(
        list(formula.rhs=main.meta.var.cont,
             strata=NULL,
             descr=sprintf("Association with %s",main.meta.var.cont))        
      )
    })
        
    plot.profiles.task = within(plot.profiles.task, {
      id.vars.list = list(c(main.meta.var))
      clade.meta.x.vars=c()
      do.profile=T
      do.clade.meta=F
    })
    
    heatmap.abund.task = within(heatmap.abund.task,{
      attr.annot.names=c(main.meta.var.cont)
    })
    
  })
  
})



tasks = foreach(main.meta.var.loop = c("Treatment"
)) %do%
{
                  within( task0, {
  
  main.meta.var = main.meta.var.loop
                    
  do.summary.meta = F
  
  do.tests = T
  
  descr = paste("All samples, association with",main.meta.var)
  
  test.counts.task = within(test.counts.task, {
    
    norm.count.task = within(norm.count.task, {
      method="norm.ihs.prop"
    })
    
    do.deseq2 = T
    do.adonis = T
    do.genesel = T
    do.stabsel = T
    do.glmer = F
    #do.divrich = c()
    
    do.plot.profiles.abund=T
    do.heatmap.abund=T
    
    divrich.task = within(divrich.task,{
      #n.rar.rep=4
      is.raw.count.data=T
      group.attr = main.meta.var
      counts.glm.task = within(counts.glm.task,{
        formula.rhs = main.meta.var
      })      
    })
    
    deseq2.task = within(deseq2.task, {
      formula.rhs = main.meta.var
    })
    
    genesel.task = within(genesel.task, {
      group.attr = main.meta.var
    })
    
    stabsel.task = within(stabsel.task, {
      resp.attr=main.meta.var
      args.fitfun = within(args.fitfun, {
        family="binomial"
      })
    })
    
    adonis.task = within(adonis.task, {
      
      tasks = list(
        list(formula.rhs=main.meta.var,
             strata=NULL,
             descr=sprintf("Association with %s",main.meta.var))        
      )
    })
    
    plot.profiles.task = within(plot.profiles.task, {
      id.vars.list = list(c(main.meta.var))
      clade.meta.x.vars=c(main.meta.var.cont)
      do.profile=T
      do.clade.meta=T
    })
    
    heatmap.abund.task = within(heatmap.abund.task,{
      attr.annot.names=c(main.meta.var.cont,main.meta.var)
    })
    
  })
  
})

}

return (c(list(task1),tasks))
}




## number of cores to use on multicore machines
options(mc.cores=1)
options(boot.ncpus=1)
## parallel backend
options(boot.parallel="snow")
library("BiocParallel")
register(SnowParam(1))


## location of MGSAT code
MGSAT_SRC = "C:\\Users\\MBS\\Documents\\MGSAT"
## MGSAT_SRC = "/Users/shiltsmh/work/mgsat"

source(paste(MGSAT_SRC,"dependencies.r",sep="/"),local=T)

## Uncomment next line to install packages needed by MGSAT (!!!comment it out
## in all subsequent runs once the packages have been installed!!!).
## Note: you should also pre-install Pandoc program from http://johnmacfarlane.net/pandoc/
## or using your OS package manager (if running on Linux)

#install_required_packages()

## loads dependency packages (which already must be installed)
load_required_packages()

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=F)

## set incremental.save=T only for debugging or demonstration runs - it forces 
## report generation after adding every header section, thus slowing down
## a long run. But then incremental.save=T, you can open HTML report file in
## a Web browser and refresh it periodically to see it grow.
report <- PandocAT$new(author="matthew.b.scholz@vanderbilt.edu",
                       title="MGSAT practice",
                       incremental.save=F)


res = proc.project(
  task.generator.method=gen.tasks.rsv
)

report$save()
