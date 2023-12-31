
#run this code from it's own folder

tfleaves = TRUE # to remove excess leaves
tfpoly = TRUE # to shrink polytomies
maxPoly = 4 # max polytomy size
maxLeaf = 2 # maximium leaf nodes for each parent (may exceed value based on number of annotated leaves across different functions)
functionNum = c(1,2)
z1 = readRDS('../treeSets/3pthr4treesEachwotaxon.rds')



library(geese)
source('../functions/shrinkPolytomies.r')
source('../functions/removeExcessLeafPolytomies.r')


#x = readRDS('AphyloTreePTHR42908-4Func.rds')
ruleLim = FALSE # set true to use rule limit changes
set.seed(53)


#### this set is not useful atm, more for tracking names when doing multiple sets
res = list()
c0 = 1
res[[c0]] = list()
res[[c0]]$name = c()
res[[c0]]$pthr = c()
for(i in functionNum){
    res[[c0]]$name = c(res[[c0]]$name,colnames(z1[[i]]$tip.annotation)[1])
    res[[c0]]$pthr = c(res[[c0]]$pthr,names(z1)[i])
}


annotations = data.frame()

for(i in functionNum){
	annotations = rbind(annotations,c(z1[[i]]$tip.annotation,rep(9,length(z1[[i]]$node.type))))
}
annotations = t(annotations)

rownames(annotations) = 1:length(annotations[,1])
annos = asplit(annotations,1)


edges = z1[[functionNum[1]]]$tree$edge
rootpos = 1 + length(z1[[functionNum[1]]]$tip.annotation)
#annotations = as.list(c(rep(9,length(z1[[1]]$node.type)),z1[[1]]$tip.annotation))
dups = c(rep(FALSE,rootpos -1),z1[[functionNum[1]]]$node.type == 0)
#annotations <- replicate(length(annotations), c(9, 9), simplify = FALSE)
print(rootpos)
z1 = z1[functionNum]
if(tfleaves) {
    z2 = removeLeaves(z1,2)
    z1 = z2$tree
    edges = z1[[1]]$tree$edge
    annos = z2$annos
    rootpos = min(edges[,1])
    print(rootpos)
    dups = c(rep(FALSE,rootpos -1),z1[[1]]$node.type == 0)
}


if(tfpoly && max(table(edges[,1])) > maxPoly) {
	print('start polytomie clean up')
	q = shrinkPolytomies(z1[[functionNum[1]]],annos,maxPoly)
	print('end polytomy cleanup')
}

print(max(table(q$tree$tree$edge[,1])))
edges = q$tree$tree$edge
annos = q$annos
dups = q$dups

#plottingGeese1(8,z2,1)
#plottingGeese2(8,z2,1)

rootpos = min(edges[,1])

bmodel <- new_geese(
  annotations = annos,
  geneid = c(edges[, 2], rootpos),
  parent = c(edges[, 1], -1),
  duplication = dups
  )
#parse_polytomies(bmodel)

# Adding the model terms
term_gains(bmodel,0:1,duplication = 0)
term_gains(bmodel,0:1,duplication = 1)
term_loss(bmodel,0:1,duplication = 0)
term_loss(bmodel,0:1,duplication = 1)
#term_maxfuns(bmodel, 0, 1,duplication = 2) #TODO but try seperate with 1 and 2 later


if(ruleLim == TRUE){
 k = 2
rule_limit_changes(bmodel, 0, 0, k)
rule_limit_changes(bmodel, 1, 0, k)
rule_limit_changes(bmodel, 2, 0, k)
rule_limit_changes(bmodel, 3, 0, k)
rule_limit_changes(bmodel, 4, 0, k)
rule_limit_changes(bmodel, 5, 0, k)
rule_limit_changes(bmodel, 6, 0, k)
rule_limit_changes(bmodel, 7, 0, k)
}

#rule_limit_changes(bmodel, 0, 0, 2, TRUE)
#rule_limit_changes(bmodel, 1, 0, 2, FALSE)

# We need to initialize to do all the accounting
init_model(bmodel)



res[[c0]]$mcmc <- geese_mcmc(
  bmodel,
  nsteps  = 20000,
  kernel  = fmcmc::kernel_ram(warmup = 2000), 
  prior   = function(p) dlogis(p, scale = 2, log = TRUE)
  )
mcmc = res[[c0]]$mcmc 
par_estimates <- colMeans(window(mcmc, start = 15000))

#params[[maxPoly]] = par_estimates
#saveRDS(params,'paramSetgeese3to7set.rds')