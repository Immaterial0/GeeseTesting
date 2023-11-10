#run this code from it's own folder

tfleaves = TRUE # doesn't work for flocks atm I think
tfpoly = TRUE
maxPoly = 4 # max polytomy size
maxLeaf = 2 # maximium leaf nodes for each parent (may exceed value based on number of annotated leaves across different functions)
stepnum = 20000
ruleLim = TRUE # set true to use rule limit changes
randomizeit = TRUE
timecounter = c()
scale = FALSE

set.seed(53)
set.seed(399292)
annosunderpoly = list()
library(geese)
source('../functions/shrinkPolytomies.r')
source('../functions/removeExcessLeafPolytomies.r')

randres = list()
#for(randomlist in 1:10){
randomlist = 3

annoslist = list()
annolistcount = 1
flock <- new_flock()
#treelocs = c('../treeSets/3pthr4treesEachwotaxon.rds','../treeSets/3pthr4treesEachwotaxon.rds','../treeSets/3pthr4treesEachwotaxon.rds','../treeSets/AphyloTreePTHR42908-4Func.rds')
treeNums = list()
treeNums[[1]] = c(seq(1,33,2))[1:randomlist]#,23)
treeNums[[2]] = treeNums[[1]]+1
#treeNums[[3]] = c(9:12)
#treeNums[[4]] = c(1:4)
#for(i in 1:(randomlist+1)){
 # treeNums[[i]] = c(i*2)
#}
treelocs = rep('../treesets/pthrPTHR23255_PTHR24416forflock33trees2p2nwotax.rds',length(treeNums))




for(flocknum in 1:(length(treelocs))){
#x = readRDS('AphyloTreePTHR42908-4Func.rds')
z1 = readRDS(treelocs[flocknum])


res = list()
c0 = flocknum
functionNum = treeNums[[flocknum]]

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

z1 = z1[functionNum]
if(tfleaves) {
    z2 = removeLeaves(z1,2)
    z1 = z2$tree
    edges = z1[[1]]$tree$edge
    annos = z2$annos
    rootpos = min(edges[,1])
    dups = c(rep(FALSE,rootpos -1),z1[[1]]$node.type == 0)
}



if(tfpoly && max(table(edges[,1])) > maxPoly) {
	#print('start polytomie clean up')
	q = shrinkPolytomies(z1[[1]],annos,maxPoly)
	#print('end polytomy cleanup')
    edges = q$tree$tree$edge
    annos = q$annos 
    dups = q$dups   
}


#plottingGeese1(8,z2,1)
#plottingGeese2(8,z2,1)

rootpos = min(z1[[1]]$tree$edge[,1])

print('annos info')
print(length(annos) == max(edges[,1]))
print(length(annos) == length(dups))
print(length(annos))

annoslist[[annolistcount]] = annos
annolistcount = annolistcount + 1

add_geese(
  flock,
  annotations = annos,
  geneid = c(edges[, 2], rootpos),
  parent = c(edges[, 1], -1),
  duplication = dups
  )
#parse_polytomies(bmodel)
}


# Adding the model terms
term_gains(flock,c(0:(length(treeNums[[1]])-1)),duplication = 2)
# term_gains(flock,0:c(0:length(treeNums[[1]])),duplication = 1)
term_loss(flock,c(0:(length(treeNums[[1]])-1)),duplication = 2)
#term_loss(flock,0:c(0:length(treeNums[[1]])),duplication = 1)
term_maxfuns(flock, 0, 1,duplication = 2) #TODO but try seperate with 1 and 2 later
# max funs penalizes if more than 1 function?

k = 2
if(ruleLim){
  for(rulcount in 0:(length(treeNums[[1]])*2)) rule_limit_changes(flock, rulcount, 0, k, duplication = 2) # should last one be labeled duplication?
}

#rule_limit_changes(bmodel, 0, 0, 2, TRUE)
#rule_limit_changes(bmodel, 1, 0, 2, FALSE)

# We need to initialize to do all the accounting
init_model(flock)

r1 = rep(0.1,nterms(flock)-nfuns(flock))
r2 = rep(0.5,nfuns(flock))


ans_mle = geese_mle(flock, hessian = TRUE, ncores = 4, control = list(maxit = 100))

time3 = system.time(mcmc <- geese_mcmc(
  flock,
  nsteps  = stepnum,
  kernel  = fmcmc::kernel_adapt(warmup = 1000,lb = -10, ub = 10), 
  prior   = function(p) dlogis(p, scale = 2, log = TRUE),#function(p) dlogis(p, scale = 2, log = TRUE),
  initial = ans_mle$par,
  ncores = 8
  )
) 
timecounter = c(timecounter,time3)
#kernal_ram was previous kernal being used


par_estimates2 <- colMeans(window(mcmc, start = 3/4*stepnum))



ans_pred <- predict_flock(
  flock, par_estimates2,
  leave_one_out = TRUE,
  only_annotated = TRUE
  ) |>
  lapply(do.call, what = "rbind") |>
  do.call(what = rbind)


ann_obs = data.frame()
for(i in 1:length(annoslist)){
  ann_obs <- rbind(ann_obs,
    do.call(rbind, annoslist[[i]]))
}
ann_obs = as.matrix(ann_obs)

library(aphylo)

(ans <- prediction_score(as.matrix(c(ans_pred)), as.matrix(c(ann_obs)))) #vectorizing seems to fix issue with predictions
auc = ans$auc$auc + (1 + ans$auc$tpr[1]) * (1 - ans$auc$fpr[1]) / 2
library(Metrics)
mae = mae(c(ans_pred)[!(c(ann_obs) == 9)],c(ann_obs)[!(c(ann_obs) == 9)])
auc2 = pROC::auc(c(ann_obs)[!(c(ann_obs) == 9)],c(ans_pred)[!(c(ann_obs) == 9)])

print(paste('auc from aphylo (corrected pre upgrade) : ',auc))
print(paste('mae from metrics library : ',mae))
print(paste('auc from pROC library : ',auc2))
#(ans <- prediction_score(as.matrix(ans_pred), as.matrix(ann_obs))) # this doesn't work, need to linearize it
 
randres[[randomlist]] = list()
randres[[randomlist]]$auc = auc
randres[[randomlist]]$mae = mae
randres[[randomlist]]$auc2 = auc2
randres[[randomlist]]$time = time3
 