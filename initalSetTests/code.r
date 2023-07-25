
#this is used for looking at a single panther tree with multiple functions

treeLoc = 'AphyloTreePTHR42908-4Func.rds'
z1 = readRDS(treLoc)
set.seed(31)

#add the annotations into a dataframe
#requires internal node annotations which are all set as 9 here
annotations = data.frame()
for(i in c(2,3)){
	annotations = rbind(annotations,c(rep(9,length(z1[[1]]$node.type)),z1[[1]]$tip.annotation))
}
annos = asplit(t(annotations),1) #transpose to correct orientation

edges = z1[[1]]$tree$edge
rootpos = 1 + length(z1[[1]]$tip.annotation)
#annotations = as.list(c(rep(9,length(z1[[1]]$node.type)),z1[[1]]$tip.annotation))
dups = c(rep(FALSE,rootpos -1) , z1[[1]]$node.type == 0)
#annotations <- replicate(length(annotations), c(9, 9), simplify = FALSE)

bmodel <- new_geese(
  annotations = annos,
  geneid = c(edges[, 2], rootpos),
  parent = c(edges[, 1], -1),
  duplication = dups
  )


# Adding the model terms
term_gains(bmodel,0:1,duplication = 0)
term_gains(bmodel,0:1,duplication = 1)
term_loss(bmodel,0:1,duplication = 0)
term_loss(bmodel,0:1,duplication = 1)
term_maxfuns(bmodel, 0, 1,duplication = 2) #TODO but try seperate with 1 and 2 later

#rule_limit_changes(bmodel, 0, 0, 2, TRUE)
#rule_limit_changes(bmodel, 1, 0, 2, FALSE)

init_model(bmodel)

print(bmodel)
geese_mle(bmodel, hessian = TRUE)
 


