
#this is used for looking at a single panther tree with multiple functions
library(geese)
x = readRDS('../treeSelection/TestTreePTHR11848.rds')

set.seed(31)
#z1 = x[2:5]
z1 = x

annotations = data.frame()
for(i in c(1,2)){
	annotations = rbind(annotations,c(rep(9,length(z1[[1]]$node.type)),z1[[1]]$tip.annotation))
}
annos = asplit(t(annotations),1)

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
term_overall_changes(bmodel, TRUE)


 
rule_limit_changes(bmodel, 0, 0, 3)
rule_limit_changes(bmodel, 1, 0, 3)
rule_limit_changes(bmodel, 2, 0, 3)
rule_limit_changes(bmodel, 3, 0, 3)
rule_limit_changes(bmodel, 4, 0, 3)
rule_limit_changes(bmodel, 5, 0, 3)
rule_limit_changes(bmodel, 6, 0, 3)
rule_limit_changes(bmodel, 7, 0, 3)


#rule_limit_changes(bmodel, 0, 0, 2, TRUE)
#rule_limit_changes(bmodel, 1, 0, 2, FALSE)

# We need to initialize to do all the accounting
init_model(bmodel)
#> Initializing nodes in Geese (this could take a while)...
#> ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| done.

print(bmodel)
geese_mle(bmodel, hessian = TRUE)
 



