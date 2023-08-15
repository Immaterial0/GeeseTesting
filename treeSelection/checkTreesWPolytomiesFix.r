library(geese)
library(stringr)



#still need to use to update dup list

#method 1 : convert polytomies into polytomies of size s max
shrinkPolytomies = function(edgeLs,annoLs,s){ # data frame with 2 columns, vector, integer


	#get list of all polytomies that are greater than s
	resObj = list()
	resObj$edges = edgeLs
	resObj$annos = annoLs
	if(s < 2){
		Print('Polytome size constraint too small, skipping split')
		return(resObj)
	}
 
	
	tab = table(edgeLs[,1])
	n = names(tab)[tab > s]
	print('check1')
	if(length(n) == 0){
		return(resObj)
	}
	#for each node get the list of children
	maxNode = length(annoLs)
	newNodes = 0
	toDel = c()
	newEdges = data.frame()
	counter = 0
	countCounter = 0
	for(nV in n){
		#print(counter)
		#print(nV)
		counter = counter + 1
		ntf = edgeLs[,1] == nV
		e = edgeLs[ntf,2]
		#print(length(e))
		toDel = c(toDel,c(1:length(ntf))[ntf])
		#get the size of polytomies relative to s to find number of new nodes to create
		sMult = length(e) / s
		sMult = floor(sMult)
		if(length(e) - s*sMult > s-sMult) sMult = sMult + 1
		newNodes = newNodes + sMult
		#assign most possible to parent node while limiting size of polytomies to s
		if(sMult > s) print('Recursive')
		numOnParent = s-sMult
		count = 0 
		if(numOnParent > 0){
			for(i in 1:numOnParent){
				count = count + 1
				if(is.na(e[count])) print(paste(nV,count))
				newEdges = rbind(newEdges,c(nV,e[count]))
			}
		}
		# equalize the number of children on all new nodes
		sMultNodes = rep(0,length(sMult))
		for(i in 1:sMult){
			maxNode = maxNode + 1
			newEdges = rbind(newEdges,c(nV,maxNode))
			sMultNodes[i] = maxNode
		}

		while(count < length(e)){
			for(i in 1:sMult){
				count = count + 1
				if(count > length(e)) break
				if(is.na(e[count])) print(paste(nV,count,'type2'))
				newEdges = rbind(newEdges,c(sMultNodes[i],e[count]))
			}
		}
		#countCounter = countCounter + count
	}
	#delete all previous children
	#print(toDel)
	edgeLs = edgeLs[-toDel,]

	newEdges = as.matrix(newEdges)
	mode(newEdges) = 'integer'
	newEdges = unname(newEdges)
	#print(head(newEdges))
	#print(head(edgeLs))
	#print(colnames(newEdges))
	#print(colnames(edgeLs))

	if(!is.null(colnames(edgeLs))) colnames(newEdges) = colnames(edgeLs)

	edgeNew = rbind(edgeLs,newEdges)
	for(i in (1+length(annoLs)):(length(annoLs) + newNodes)) annoLs[[i]] = rep(9,length(annoLs[[1]]))
	#annoNew = c(annoLs,rep(9,newNodes)) #assumed 9 is the unknown anno
	resObj$edges = edgeNew
	resObj$annos = annoLs

	if(sMult > s) resObj = shrinkPolytomies(resObj$edges,resObj$annos,s)
	
	
	return(resObj)


}













x = readRDS('3pthr4treesEachwotaxon.rds')
ruleLim = FALSE # set true to use rule limit changes
set.seed(31)
#z1 = x[2:5]
z1 = x
 
for(ca in 0:2){

res = list()
c0 = 1
c1 = 1
c2 = 2

for(c1 in (1+ca*4):((ca+1)*4-1)){
for(c2 in (c1+1):((ca+1)*4)){

print(c0)
res[[c0]] = list()
res[[c0]]$name1 = colnames(z1[[c1]]$tip.annotation)[1]
res[[c0]]$name2 = colnames(z1[[c2]]$tip.annotation)[1]
res[[c0]]$pthr1 = names(z1)[c1]
res[[c0]]$pthr2 = names(z1)[c2]

annotations = data.frame()

for(i in c(c1,c2)){
	annotations = rbind(annotations,c(z1[[i]]$tip.annotation,rep(9,length(z1[[i]]$node.type))))
}
annos = asplit(t(annotations),1)

edges = z1[[c1]]$tree$edge
rootpos = 1 + length(z1[[c1]]$tip.annotation)
#annotations = as.list(c(rep(9,length(z1[[1]]$node.type)),z1[[1]]$tip.annotation))
dups = c(rep(FALSE,rootpos -1),z1[[c1]]$node.type == 0)
#annotations <- replicate(length(annotations), c(9, 9), simplify = FALSE)

if(max(table(edges[,1])) > 8) {
	print('start polytomie clean up')
	q = shrinkPolytomies(edges,annos,8)
	print('end polytomy cleanup')
}
edges = q$edges
annos = q$annos

bmodel <- new_geese(
  annotations = annos,
  geneid = c(edges[, 2], rootpos),
  parent = c(edges[, 1], -1),
  duplication = c(dups,rep(FALSE,length(annos)-length(dups)))
  )


# Adding the model terms
term_gains(bmodel,0:1,duplication = 0)
term_gains(bmodel,0:1,duplication = 1)
term_loss(bmodel,0:1,duplication = 0)
term_loss(bmodel,0:1,duplication = 1)
term_maxfuns(bmodel, 0, 1,duplication = 2) #TODO but try seperate with 1 and 2 later
term_overall_changes(bmodel, TRUE)

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


init_model(bmodel)
#print(bmodel)



res[[c0]]$mle = geese_mle(bmodel, hessian = TRUE)



c0 = c0 + 1



}}


res2 = data.frame()
for(i in 1:length(res)){
	res2 = rbind(res2,c(res[[i]]$pthr1,res[[i]]$pthr2,res[[i]]$name1,res[[i]]$name2,res[[i]]$mle$par))
}
colnames(res2)[5:(dim(res2)[2]-2)] = names(bmodel)
 write.csv(res2,str_c('testSetGeese-',res2[[1]]$pthr1,'-all2FuncPairs.csv'))

}


for(i in 1:length(z1)){
print(colnames(z1[[i]]$tip.annotation)[1])
}

