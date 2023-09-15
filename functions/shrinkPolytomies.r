
#method 1 : convert polytomies into polytomies of size s max
shrinkPolytomies = function(tree,annoLs,s){ # data frame with 2 columns, vector, integer
	terminaltoParent = TRUE # make this false if don't want it, will sort edges so add those with lower values first (hence the leaf nodes) to the original node
	edgeLs = tree$tree$edge
	#get list of all polytomies that are greater than s
	resObj = list()
	resObj$tree$edge = edgeLs
	resObj$annos = annoLs
	if(s < 2){
		Print('Polytome size constraint too small, skipping split')
		return(tree)
	}
 
	rootpos = 1 + length(tree$tip.annotation)
	tab = table(edgeLs[,1])
	n = names(tab)[tab > s]
	#print('check1')
	if(length(n) == 0){
		return(resObj)
	}
	#for each node get the list of children
	origMax = length(annoLs)
	maxNode = length(annoLs)
	newNodes = 0
	toDel = c()
	newEdges = data.frame()
	counter = 0
	countCounter = 0
	#print('check2')
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
		if(terminaltoParent){
			e = sort(e)
		}
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


	#print('check3')
	newEdges = as.matrix(newEdges)
	mode(newEdges) = 'integer'
	newEdges = unname(newEdges)
	#print(head(newEdges))
	#print(head(edgeLs))
	#print(colnames(newEdges))
	#print(colnames(edgeLs))

	if(!is.null(colnames(edgeLs))) colnames(newEdges) = colnames(edgeLs)

	#get branch distance from each child to parent (of ones being altered)
	#then look through each set that share a parent and get smallest distance
	#then subtract half that distance from all in set
	#then add that half distance from parent to upper parent
	#update any with parent greater than maxNode

	#find the branches for already existing nodes and see what original lengths were
	# then find the ones with parents that are greater than max node (hence are new nodes)
	# find min value in this set and subtract half from all new branches going to that node (to put new node halfway between the closest node to parent)
	# add that half from new node to its parent
	medge = match(newEdges[,2],edgeLs[,2])
	mbranch = tree$tree$edge.length[medge]
	for(i in 1:length(mbranch)){
		if(is.na(mbranch[i])) mbranch[i] = 0
	}
	mbranchparent = newEdges[,1]
	tf = mbranchparent > origMax
	unq = unique(mbranchparent[tf])
	minvallist = c()
	for(i in unq){
		w = which(mbranchparent == i)
		minval = min(mbranch[w]) / 2
		mbranch[w] = mbranch[w] - minval
		minvallist = c(minvallist,minval)
		
	}
	#print('check4')
	m = match(unq,newEdges[,2])
	mbranch[m] = minvallist
	edgeLs = edgeLs[-toDel,]
	tree$tree$edge.length = tree$tree$edge.length[-toDel]
	edgeNew = rbind(edgeLs,newEdges)
	for(i in (1+length(annoLs)):(length(annoLs) + newNodes)) annoLs[[i]] = rep(9,length(annoLs[[1]]))
	w = which(edgeNew[,2] > maxNode)
	parentNew = edgeNew[w,2]
	tree$node.type = c(tree$node.type,tree$node.type[parentNew - maxNode])
	tree$tree$node.label = c(tree$tree$node.label,rep('New',newNodes))
	tree$tree$edge.length = c(tree$tree$edge.length,mbranch)
	tree$tree$edge = edgeNew
	#resObj$edges = edgeNew
	resObj$annos = annoLs
	resObj$tree = tree
	if(sMult > s) resObj = shrinkPolytomies(resObj$tree,resObj$annos,s)
	return(resObj)
}


