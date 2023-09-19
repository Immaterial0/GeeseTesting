
removeLeaves = function(trees,s=2){
    # s = max size of polytomies
    rootpos = min(trees[[1]]$tree$edge[,1])
    edges = trees[[1]]$tree$edge
    annotations = data.frame()
    for(i in 1:length(trees)){
        annotations = rbind(annotations,c(trees[[i]]$tip.annotation,rep(9,length(trees[[i]]$node.type))))
    }
    annotations = t(annotations)
    #print(dim(annotations))
    rownames(annotations) = 1:length(annotations[,1])
    annosNums = 1:length(annotations[,1])
    m = match(edges[,2],annosNums) # note edges is 1 less length than annos as doesn't include root in child terms
    edges1 = data.frame(c(1:length(edges[,1])   ),edges,rowSums(annotations[m,]),trees[[1]]$tree$edge.length) #trees is a fixed number here
    colnames(edges1) = c('num','parent','child','rowsums','edgelength')
    pret = edges1[edges1$child < rootpos,]
    t = table(pret$parent) # table of parents for which children are leaf nodes
    tvals = names(t)[t > s] # gives the names of parents that have more leaf nodes as children than s
    removeList = data.frame()
    evCount = 0
    for(i in 1:length(tvals)){
        tarr = pret[pret$parent == tvals[i] ,] # leaf node children of a given parent
        evRows = tarr$rowsums != 9*dim(annotations)[2] # finding out which leafs to keep due to being annotated
        if(sum(evRows) < s) dontRemove = s - sum(evRows)
        else dontRemove = 0
        #get the list of leaves to remove and add to a list based on shortest branch
        tarr2 = tarr[!evRows,] # getting leaf node children of given parent that aren't annotated for any function under consideration
        evCount = evCount + sum(evRows) # counter to keep track of children with evidence in the set of leaves of polytomies greater than size s
        if(dontRemove > 0){
            for(d in 1:dontRemove){
                tarr2 = tarr2[-which.min(tarr2$edgelength),] # keeping the shortest edge length children that aren't annotated
            }
        }
        removeList = rbind(removeList,tarr2)
    }
    # remome leaves from both edges and annotations
    # can use $num for edges, and child for annotations I think
    #is it correct not to keep unknowns if multiple evidence rows?
    #####################
    removedNum = edges[removeList$num,2]
    edges2 = edges[-removeList$num,]
    annotations2 = annotations[-removeList$child,]
    annos2 = asplit(annotations2,1)
    #check table of tvals at end make sure two children each

    # create dataframe linking 1:number of annos to the values still in the edge list, then replace all values
    fix1 = 1:length(annos2)
    fix2 = sort(unique(c(edges2[,1],edges2[,2])))
    fixdf = data.frame(fix1,fix2)
    m1 = match(edges2[,1],fix2)
    m2 = match(edges2[,2],fix2)
    edges2[,1] = fix1[m1]
    edges2[,2] = fix1[m2]
    edges = edges2
    annos = annos2
    for(i in 1:length(trees)){
        trees[[i]]$tree$edge = edges
        trees[[i]]$tree$edge.length = trees[[i]]$tree$edge.length[-removeList$num]
        trees[[i]]$tree$tip.label = trees[[i]]$tree$tip.label[-removedNum]
        trees[[i]]$tree$Nnode  = length(trees[[i]]$tree$tip.label) + length(trees[[i]]$tree$node.label)
    }
    resObj = list()
    resObj$annos = annos
    resObj$tree = trees
    return(resObj)
}