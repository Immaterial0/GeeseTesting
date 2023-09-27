
if(!exists('flock4of4')) flock4of4 = readRDS('./4setsof4parestimatesgeese.rds')

if(FALSE){ Comment - as r designers wont add multi line comments
Gains 0 at speciation  -3.66793652
Gains 1 at speciation  -0.06226243
Gains 0 at duplication -0.39519106
Gains 1 at duplication -0.10060788
Loss 0 at speciation   -4.05708299
Loss 1 at speciation    1.36502243
Loss 0 at duplication   0.77646939
Loss 1 at duplication  -0.51967929
Root 1                  0.20829889
Root 2                  0.87044720
Root 3                  0.29086585
Root 4                 -0.18299463  
}



probsum = 0
res = c()
bitlist = data.frame()
for(j in 0:(2^1-1)){
bitsa = as.numeric(intToBits(j)[1:2])
dfa =  matrix(bitsa, ncol=1, byrow=TRUE)
for(i in 0:(2^4-1)){
 bits = as.numeric(intToBits(i)[1:4])
 df = matrix(bits, ncol=2, byrow=TRUE)
 q = transition_prob(bmodel, flock4of4[c(1:8)],FALSE,dfa,df)
 res = c(res,q)
 probsum = probsum + q
 bitlist = rbind(bitlist,c(bitsa,bits))
}
}
resdf = data.frame(bitlist, res)
