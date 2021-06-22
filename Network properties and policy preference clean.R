rm(list = ls())
setwd("~/Google Drive/Dissertation/Analysis_network properties and policy preference")
library(igraph)
nodes <- read.csv("Node.csv", header=T, as.is=T, row.names=1)
matrix <- read.csv("Matrix_bi_monthly.csv", header=T, row.names=1)
head(nodes)

head(matrix)
matrix <- as.matrix(matrix)
dim(matrix)
net <- graph_from_incidence_matrix(matrix)
nodes <- nodes[names(V(net)),]
#links <- as_edgelist(net_temp, names = TRUE)
#head(links)
#net <- graph_from_data_frame(d=links, vertices=nodes, directed=FALSE) #use this to add attribute
E(net) #edge
V(net) #vertices
vertex_attr(net) <- c(vertex_attr(net), as.list(nodes))
vertex_attr(net) #check attribute
V(net)$type
#names(V(net))[!names(V(net)) %in% nodes$id]
##try plot the graph
V(net)$color <- c("lightsteelblue2", "tomato")[V(net)$type+1]
V(net)$shape <- c("circle", "square")[V(net)$type+1]
V(net)$size	<- 0.1*degree(net) + 3
#coord = layout_with_graphopt(net, charge = 1e-1)
plot((delete.vertices(net, degree(net)==0)), 
     vertex.label=NA, layout=layout_with_kk)

##generate static model
library(ergm)
matrix <- read.csv("Matrix_bi_monthly.csv", header=T, row.names=1)
matrix <- as.matrix(matrix)
net_ergm <- network(matrix, directed=FALSE, loops=FALSE, bipartite=TRUE, matrix.type="bipartite")
class(net_ergm)
nodeInfo <- read.csv("Node_attr.csv",header=TRUE,stringsAsFactors=FALSE)
net_ergm
network.density(net_ergm, na.omit=TRUE, discount.bipartite=TRUE)
network.vertex.names(net_ergm)

nodeInfo
nodeInfo$Q13H
net_ergm%v%"sector" <- nodeInfo$sector
#net_ergm%v%"Ground" <- nodeInfo$Ground
net_ergm%v%"Q13A" <- nodeInfo$Q13A
net_ergm%v%"Q13B" <- nodeInfo$Q13B
net_ergm%v%"Q13C" <- nodeInfo$Q13C
net_ergm%v%"Q13D" <- nodeInfo$Q13D
net_ergm%v%"Q13E" <- nodeInfo$Q13E
net_ergm%v%"Q13F" <- nodeInfo$Q13F
net_ergm%v%"Q13G" <- nodeInfo$Q13G
net_ergm%v%"Q13H" <- nodeInfo$Q13H
net_ergm%v%"Q13I" <- nodeInfo$Q13I
net_ergm%v%"Q13J" <- nodeInfo$Q13J
net_ergm%v%"Q13K" <- nodeInfo$Q13K
net_ergm%v%"Q13L" <- nodeInfo$Q13L
net_ergm%v%"Q13M" <- nodeInfo$Q13M
net_ergm%v%"Q13_A14" <- nodeInfo$Q13_A14
net_ergm%v%"Q13_A15" <- nodeInfo$Q13_A15
net_ergm%v%"Q13_A16" <- nodeInfo$Q13_A16
get.vertex.attribute(net_ergm, 'Q13A') ##this command has problem now
list.vertex.attributes(net_ergm) ##this command has problem now
get.network.attribute(net_ergm, 'bipartite', unlist = FALSE)

net_ergm%v%"Q13D"
###determine the parameters
set.seed(111)
MLE.alpha = numeric(11)
MLE.beta = numeric(11)
for (i in 0:10){
  print(i)
  alpha = i/10
  beta = i/10
  model.alpha = try(ergm(net_ergm ~ edges + b1nodematch('Q13A', alpha=alpha)))
  model.beta = try(ergm(net_ergm ~ edges + b1nodematch('Q13A', beta=beta)))
  MLE.alpha[i+1] = model.alpha$mle.lik
  MLE.beta[i+1] = model.beta$mle.lik
}
MLE.beta
MLE.alpha
###picture
depth = seq(0,1,0.1)
plot(depth, MLE.beta, pch=19, type='b', col='red', xlab='alpha, beta', ylab='Loglikelood')
lines(depth, MLE.alpha, pch=15, type='b', col='green')
#legend
###end here

summary(net_ergm ~ edges + b1nodematch('sector', beta=0.5, diff=TRUE) + 
          b1nodematch('Q13A', beta=0.5) + 
          b1nodematch('Q13B', beta=0.5) + 
          b1nodematch('Q13C', beta=0.5) + 
          b1nodematch('Q13D', beta=0.5) + 
          b1nodematch('Q13E', beta=0.5) + 
          b1nodematch('Q13F', beta=0.5) + 
          b1nodematch('Q13G', beta=0.5) + 
          b1nodematch('Q13H', beta=0.5) + 
          b1nodematch('Q13I', beta=0.5) + 
          b1nodematch('Q13J', beta=0.5) + 
          b1nodematch('Q13K', beta=0.5) + 
          b1nodematch('Q13L', beta=0.5) + 
          b1nodematch('Q13M', beta=0.5) + 
          b1nodematch('Q13_A14', beta=0.5) + 
          b1nodematch('Q13_A15', beta=0.5) + 
          b1nodematch('Q13_A16', beta=0.5))
##model fit
set.seed(111)
model.attr <- ergm(net_ergm ~ edges + b1nodematch('sector', beta=0, diff=TRUE) + 
                     b1nodematch('Q13A', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13B', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13C', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13D', beta=0, diff=TRUE, levels=2) + 
                     b1nodematch('Q13E', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13F', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13G', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13H', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13I', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13J', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13K', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13L', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13M', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13_A14', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13_A15', beta=0, diff=TRUE, levels=3) + 
                     b1nodematch('Q13_A16', beta=0, diff=TRUE, levels=3))
summary(model.attr)
names(model.attr)
mcmc.diagnostics(model.attr)
model.attr$mle.lik

##gooddness of fit
#help(gof)
gof.deg <- gof(model.attr ~ model)
gof.deg
plot(gof.deg)
