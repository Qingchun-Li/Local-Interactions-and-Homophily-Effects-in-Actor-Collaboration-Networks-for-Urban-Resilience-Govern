rm(list = ls())
setwd("~/Google Drive/Dissertation/Analysis_network properties and policy preference")
library(igraph)
library(statnet)
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

#calculate motif
library(bmotif)
motif_no <- mcount(matrix, normalisation=FALSE, mean_weight = FALSE, standard_dev = FALSE)
motif_position <- node_positions(matrix, weights_method = "none", weights_combine = "none")[1:197,1:16] 
## I choose the first 16 types of motif and 197 is the number of columns
motif_links <- link_positions(matrix, six_node = FALSE, weights = FALSE, normalisation = "none")[,2:9]
write.csv(motif_position, file = "motif_position.csv")
write.csv(motif_no, file = "motif_no.csv")
write.csv(motif_links, file = "motif_link.csv")
##generate static model
library(ergm)
matrix <- read.csv("Matrix_bi_monthly.csv", header=T, row.names=1)
matrix <- as.matrix(matrix)
net_ergm <- network(matrix, directed=FALSE, loops=FALSE, bipartite=TRUE, matrix.type="bipartite")
net_ergm ##generated network
network.density(net_ergm, na.omit=TRUE, discount.bipartite=TRUE)
network.vertex.names(net_ergm)

##try some other terms related to structure
#a <- net_ergm ~ edges + b1star(2) + b1star(3) + b2star(2) + b2star(3) + threetrail + cycle(4)
#summary(a)
# this results could test the above motif result


###generate random bipartite network with same degree distribution

##another package
library(networksis)
library(statnet)
sim.structure <- simulate.sisnetwork(net_ergm, nsim = 1001, seed=3, save.networks=TRUE)
mcount(as.matrix.network(sim.structure_test$networks), normalisation=FALSE, mean_weight = FALSE, standard_dev = FALSE)
summary(sim.structure)
names(sim.structure)

edges.no = numeric(999)
b1star2 = numeric(999)
b1star3 = numeric(999)
b2star2 = numeric(999)
b2star3 = numeric(999)
threetrail.no = numeric(999)
cycle.no = numeric(999)
for (i in 2:1000) {
  print(i)
  motif_sim <- mcount(as.matrix.network(sim.structure$networks[[i]]), normalisation=FALSE, mean_weight = FALSE, standard_dev = FALSE)
  edges.no[i-1] = motif_sim[[1,3:3]]
  b1star2[i-1] = motif_sim[[2,3:3]]
  b2star2[i-1] = motif_sim[[3,3:3]]
  b1star3[i-1] = motif_sim[[4,3:3]]
  threetrail.no[i-1] = motif_sim[[5,3:3]]
  cycle.no[i-1] = motif_sim[[6,3:3]]
  b2star3[i-1] = motif_sim[[7,3:3]]
}
library(psych)
mean(edges.no)
sd(edges.no)
mean(b1star2)
sd(b1star2)
mean(b2star2)
sd(b2star2)
mean(b1star3)
sd(b1star3)
mean(b2star3)
sd(b2star3)
mean(threetrail.no)
sd(threetrail.no)
mean(cycle.no)
sd(cycle.no)
edges.no
b1star2
b2star2
b1star3
b2star3
threetrail.no
cycle.no

fig_datafram <- data.frame(
  name <- c('Edges', 'Two stars:SR2', 'Two stars:SP2', 'Three stars:SR3', 'Three stars:SP3', 'Three trails', 'Cycles'),
  Mean_value <- c(mean(Bedges.no), mean(Bb1star2), mean(Bb2star2), mean(Bb1star3), mean(Bb2star3), mean(Bthreetrail.no), mean(Bcycle.no)),
  line <- c(1414, 13635, 18302, 231964, 126400, 430999, 49626)
)
Network_configuration  = factor(fig_datafram$name, levels=c('Edges', 'Two stars:SR2', 'Two stars:SP2', 'Three stars:SR3', 'Three stars:SP3', 'Three trails', 'Cycles'))
fig_datafram
motif_datafram
library(ggplot2)
# line
ggplot(fig_datafram) +
  geom_bar(aes(x=Network_configuration, y=Mean_value), stat="identity", size=1, color='darkblue', fill="white", alpha=0.5) + 
  geom_line(aes(x=Network_configuration, y = line), size = 1.5, color="red", group = 1) + 
  scale_y_continuous(labels = scales::label_comma(), limits=c(0, 550000)) + 
  theme(axis.text=element_text(size=13),
          axis.title=element_text(size=16,face="bold"))

  
###visualization to compare the number of network motif with generated random network
library("Hmisc")
library(ggplot2)
hist(cycle.no, prob=TRUE, breaks=12, col="#E5E5E5", border=0, 
     xlim=c(40000,50000))
lines(density(cycle.no), col="red", lwd=2)
lines(x=c(49626, 49626), y=c(0, 4e-4), col='darkred', lwd=2)


hist(threetrail.no, prob=TRUE, breaks=12, col="#E5E5E5", border=0)
lines(density(threetrail.no), col="red", lwd=2)
lines(x=c(629503, 629503), y=c(0, 8e-5), col='darkred', lwd=2)
