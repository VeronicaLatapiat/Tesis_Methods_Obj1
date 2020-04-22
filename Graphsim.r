#install.packages("igraph")
#if(! require("devtools")){
#  install.packages("devtools")
  library("devtools")
#}

#devtools::install_github("TomKellyGenetics/graphsim")

#install.packages("gplots")

library("igraph")
library("graphsim")
library("gplots") 
library("scales")

  
graph_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"), c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph <- graph.edgelist(graph_edges, directed = TRUE)

#plot graph structure (Figure 1)
plot_directed(graph, border.node=alpha("black", 0.25), state=c(1, 1, 1, 1, 1, 1, 1, 1), col.arrow = c(alpha("navyblue", 0.5), alpha("navyblue", 0.5))[state], fill.node="lightblue")

#plot graph structure (Figure 1-inhibited)
plot_directed(graph, border.node=alpha("black", 0.25), state=c(1, 1, -1, 1, 1, 1, -1, 1), col.arrow = c(alpha("navyblue", 0.5), alpha("navyblue", 0.5))[state], fill.node="lightblue")
plot_directed(graph, state ="activating", layout = layout.kamada.kawai, cex.node=3, cex.arrow=5, arrow_clip = 0.2)

E(graph)$state <- c(1, 1, 1, 1, 1, 1, 1, 1)
E(graph)$state <- c(1, 1, -1, 1, 1, 1, -1, 1)
#plot_directed(graph)
  
#adjacency matrix
adj_mat <- make_adjmatrix_graph(graph)
adj_mat
#plot adjacency matrix
heatmap.2(make_adjmatrix_graph(graph), scale = "none", trace = "none",col = colorpanel(3, "grey75", "white", "blue"),colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
  
#relationship matrix
dist_mat <- make_distance_graph(graph, absolute = FALSE)
#plot relationship matrix
heatmap.2(make_distance_graph(graph, absolute = FALSE), scale = "none", trace = "none", col = bluered(50), colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
  
  
#sigma matrix directly from graph
sigma_mat <- make_sigma_mat_dist_graph(graph, 0.8, absolute = FALSE)
#show shortest paths of graph
shortest_paths <- shortest.paths(graph)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph, 0.8, absolute = FALSE), scale = "none", trace = "none", col = bluered(50), colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
  
#generate expression data directly from graph
#simulated expression data con n=100
  
expr <- generate_expression(100, graph, cor = 0.8, mean = 0, comm = FALSE, dist = TRUE, absolute = FALSE, state = state)
expr
write.table(expr, file = "simulated_expression_data_CDK5.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
  
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50), colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
  
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none", col = bluered(50), colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
  
####### PATH #######


#path
graph_edges_path <- rbind(c("CAST","CDK5R1"),c("Calpain1-Calpain2","CDK5R1"),c("Myr-CDK5R1","CDK5R1"),c("CDK5R1","CDK5-p25"),c("CDK5","CDK5-p25"),c("CDK5-p25","FOXO3"),c("FOXO3","FOXO3m"),c("FOXO3m","SOD2"),c("FOXO3m","BCL2L11"), c("FOXO3m","FASLG"))
graph_path <- graph.edgelist(graph_edges_path, directed = TRUE)

#plot path 
plot_directed(graph_path, border.node=alpha("black", 0.25), state=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1), col.arrow = c(alpha("navyblue", 0.5), alpha("navyblue", 0.5))[state], fill.node="lightblue")
plot_directed(graph_path, border.node=alpha("black", 0.25), state=c(1, 1, 1, -1, -1, 1, 1, 1, 1, 1), col.arrow = c(alpha("navyblue", 0.5), alpha("navyblue", 0.5))[state], fill.node="lightblue")


E(graph_path)$state <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
# inhibido-- path
E(graph_path)$state <- c(1, 1, 1, -1, -1, 1, 1, 1, 1, 1)

#adjacency matrix
adj_mat <- make_adjmatrix_graph(graph_path)
adj_mat
#plot adjacency matrix
heatmap.2(make_adjmatrix_graph(graph_path), scale = "none", trace = "none",col = colorpanel(3, "grey75", "white", "blue"),colsep = 1:length(V(graph_path)), rowsep = 1:length(V(graph_path)))

#relationship matrix
dist_mat <- make_distance_graph(graph_path, absolute = FALSE)
#plot relationship matrix
heatmap.2(make_distance_graph(graph_path, absolute = FALSE), scale = "none", trace = "none", col = bluered(50), colsep = 1:length(V(graph_path)), rowsep = 1:length(V(graph_path)))


#sigma matrix directly from graph
sigma_mat <- make_sigma_mat_dist_graph(graph_path, 0.8, absolute = FALSE)
#show shortest paths of graph
shortest_paths <- shortest.paths(graph_path)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph_path, 0.8, absolute = FALSE), scale = "none", trace = "none", col = bluered(50), colsep = 1:length(V(graph_path)), rowsep = 1:length(V(graph_path)))

#generate expression data directly from graph
#simulated expression data con n=100

expr <- generate_expression(100, graph_path, cor = 0.8, mean = 0, comm = FALSE, dist = TRUE, absolute = FALSE, state = state)
expr
write.table(expr, file = "simulated_expression_data_CDK5_inh.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50), colsep = 1:length(V(graph_path)), rowsep = 1:length(V(graph_path)))

#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none", col = bluered(50), colsep = 1:length(V(graph_path)), rowsep = 1:length(V(graph_path)))
