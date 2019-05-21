library(TreeSim)

rootDir = "/data/Suha/birth-death" #path top the root directory
m = c(20,40,60,80,100) #number of taxa

for (i in m) {
	f = paste("simTrees",i,"taxa.txt", sep = "")
	for (j in 1:3960){
	  Lambda = runif(1, min = 0, max = 1) #birth rate
	  Mu = runif(1, min = 0, max = Lambda) #death rate
	  fraction = runif(1, min = 0, max = 1) #fraction of sampled taxa
	  lis = sim.bd.taxa(i, 1, Lambda, Mu, fraction,  complete = FALSE, stochsampling = FALSE)
	  write.tree(lis[[1]], file.path(rootDir, f),  append = TRUE)
	}
}
