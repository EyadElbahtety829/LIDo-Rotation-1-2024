library(ape)
library(TransPhylo)
library(ggtree)
library(tidyverse)
library(lubridate)
library(coda)
library(lattice)
library(viridisLite)
library(rasterVis)
library(phytools)
###
#Setting working directory
setwd("/Users/eyad/Desktop/LIDo LSHTM/LSHTM 1st Rotation/Data/TreeTime files/n=100 Alpha/2023-12-17_treetime")

#Reading tree file nexus and newick
nex <- ape::read.nexus("timetree.nexus.tree") 

#Fixing branch length zero error

t <- multi2di(nex) # remove multifurcations to avoid negative-length branches error.
t$edge.length <- pmax(t$edge.length,1/365) # use a day as minimum branch length

final_tree <- drop.tip.phylo(t, "NC_045512.2")
#Plotting tree for visualization
ggtree(final_tree, mrsd = "2021/06/11") + 
  theme_tree2() +
  geom_rootedge()

#creating phylo tree from phylo object:v23-4 

covptree <- ptreeFromPhylo(tr = final_tree, dateLastSample = 2021.44)
plot(covptree)

#inferring transmission tree from ptree:

set.seed(0) 
mcmc_tree <- inferTTree(ptree = covptree, mcmcIterations = 1e5,
                        thinning = 1,
                        #w.mean = 0.013, w.std = 0.0041, #for delta variant
                        #w.mean = 0.015, w.std = 0.0048, #for alpha variant
                       startOff.r = 1.22, updateOff.p = T
                       #1.38 for delta, 1.22 for alpha
                        ) 
plot(mcmc_tree) +
  print(mcmc_tree)

#confirming mcmc analysis results

mcmc = convertToCoda(mcmc_tree)
effectiveSize(mcmc)

#inference of results - medoid of transmission trees

med = medTTree(mcmc_tree, burnin = .5)

#plotting coloring of dated phylogeny

plotCTree(med, cex = .2, showStars = T) + 
  title(main = "Colored Dated Phylogeny for Alpha Variants")


#plotting transmission tree

plotTTree2(extractTTree(med), showLabels = T, showMissingLinks = 2,
                                 cex = .1) +#type = "detailed",
    # w.shape = 10, w.scale =.0013, cex = 0.5) +
  title(main = "Transmission Tree for Alpha Variants") 


 #plotting probability matrix of who infected whom

matWIW = computeMatWIW(mcmc_tree, burnin = .5)
col <- viridis::magma(100)
levelplot(matWIW, 
          xlab = "Source Case", ylab = "Recipient Case", main = "Probability matrix of transmission pairs - Alpha", 
          col.regions = col, 
          scales = list(y = list(cex = .2), x = list(rot = 90, cex = .2), relation = "same",tck = c(1,0)),
)

#building intermediairy matrix

matDist = computeMatTDist(mcmc_tree)
levelplot(matDist, xlab = "Source", ylab = "Recipient", main = "Distribution matrix for transmission pairs - Alpha",
          col.regions = col,
          scales = list(y = list(cex = .1), x = list(rot = 90, cex = .1), relation = "same",
                        tck = c(1,0)))

#inferring transmission properties for each individual

tim = getInfectionTimeDist(mcmc_tree, k = readline(prompt = "Please enter a case ID: "), 
      show.plot = T) 
off = getOffspringDist(mcmc_tree,  k = readline(prompt = "Please enter a case ID: "), 
      show.plot = T) 

#MCMC realized distribution, generation times, and sampling times


getIncidentCases(mcmc_tree,show.plot = T, numBins = 30) +
  title(main = "Incident Cases - Delta")
getGenerationTimeDist(mcmc_alpha,show.plot = T, numBins = 30) +
  title(main = "Generation Time Distribution - Delta")
getSamplingTimeDist(mcmc_tree,show.plot = T, numBins = 30) 
