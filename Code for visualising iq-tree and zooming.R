library(phangorn)
library(ape)
library(phytools)
library(grDevices)
library(adegenet)
library(ggtree)
library(treeio)
library(dplyr)
library(tidytree)
library(phylobase)
library(cluster)
library(factoextra)
library(rbiom)

setwd("")

#Changing sample names

taxaid <- read.csv("(.csv file with desired new sequence names and old names)")

newnames <- data.frame(
  ERR = taxaid$ERR,
  Old = taxaid$Formula)

newnames

#Pairwise Matrix generation

alignedseq <- ape::read.FASTA("/Users/eyad/Desktop/LIDo LSHTM/LSHTM 1st Rotation/Data/aligned.consensus.masked.fasta")
alignedseq
finalseq <- updateLabel(alignedseq, old = newnames$Old, new = newnames$ERR)
class(finalseq)
seq.raw <- ape::dist.dna(alignedseq, model = "raw", pairwise.deletion = T)
pairwise <- round(1-(seq.raw)*100, digits = 3)

fviz_dist(dist.obj = pairwise, order = F, show_labels = T)

# Import IQ-Tree tree file and root it to NC_045512.2

myiqtree <- ape::read.tree("IQ_Tree file.tree")
myiqtree
ape::is.rooted(myiqtree)
myroot <- which(myiqtree$tip.label == "NC_045512.2")
treerooted <- phytools::reroot(myiqtree, myroot)
ape::is.rooted(treerooted)
treerooted

treerooted <- rename_taxa(treerooted, newnames, key = 2, value = 1) # Renaming Taxa

ape::write.tree(treerooted, file = "rooted.aligned.consensus.masked.fasta.treefile") # Save Tree

plotTree(treerooted,fsize = 0.5, lwd = 0.9, node.numbers = F, ftype = "off", pts = T,
         color = "darkgreen") 

# Adding HOCO Numbers to match Sequence Data for Isolation of individual clades

Hocofile <- read.csv("~/Desktop/LIDo LSHTM/LSHTM 1st Rotation/supplemental_data.csv")

hoco <- data.frame(
  ERR = Hocofile$ERR,
  Hoco = Hocofile$HOCONUMBER
)

generate_tree <- function(tree_obj, tips=NULL){
  plot1 <- ggtree(tree_obj, tips = tips) + #grep(readline(prompt = "Enter a household number: "), 
                                                      #treerooted$tip.label))) +
    geom_nodepoint(color="#0000FF", alpha=1/2, size=1) +
    geom_tippoint(color="#FF3333", alpha = 1, shape=20, size=3) +
    theme_tree2(plot.margin = margin(10, 10, 6, 6)) +
    geom_tiplab(align = F, linesize = .5, size =3, 
                hjust =0, offset = 0.0000000001) + #adjust hjust accordingly. - pushes it to higher x values (right), positive towards lower x values (left).
    geom_treescale(x = 0.00000008, y = 5.5 ,fontsize = 4, linesize = 0.5, offset = 0.2) +
    xlim(0, 7e-04) + #adjust xlim accordingly
    ylim(0, 6) + #adjust ylim accordingly. Adjust y axis in the geom_treescale geom.
    coord_cartesian(clip = "off") +
    geom_nodelab(aes(subset = as.numeric(label) > 10), size = 3) #+
    #labs(title = readline(prompt  = "Enter a title for the tree: ")) 

return(plot1)
}

#Plotting the Phylogenetic Tree using ggtree
ggtree(treerooted, branch.length = -10) + 
  scale_color_continuous(low= "darkgreen", high= "red")+
  geom_nodepoint(color="#b5e521", alpha=1/4, size=3) +
 geom_tippoint(color="#FDAC4F", alpha = 1/3, shape=8, size=2) +
  theme_tree2(plot.margin = margin(10, 10, 6, 6)) +
  geom_tiplab(align = T, linesize = .5, size = 0.7, 
              hjust =-0.1, offset = 0.0001, nudge_y = 0, parse = F) +
  geom_treescale(x = 0.0008, y = 130 ,fontsize = 4, linesize = 0.5, offset = 2) +
  xlim(0, 0.0020) +
  coord_cartesian(clip = "off") +
  geom_nodelab(aes(subset = as.numeric(label) > 95), size = 3)



# Plotting individual households

ggtree(subtree(treerooted, tips = grep(readline(prompt = "Enter a household number: "), treerooted$tip.label))) +
  geom_nodepoint(color="#0000FF", alpha=1/2, size=1) +
  geom_tippoint(color="#FF3333", alpha = 1, shape=20, size=3) +
  theme_tree2(plot.margin = margin(10, 10, 6, 6)) +
  geom_tiplab(align = F, linesize = .5, size =3, 
              hjust =0, offset = 0.0000000001) + #adjust hjust accordingly. - pushes it to higher x values (right), positive towards lower x values (left).
  geom_treescale(x = 0.00000008, y = 5.5 ,fontsize = 4, linesize = 0.5, offset = 0.2) +
  xlim(0, 7e-04) + #adjust xlim accordingly
  ylim(0, 6) + #adjust ylim accordingly. Adjust y axis in the geom_treescale geom.
  coord_cartesian(clip = "off") +
  geom_nodelab(aes(subset = as.numeric(label) > 10), size = 3) +
  labs(title = readline(prompt  = "Enter a title for the tree: ")) 

# Plotting a subset of the phylogeny tree with color coding

## subset the tree

p <- tree_subset(treerooted, node = grep('ERR6114168', treerooted$tip.label), # grep a sequence from household of interest
                 levels_back = 10) # go n nodes back till all household members are in the subset tree

## create a list with the tip labels in the subset tree
## order must be the same as the tip order in the subset tree

view(p) # to get the tip order

branches <- list(ERR6115673 = 20, ERR6114753 = 19, ERR6117249 = 21,
                 ERR6117213 = 22, ERR5989911 = 18, ERR5990034 = 23,
                 ERR6061736 = 24, ERR6115144 = 15, ERR6114168 = 14,
                 ERR6113712 = 13, ERR5971409 = 12, ERR6116969 = 17,
                 ERR6111499 = 16,  ERR6060085 = 25, ERR6115600 = 29,
                 ERR6113216 = 28, ERR6063288 = 27, ERR6062614 = 26,
                 ERR6063683 = 30, ERR5992079 = 8, ERR5990363 = 7, 
                 ERR5991896 = 6, ERR5989497 = 5, ERR5991540 = 4,
                 ERR5991500 = 3, ERR5989091 = 2, ERR5971008 = 1,
                 ERR5989848 = 10, ERR5971582 = 9, ERR5990182 = 11,
                 ERR5971152 = 31)

## group the tree branches with the tree

tree <- groupOTU(p, branches)
tree2 <- groupOTU(p, col$HOCO)

## plot the tree with specific color corresponding to all members of same household

ggtree(tree)+
  geom_tiplab(aes(color = group), show.legend = F, #aes color = group colors the labels as grouped in groupOTU
              align = F, linesize = .5, size = 3, 
              hjust =0.7, offset = 0, nudge_y = 0,
              nudge_x = .0000089) +
  scale_color_manual(values = c(ERR5971409 = 'red', ERR5989911 = 'red', ERR5971152 = 'red', ERR5990034 = 'red',
                                ERR5971008 = 'blue',  ERR5971582='blue',ERR5989848='blue', ERR5990182='blue',
                                ERR5989497 = 'green', ERR5990363 = 'green', ERR5991896 = 'green',ERR5992079 = 'green',
                                ERR6060085 = 'violet', ERR6061736 = 'violet',ERR6063683 = 'violet', ERR6117213='violet',ERR6117249 = 'violet',
                                ERR6113712 = 'black', ERR6062614 = 'purple', ERR6063288 = 'purple', ERR6113216 = 'purple', ERR6115600 = 'purple',
                                ERR6114168 = 'brown',ERR6114753='brown', ERR6115144 = 'brown', ERR6115673='brown',
                                ERR5989091 = 'orange', ERR5991500 = 'orange', ERR5991540 = 'orange',
                                ERR6111499 = 'lightsalmon', ERR6116969 = 'lightsalmon')) +
  guides(color = guide_legend(title = "Households")) # removing the legend

####
# End of script