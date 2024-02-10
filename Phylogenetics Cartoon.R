setwd("~/Desktop")

library(ape)

set.seed(1)
tree <- rtopology(n = 3, tip.label = c("B", "A", "C"), equiprob = T, rooted = T)
par(mar = c(.5,.5,.5,.5))

# setting format and layout parameters of plot

tiff("transmission cartoon.tiff", width = 1080, height = 1080, units = "px")
layout(matrix(c(1,2,3,4), 2, 2, byrow = T))
par(adj = 0, cex.main = 2)

# A infecting B and C
plot(tree, node.color = c("Blue","Green", "Red", "Red", "Red", "Blue"), 
     use.edge.length = F, edge.width = 3, label.offset = .001, cex = 2,
     root.edge = T
     ) 
  title(main = "I)")
  

# A infecting B, B infecting C

plot(tree, node.color = c("Red" , "Green", "Green", "Red", "Green"),
     edge.color = c("Red", "Blue", "Green", "Red"), 
     use.edge.length = F, edge.width = 3, root.edge = T, label.offset = .001,
cex = 2) 
  title(main = "II)")
  

# A or B infecting C
plot(tree, node.color = c("Red" , "Red", "Blue", "Red", "Red"),
     use.edge.length = F, edge.width = 3, root.edge = T, label.offset = .001,
cex = 2) 
  title(main = "III)")

# A shares common ancestor with B and C

plot(tree, node.color = c("Red" , "Green", "Blue", "Orange", "Orange"),  
     use.edge.length = F, edge.width = 3, root.edge = T, label.offset = .001,
cex = 2) 
  title(main = "IV)")
  
  
dev.off()
 
