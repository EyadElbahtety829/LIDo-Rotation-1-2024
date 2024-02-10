setwd("~/Desktop")

tiff("transphylo cartoon.tiff", width = 1080, height = 1080, units = "px")
layout(matrix(c(1,2,3,4), 2, 1, byrow = T))
par(adj = 0, cex.main = 2)

set.seed(0)
sim = simulateOutbreak(nSampled = 3, dateStartOutbreak = 2019, dateT = 2023)
plot(sim) + 
  title(main = "I)")


sim = simulateOutbreak(nSampled = 3, dateStartOutbreak = 2019, dateT = 2023)
plot(sim) + 
  title(main = "II)")

dev.off()


