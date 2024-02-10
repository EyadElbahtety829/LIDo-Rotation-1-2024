library(tidyverse)
library(ggplot2)
library(egg)

setwd("/Users/eyad/Desktop/LIDo LSHTM/LSHTM 1st Rotation/Data")

df <- readxl::read_excel("Serial Intervals Alpha & Delta.xlsx")

#create separate df for each variant

delta <- transmute(df,
                   Delta_Gene = df$Delta_Gene,
                   Delta_Epi = df$Delta_Epi)

alpha <- transmute(df,
                   Alpha_Gene = df$Alpha_Gene,
                   Alpha_Epi = df$Alpha_Epi)
delta <- drop_na(delta)
alpha <- drop_na(alpha)

#plotting both groups of each variant spearately
## Genomic distribution## Epi distribution

tiff("Serial Intervals Hist.tiff", width = 1080, height = 1080, units = "px")
layout(matrix(c(1:4), 2, 2, 
              byrow = T))

# alpha_gene
alpha_gene <-  ggplot(alpha, x = `Alpha_Gene`) +
  geom_histogram(aes(x = `Alpha_Gene`, y = ..density..),
                 bins =20, colour = 1, fill = 'gray') +
  geom_density(aes(x = `Alpha_Gene`), color = 'red') +
  xlab("Days of recipient symptom onset \n start date from source case") +
  ylab("Density") +
  theme_classic() +
  scale_x_continuous(breaks = c(-9:12)) +
  geom_vline(aes(xintercept = mean(`Alpha_Gene`)), linetype = "dashed", col = "blue",
             lwd = 1) +
  ggtitle("I)")

# alpha_epi
alpha_epi <- ggplot(alpha, x = `Alpha_Epi`) +
  geom_histogram(aes(x = `Alpha_Epi`, y = ..density..),
                 bins =20, colour = 1, fill = 'gray') +
  geom_density(aes(x = `Alpha_Epi`), color = 'red') +
  xlab("Days of recipient symptom onset \n start date from source case") +
  ylab("Density") +
  theme_classic() +
  scale_x_continuous(breaks = c(-9:12)) +
  geom_vline(aes(xintercept = mean(`Alpha_Epi`)), linetype = "dashed", col = "blue",
             lwd = 1) +
  ggtitle("II)")

# delta_gene
delta_gene <- ggplot(delta, x = `Delta_Gene`) +
  geom_histogram(aes(x = `Delta_Gene`, y = ..density..),
                 bins =20, colour = 1, fill = 'gray') +
  geom_density(aes(x = `Delta_Gene`), color = 'red') +
  xlab("Days of recipient symptom onset \n start date from source case") +
  ylab("Density") +
  theme_classic() +
  scale_x_continuous(breaks = c(-9:12)) +
  geom_vline(aes(xintercept = mean(`Delta_Gene`)), linetype = "dashed", col = "blue",
             lwd = 1) +
  ggtitle("III)")

# delta_epi
delta_epi <- ggplot(delta, x = `Delta_Epi`) +
  geom_histogram(aes(x = `Delta_Epi`, y = ..density..),
                 bins =20, colour = 1, fill = 'gray') +
  geom_density(aes(x = `Delta_Epi`), color = 'red') +
  xlab("Days of recipient symptom onset \n start date from source case") +
  ylab("Density") +
  theme_classic() +
  scale_x_continuous(breaks = c(-2:12)) +
  geom_vline(aes(xintercept = mean(`Delta_Epi`)), linetype = "dashed", col = "blue",
             lwd = 1) +
  ggtitle("IV)")
           
ggarrange(alpha_gene, alpha_epi, delta_gene, delta_epi)

dev.off()

