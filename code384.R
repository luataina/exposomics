library(ggplot2)
library(openxlsx)
library(cowplot)
library(tidyverse)
library(ggpubr)

# set working directory to the folder containing the excel files
# setwd() # uncomment this

setwd("c:/Users/LuaReis/Desktop/")
file <- "rstudioLPcompound08052023.xlsx"

plate <- read.xlsx(file, 1)
conc <- read.xlsx(file, 2)
conc <- log10(conc)

#out e um objeto dentro da funcao, logo e possivel colocar objetos dentro de funcoes
extractvalues <- function(plate, range1, range2, range3, range4, range5, range6) {
  out <- rbind(as.matrix(plate[range1, range2]), as.matrix(plate[range3, range4]), as.matrix(plate[range5, range6]))
  out <- apply(out, c(1,2), as.numeric)
  
  rownames(out) <- NULL
  colnames(out) <- NULL
  out
}


compound1 <- extractvalues(plate, 3:5,4:14, 21:23, 4:14, 39:41, 4:14)
#apply e uma funcao, que aplica uma funcao em uma matrix ou data frame, 
#especificamente nas linhas (1) ou colunas (2)

mean_compound1 <- apply(compound1, 2, mean) %>% as.matrix %>% t
sd_compound1 <- apply(compound1, 2, sd) %>% as.matrix %>% t
statsdescpt_compound1 <- rbind(mean_compound1,sd_compound1)
colnames(statsdescpt_compound1) <- conc[,1]


compound2
compound3
compound4
compound5
compound6
compound7
compound8



# RLU / area normalised
values <- luciferace / covered_area

# RLU / area - avg DMSO
values <- values - mean(values[8,7:12])

# RLU / area scaled: min = avg DMSO, max = avg E2
values <- ( values - mean(values[8, 7:12]) ) / ( mean(values[8, 1:6]) - mean(values[8, 7:12]) ) * 100

# generate dataset for plot, assumes 2 compounds
compounds <- colnames(conc)

data <- data.frame(conc = unlist(conc))
data$compounds <- rep(colnames(conc), each = 7)
rownames(data) <- NULL
data$values <- c(apply(values[1:7, 1:6], 1, mean), apply(values[1:7, 7:12], 1 , mean))
data$sd <- c(apply(values[1:7, 1:6], 1, sd), apply(values[1:7, 7:12], 1 , sd))
data$sd <- ifelse(data$sd < 4, NA, data$sd)
#data$shape <- rep(c(21, 22), each = 7)

p <- ggplot(data, aes(factor(conc), values, shape = compounds, fill = compounds, group = compounds)) +
  geom_line(aes(color = compounds)) + 
  geom_errorbar(aes(ymin = values - sd, ymax = values + sd), width = 0.2) +
  geom_hline(yintercept = 80, alpha = 0.5, lty = 2) +
  geom_point(size = 4, color = "black") +
  #geom_point(size = 3) +
  theme_pubclean(base_size = 16) +
  labs(x = "log concentrations (M)", y = "RLU (% E2)", title = "RLU / Area", fill = "", shape = "", color = "") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  scale_shape_manual(values = c(21, 22)) 
  

# save pdf
pdfname <- paste0(sub(".xlsx", "", file), ".pdf")
pdf(pdfname, width = 10)
p
dev.off()
