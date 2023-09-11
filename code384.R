library(ggplot2)
library(openxlsx)
library(cowplot)
library(tidyverse)
library(ggpubr)

setwd("c:/Users/LuaReis/Documents/Postdoc/Stockholm university/Exposomics/Estrogen/Results/Luciferase/")
file <- "rstudioLPcompound08052023.xlsx"

plate <- read.xlsx(file, 1)
conc <- read.xlsx(file, 2) %>% log10()

#out e um objeto dentro da funcao, logo e possivel colocar objetos dentro de funcoes
#apply e uma funcao, que aplica uma funcao em uma matrix ou data frame, 
#especificamente nas linhas (1) ou colunas (2)

#criou-se uma funcao com 2 parametros, 1) a tabela excel com os dados e 2) o subgrupo
#de dados dentro dessa tabela, chamado de range1-6
extractvalues <- function(plate, range1, range2, range3, range4, range5, range6) {
  out <- rbind(as.matrix(plate[range1, range2]), as.matrix(plate[range3, range4]),
  as.matrix(plate[range5, range6]))
  out <- apply(out, c(1,2), as.numeric)
  rownames(out) <- NULL
  colnames(out) <- NULL
  out
}

e2 <- extractvalues(plate, 16:17,4:14, 34:35, 4:14, 52:53, 4:14)
dmso <- extractvalues(plate, 16:17,16:18, 34:35, 16:18, 52:53, 16:18) %>% mean()

#funcao com a media e desvio padrao de cada concentracao dos compostos
descre.stats <- function(x) {
  mean.x <- apply(x, 2, mean) %>% as.matrix %>% t
  sd.x <- apply(x, 2, sd) %>% as.matrix %>% t
  meansd.x <- rbind(mean.x,sd.x)
  colnames(meansd.x) <- conc[,9]
  return(as.matrix(meansd.x))
}

#criar objeto inducao e funcao if(any) para processar os dados somente se eles forem >4
induction <- descre.stats(e2)
induction[1, ] <- induction[1, ] / dmso

if(any(induction[1,]>4)) {

controle2 <- extractvalues(plate, 16:17,4:14, 34:35, 4:14, 52:53, 4:14) %>% - dmso  
compound1 <- extractvalues(plate, 3:5,4:14, 21:23, 4:14, 39:41, 4:14)  %>% - dmso 
compound2 <- extractvalues(plate, 6:8,4:14, 24:26, 4:14, 42:44, 4:14) %>% - dmso
compound3 <- extractvalues(plate, 9:11,4:14, 27:30, 4:14, 45:48, 4:14) %>% - dmso
compound4 <- extractvalues(plate, 12:14,4:14, 31:33, 4:14, 49:51, 4:14) %>% - dmso
compound5 <- extractvalues(plate, 3:5,15:25, 21:23, 15:25, 39:41, 15:25) %>% - dmso
compound6 <- extractvalues(plate, 6:8,15:25, 24:26, 15:25, 42:44, 15:25) %>% - dmso
compound7 <- extractvalues(plate, 9:11,15:25, 27:30, 15:25, 45:48, 15:25) %>% - dmso
compound8 <- extractvalues(plate, 12:14,15:25, 31:33, 15:25, 49:51, 15:25) %>% - dmso

#normalizando min/max relativo ao E2
maxe2 <- colMeans(controle2) %>% max()
mine2 <- colMeans(controle2) %>% min()

controle2 <- ((controle2-mine2)/(maxe2-mine2))*100
compound1 <- ((compound1-mine2)/(maxe2-mine2))*100
compound2 <- ((compound2-mine2)/(maxe2-mine2))*100
compound3 <- ((compound3-mine2)/(maxe2-mine2))*100
compound4 <- ((compound4-mine2)/(maxe2-mine2))*100
compound5 <- ((compound5-mine2)/(maxe2-mine2))*100
compound6 <- ((compound6-mine2)/(maxe2-mine2))*100
compound7 <- ((compound7-mine2)/(maxe2-mine2))*100
compound8 <- ((compound8-mine2)/(maxe2-mine2))*100


# RLU / area scaled: min = avg DMSO, max = avg E2
values <- ( values - mean(values[8, 7:12]) ) / ( mean(values[8, 1:6]) - mean(values[8, 7:12]) ) * 100
  
# generate dataset for plot, assumes 2 compounds
compounds <- colnames(conc)
  
data <- data.frame(conc = unlist(conc))
#creates a data frame named data. It seems like conc is a list or vector containing concentration values. 
#The unlist() function is used to convert this list into a single vector
#and the vector is assigned as a column named "conc" in the data data frame.

data$compounds <- rep(colnames(conc), each = 7)
#a new column named "compounds" is added to the data data frame.
#This column is populated with the names of compounds, which are assumed to be the column names of the conc data. 
#The colnames(conc) returns the names of the columns, and the rep() function repeats these column names 7 times each to match the length of the data.

rownames(data) <- NULL
#This line removes the row names from the data data frame. 
#Row names are typically used to uniquely label rows in a data frame, but in this case, the row names are being reset to NULL.
#This means that rows will be indexed using default integer values.

data$values <- c(apply(values[1:7, 1:6], 1, mean), apply(values[1:7, 7:12], 1 , mean))
#A new column named "values" is added to the data data frame. 
#This column is populated with mean values calculated row-wise from the values data.
#It seems that values is expected to be a matrix or data frame containing some values. 
#The apply() function is used to apply the mean() function row-wise
#calculating the mean for the first 6 columns (1:6) and the next 6 columns (7:12) separately for the first 7 rows

data$sd <- c(apply(values[1:7, 1:6], 1, sd), apply(values[1:7, 7:12], 1 , sd))
#Similarly, a new column named "sd" is added to the data data frame.
#This column is populated with standard deviation values calculated row-wise from the values data. 
#The apply() function with the sd() function is used to calculate standard deviations for the same ranges of columns as before.

data$sd <- ifelse(data$sd < 4, NA, data$sd)
#This line applies a condition to the sd column of the data data frame. If the standard deviation value in a row is less than 4, it is replaced with NA (indicating missing data). 
#This could be done to filter out low-standard-deviation values.

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
  
}


