PROJECT SUMMARY

Most toxicological studies investigate single compounds, overlooking the combined effects of multiple exposures that more accurately reflect real-world conditions. Many such compounds can function as endocrine-disrupting chemicals (EDCs), disturbing hormonal homeostasis. Because the chemical exposome—the total chemical burden accumulated by an individual—differs widely among people, its biological significance remains uncertain. To explore these effects, using a high-throughput screening (HTS) 384-well plate version of the OECD No. 455 assay, we analyzed 16 human exposomes and their 24 persistent organic pollutant (POP) components for their influence on estrogen receptor (ER) activity, based on blood profiles from participants in the Swedish Västerbotten Intervention Programme. 
The data were analyzed using a linear mixed-effects model, and Dunnett’s multiple comparison was applied as a post-hoc test.

The code was used in this paper: High-throughput screening of estrogen receptor activity for personalized mixtures of persistent organic pollutants detected in blood of Swedish adults. 

PROJECT STRUCTURE

The dataset is composed of human VM7Luc4E2 cells exposed to 24 POPs compounds and 16 human exposomes. The molar concentration of the compounds ranged from 1 nM to 10 µM. Experiments were performed in 3 batches. The dataset was curated in an Excel file. 

Code/analysis.R contains the code for statistical analysis (linear mixed-effects model, ANOVA, and post-hoc) and will generate:

(1) plots of all compounds (img/plots.pdf)

(2) the statistical output (data/fit_linear_mixed_model_output.xlsx). 

All analyses were performed in R.

