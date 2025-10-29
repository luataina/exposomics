# --- pipeline --------------------------------------------------------------- #

# import curated excel file, expect 3 sheets: agonist, antagonist, viability
# files to use: misc/mixtures__rmix.xlsx & misc/single_cmpd_rmix.xlsx

# log transform the RLT column
# run mixed linear model, run post-hoc test, save output to excel 
# and exponent the measures (reverse from log transformation)

# save anova, linear model and post hoc summaries into an excel
# plot each compound

# --- load environment ------------------------------------------------------- #

library(openxlsx)
library(writexl)
library(tidyverse)
library(ggpubr)   # signficance bars
library(scales)

# statisical packages
library(car)      # Anova
library(lme4)     # lm model
library(multcomp) # post-hoc test
library(broom)    # post-hoc test

importer <- function(xlsxfile){
  xlsxfile <- file.path("misc", xlsxfile)
  
  # expect 3 sheets: agonist, antagonist, viability
  
  sheets <- c("agonist", "antagonist", "viability")
  dat <- lapply(sheets, function(sheet) read.xlsx(xlsxfile, sheet = sheet))
  names(dat) <- sheets

  return(dat)
}

fix_columns <- function(dat){
  cols <- names(dat) %in% c("cmpd", "treat", "n")
  dat[, !cols] <- sapply(dat[, !cols], function(x) log2(as.numeric(x)))

  # remove dashes in compounds and white space
  dat[, "cmpd"] <- sapply(dat[, "cmpd"], function(x){
    x <- gsub("-", "", x)
    x <- gsub(" ", "", x)
    return(x)
  })

  return(dat)
}

run_glm_analysis <- function(dat){
  # compounds to test, if present in all 3 sheets
  compounds <- lapply(dat, function(x){
    x <- unique(x$cmpd)
    sort(x)
  }) %>% Reduce(intersect, .)

  analysis <- list()
  for (compound in compounds){

    print(paste("Testing compound:", compound))

    lm <- lapply(dat, function(x){
      out <- subset(x, cmpd %in% compound)
      out$treat <- factor(out$treat)
      out$n <- factor(out$n)
      
      lmer(RLT ~ treat + n + (1 | treat : n), data = out, REML = TRUE)
    })

    # save post-hoc output
    posthoc <- list()
    for (parameter in names(lm)) {
      posthoc[[parameter]] <- data.frame(tidy(glht(lm[[parameter]]))[, -2], parameter = parameter)
    }
    posthoc <- Reduce(function(x,y) rbind(x,y), posthoc)
    posthoc$compound <- compound

    # save anova output
    anova <- list()
    for (parameter in names(lm)) {
      anova[[parameter]] <- calc_anova(lm[[parameter]], parameter)
    }
    anova <- Reduce(function(x,y) rbind(x,y), anova)
    anova$compound <- compound

    # save fold change + CI output
    l1 <- list()
    for (parameter in names(lm)) {
      l1[[parameter]] <- calc_CI(lm[[parameter]], parameter)
    }
    data <- Reduce(function(x,y) rbind(x,y), l1)
    rownames(data) <- NULL
    data$compound <- compound
  
    out <- list()
  
    out$data <- data
    out$anova <- anova
    out$posthoc <- posthoc

    analysis[[compound]] <- out
  }
  return(analysis)
}

# used inside run_glm_analysis()
calc_anova <- function(res, parameter){
  aov <- Anova(res, type = "III")
  out <- data.frame(variable = rownames(aov), aov, parameter = parameter)
  rownames(out) <- NULL
  
  out
}

# used inside run_glm_analysis()
calc_CI <- function(res, parameter){
  res <- summary(res)
  rows <- grepl("treat", rownames(res$coefficients))
  treatments <- sub("treat", "", rownames(res$coefficients)[rows])
  
  fold_change <- exp(res$coefficients[rows, 1])
  pval <- 2*(1 - pnorm(abs(res$coefficients[rows, 3])))
  upper_ci <- exp(res$coefficients[rows, 1] + 1.96*res$coefficients[rows, 2])
  lower_ci <- exp(res$coefficients[rows, 1] - 1.96*res$coefficients[rows, 2])
  
  data <- data.frame(treatment = as.numeric(treatments),
                     fold_change,
                     p_value = pval,
                     upper_ci,
                     lower_ci,
                     parameter = parameter)
  return(data)
}

plotter <- function(data, compound, pltsize = 9, pointsize = 2.5, sign_ind = 3){
  # table for plotting significance indicators: *, **, ***, ns = nothing
  statdata <- data.frame(
    parameter = data$parameter,
    y.position = data$upper_ci * 1.025, # needed to see sign-indicator
    group1 = data$treatment,
    group2 = data$treatment,
    p_signif = sapply(data$p_value, function(x){
      if (x < 0.001) return("***")
      if (x < 0.01) return("**")
      if (x <= 0.05) return("*")
      if (x > 0.05) return(NA)
  }))

  ggplot(data, aes(treatment, fold_change)) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5, 
               color = "darkorange3") +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    geom_line(aes(group = 1), alpha = 0.7, linewidth = 0.5) +
    geom_point(aes(fill = parameter), shape = 21, size = pointsize) +
    facet_wrap(. ~ parameter, scales = "free") + 

    scale_x_log10() +
    scale_y_continuous(labels = comma, limits = c(0, 30)) +

    theme_pubr(pltsize) + 
    theme(
      legend.position = "none", 
      plot.title = element_text(hjust = 0.5),
      strip.text.x = element_text(size = 10)) +
    labs(x = "Blood concentration (fold-change)", y = "Fold of control", title = compound) +
    stat_pvalue_manual(statdata, label = "p_signif", size = sign_ind, 
                       tip.length = 0)
  }

if (!dir.exists("data")) dir.create("data")
if (!dir.exists("img")) dir.create("img")

# --- script ----------------------------------------------------------------- #

# files to use misc/mixtures__rmix.xlsx & misc/single_cmpd_rmix.xlsx
luas_files <- c("mixtures__rmix.xlsx", "single_cmpd_rmix.xlsx")
files <- list.files("misc")
files <- files[files %in% luas_files]

# iterate over files found in misc/
for (file in files){
  dat <- importer(file)

  # make measures numeric and log transform
  dat <- lapply(dat, fix_columns)

  # run linear model analysis
  dat <- run_glm_analysis(dat)

  # Note, error message: "boundary (singular) fit: see help('isSingular')"
  # means = model was fit, but random effects were small
  # source: https://stackoverflow.com/questions/60028673/lme4-error-boundary-singular-fit-see-issingular

  # save statistical output
  anova_summary <- lapply(dat, function(x) x$anova) %>%
    Reduce(function(x,y) rbind(x,y), .)
  linearmodel_summary <- lapply(dat, function(x) x$data) %>%
    Reduce(function(x,y) rbind(x,y), .)
  posthoc_summary <- lapply(dat, function(x) x$posthoc) %>%
    Reduce(function(x,y) rbind(x,y), .)

  xlsxfile <- list(
    "linearmodel_summary" = linearmodel_summary,
    "anova_summary" = anova_summary,
    "posthoc_summary" = posthoc_summary)

  filename <- sub(".xlsx", "", file)
  filename <- paste0(filename, "_", "GLM_output.xlsx")
  write_xlsx(xlsxfile, path = file.path("data", filename))
  
  # Make plots
  tmp <- lapply(dat, function(x) x$data) # dirty fix

  # base p-values on compound-level
  d <- anova_summary[anova_summary$variable == "treat" & anova_summary$Pr..Chisq. <= 0.05, ]

  cont <- unique(posthoc_summary$contrast)
  cont <- cont[grepl("treat", cont)]
  posthoc_sub <- subset(posthoc_summary, contrast %in% cont)

  # keep track of which are significant
  data <- lapply(tmp, function(x) cbind(x, sign = FALSE))
  for (i in seq_len(nrow(d))){
    comp <- d$compound[i]
    param <- d$parameter[i]
    
    out <- data[[comp]]
    
    contrast <- subset(posthoc_sub, parameter %in% param & compound %in% comp & adj.p.value <= 0.05)$contrast
    contrast <- sub("treat", "", contrast)
    
    rows <- out$treatment %in% contrast & out$parameter %in% param & out$compound %in% comp
    
    out$sign[rows] <- TRUE
    
    data[[comp]] <- out
  }

  # Set p-value of non-significant post hoc correlations = 1, quick-fix
  for (i in seq_len(length(data))){
    sign <- data[[i]]$sign
    data[[i]]$p_value[!sign] <- 1
  }

  compounds <- names(data)
  
  filename <- sub(".xlsx", "", file)
  filename <- paste0(filename, "_plots.pdf")
  
  pdf(file.path("img", filename), height = 3)
  
  for (compound in compounds){
    p <- plotter(data[[compound]], compound)
    print(p)
  }

  dev.off()
  
  cat(paste0(
    "\n===============================================================\n\n",
    "Summary tables & plots generated for : ", file,
    "\n\n===============================================================\n"
  ))
}
