# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# date 220919
# Author:   Andreas Halgreen Eiset, eiset@ph.au.dk
# Title:    Data set for propensity score model
# Licence:  GNU GPLv3
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


# NO NEED TO RUN: make data set with all predictors and exposure as dummy variable---------------

source("1a_loadTidy_otutax_and_matrices.R")

#detach all packages loaded via the source
lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach,
       character.only = TRUE, unload = TRUE)
#remove everything but the stated from the environment
rm(list = setdiff(ls(), "otutax_trimd3"))

tmps <- new.env()
tmps$pckg <- c("tidyverse", "rms")
lapply(tmps$pckg, library, character.only = TRUE)
sapply(tmps$pckg, packageVersion)

dta <- otutax_trimd3$data$otuPresAbsc %>%
  pivot_longer(-c(1:3),
               names_to = "smplID",
               values_to = "prsnt")

###se analysis_plan for choice of bacteria to single out

# Choose Christensenellaceae/Akkermansia/Anginosus/Dolichum/...
bact <- c("christensen", "akker", "streptococ","lactobacillac", "clostridiac", 
          "dehalobact", "proteobact", "lachnospir", "fusobact", "gallolyticus", 
          "smithii")

###what taxon_id are we talking about?
taxIDs <- unique(dta$taxon_id[grep("christensen",
                                   dta$taxmy,
                                   ignore.case = TRUE)])

###Which is family and which is genus?
dta %>%
  filter(grepl(paste(taxIDs, collapse = "|"), taxon_id)) %>%
  select(taxon_id, taxmy) %>%
  distinct(taxon_id, .keep_all = TRUE) %>%
  mutate(taxmy = gsub("^k__.*Clostridiales", #i.e. begins with "k__" and ends with "Clostridiales;"
                      "", .data$taxmy))#cut everything in front of family to allow for print

###in base R:
#unique(dta[grepl(paste(taxIDs, collapse = "|"), dta$taxon_id),
#          "taxmy"])
###"agc" is the family identifyer

##Create dummy for one of the bacteria being present or absent
###function to create the dummy
is_it_there <- function(bact.list) {
  a <- bact.list
  b <- dta %>%
    select(smplID, prsnt, taxmy) %>%
    mutate(c = if_else(
      str_detect(taxmy,
                 regex(a, ignore_case = TRUE)) &
        prsnt == "TRUE", 1, 0)) %>%
    rename({{ a }} := "c") %>%  #see https://dplyr.tidyverse.org/articles/programming.html#different-input-and-output-variable
    select(smplID, {{ a }})
  return(b)
}

#create df with one row for each id and dummies indicating whether or not bacteria is present
dta_dist_id <- map_dfc(bact, is_it_there) %>%
  select(smplID = smplID...1, {{ bact }}) %>%
  pivot_longer(cols = -smplID,
               names_to = "taxa",
               values_to = "dummy") %>%
  filter(dummy == 1) %>%
  group_by_all() %>%
  distinct(smplID) %>%  #only want one dummy per bacteria (for our analysis it dosen't matter if the same bacteria has been found several times (i.e. has several "1" creating more rows for same ID))
  ungroup() %>% #now we have only unique rows that is one row with id for each pattern of dummies
  pivot_wider(id_cols = smplID,
              names_from = taxa,
              values_from = dummy,
              names_prefix = "dum_") %>%
  replace(is.na(.), 0) %>%
  mutate_if(is.numeric, as.factor)

#colnames(dta)[-1] <- paste("dum_", colnames(dta)[-1], sep = "") #base r to add prefix

#Inspect the result
###each sample will have a number of OTUs - some will be positive
###for the family in question)
###three different ways to do the same
###two in base R
aggregate(. ~ smplID, data = dta_dist_id, FUN = table) #gives obs where all are 1

describe(dta_dist_id)
#all have lachnosper
#all have proteobact
#one does not have clostridiac
#one does not have streptococ

rm(taxIDs, is_it_there)

#Load and merge sample data
source("3_clean_sampleData.R")

tmps <- new.env()
tmps$pckg <- c("tidyverse", "rms")
lapply(tmps$pckg, library, character.only = TRUE)
sapply(tmps$pckg, packageVersion)

#tmps$mapfile <- dta[dta$dumChr == 1, "smplID"] %>%
# distinct() %>%
#mutate(smplID = sub(".*Vtm *(.*?) *.fastq.*", "\\1",
#               .$smplID))#remove leading "Vtm" and tailing ".fastq"

#smpl.dta$dumChr <- factor(ifelse(smpl.dta$smplID %in%
#                                     tmps$mapfile$smplID, 1, 0),
#                        labels = c("no", "yes"))

tmps$varnms <- colnames(dta_dist_id[!grepl("smplID", colnames(dta_dist_id))])

colnames(dta_dist_id)

if_else(colnames(dta_dist_id) %in% tmps$varnms == 1, 1, 0)

df_joind <- dta_dist_id %>%
  mutate(smplID = sub(".*Vtm *(.*?) *.fastq.*", "\\1",
                      .$smplID)) %>% #remove leading "Vtm" and tailing ".fastq"
  right_join(smpl.dta, by = "smplID")

#smpl.dta[colnames(smpl.dta) %in% tmps$varnms][is.na(
# smpl.dta[colnames(smpl.dta) %in% tmps$varnms])] <- 0 #to replace NAs with zero which has already been done in l.84

###create variable to indicate if all the "good" bact are present (1) or not (0)
###Below is all the bacteria that has been found to correlate with the health outcome,
###however R throws error if a bacteria is included that nobody have. To work around that
###I have just "#" them out
df_joind <- df_joind %>%
  mutate(dum_one_of_healthy = if_else(dum_christensen == 1 |
                                               dum_akker == 1,
                                             1, 0),
         dum_all_healthy = if_else(dum_christensen == 1 &
                                            dum_akker == 1, 
                                          1, 0),
         dum_sum_healthy = factor(paste(dum_christensen,
                                        dum_akker,
                                        sep = "")),
         dum_one_of_obese = if_else(dum_streptococ == 1 |
                                      dum_lactobacillac == 1 |
                                      dum_proteobact == 1,
                                    1, 0),
         dum_all_obese = if_else(dum_streptococ == 1 &
                                   dum_lactobacillac == 1 &
                                   dum_proteobact == 1,
                                 1, 0),
         dum_sum_obese = factor(paste(dum_streptococ,
                                      dum_lactobacillac,
                                      dum_proteobact,
                                      sep = "")),
         dum_one_of_undrnurishd = if_else(dum_proteobact == 1 |
                                            dum_fusobact == 1, #|
                                          #dum_gallolyticus == 1 |
                                          #dum_smithii == 0,
                                          1, 0),
         dum_all_undrnurishd = if_else(dum_proteobact == 1 &
                                         dum_fusobact == 1, #&
                                       #dum_gallolyticus == 1 &
                                       #dum_smithii == 0,
                                       1, 0),
         dum_sum_undrnurishd = factor(paste(dum_proteobact,
                                            dum_fusobact,
                                            #dum_gallolyticus,
                                            sep = "")),
         dum_stunt = factor(cld.gwt, 
                            levels = c(1, 3), 
                            labels = c("stunted", "non-stunted")),
         dum_all_healthy_f = factor(dum_all_healthy,
                                    levels = c(1, 0),
                                    labels = c("all_healthy", "less_healthy")),
         mat_edu_f = factor(mat.edu)
         )

rm(dta_dist_id, bact, dta, smpl.dta)

#write_rds(df_joind, path = "data_work/df_joined200706.rds")

# ALL HEALTHY ------------------------

#tutorials
#WeightIt: https://cran.r-project.org/web/packages/WeightIt/vignettes/WeightIt_A0_basic_use.html
#cobalt: https://cran.r-project.org/web/packages/cobalt/vignettes/cobalt_A0_basic_use.html


#"reset" R
graphics.off() #unload graphics

###Function to unload all loaded packages
if(length(lapply(
  names(sessionInfo()$loadedOnly),
  library,
  character.only = TRUE) == NULL) != 0) {
  
  invisible(lapply(
    names(sessionInfo()$loadedOnly),
    library,
    character.only = TRUE)
  )
  
  invisible(lapply(
    paste0('package:', names(sessionInfo()$otherPkgs)),
    detach,
    character.only = TRUE, unload = TRUE, force = TRUE)
  )
}

tmps <- new.env()
tmps$pckg <- c("tidyverse", "WeightIt", "cobalt", "rms", "Hmisc")
lapply(tmps$pckg, library, character.only = TRUE)
sapply(tmps$pckg, packageVersion)

dta_joind <- read_rds("data_work/df_joined200706.rds")

#Use RMS to do logistic model of the outcome
dd <- datadist(dta_joind)
options(datadist = "dd")

PSAnalysis <- function(outcome, weight_vars) {
  a <- paste(outcome, "~", sep = " ")
  b <- paste(weight_vars, collapse = " + ")
  reg.formula <- as.formula(paste(a, b, sep = " "))
  
  balance_table1 <- bal.tab(reg.formula,
                    data = dta_joind,
                    estimand = "ATE",
                    m.threshold = .05,
                    binary = "std")
  
  ###gives standardized mean differences for continuous variables and differences
  ###in proportion for binary variables
  
  ###at a threshold of 0.05 the sex is deemed balanced in those who has all 
  ###the healthy bacteria and those who don't. See if even better balance will 
  ###be obtained with ps weights
  wght_out <- weightit(reg.formula,
                       data = dta_joind,
                       #focal = "all_healthy",
                       method = "ps",
                       missing = "ind",
                       estimand = "ATE")
  
  #weights_out
  smry_wght <- summary(wght_out)
  
  ###Weights with low variability are desirable because they improve the 
  ###precision of the estimator. This variability is presented in several 
  ###ways: by the ratio of the largest weight to the smallest in each group, 
  ###the coefficient of variation (standard deviation divided by the mean) of 
  ###the weights in each group, and the effective sample size computed from 
  ###the weights. We want a small ratio, a smaller coefficient of variation, 
  ###and a large effective sample size (ESS). Best used when comparing 
  ###weighting methods
  
  balance_table2 <- bal.tab(wght_out, 
                            m.threshold = .05, 
                            disp.v.ratio = TRUE,
                            binary = "std") #variance ratio ~1 indicate good balance (Austin 2009)
  
  lp <- love.plot(wght_out,
            #stats = c("variance.ratio", "mean.diff"),
            #threshold = c(variance.ratio = 1, mean.diff = .1),
            binary = "std",
            #threshold = .05,
            abs = TRUE,
            var.order = "unadjusted", #if left empty/=NULL the order will be as in the original data set
            drop.distance = TRUE,
            drop.missing = FALSE
            #, shapes = c("circle", "triangle"), #it's a ggplot so all the usual parameters can be used
            #colour = c("red", "blue"),
            #var.names = new.names #create a data.frame with two columns: "old" and "new" names
  )
  
  ###estimate average treatment effect in the treated (ATT)
  m.ate <- lrm(dum_stunt ~ dum_all_healthy_f,
               data = dta_joind,
               weights = wght_out$weights,
               x = TRUE,
               y = TRUE)
  
  sumry_m.ate <- summary(m.ate)
  
  
  nom.ate <- nomogram(m.ate, dum_all_healthy = c("all_healthy", "less_healthy"),
                      fun = plogis, funlabel = "Prob[Y=1]")
  ###plogis = logit / 1 + logit = exp(log(odds) / 1 + log(odds))
  
  to_return <- list(smry_wght, 
                    balance_table1, 
                    balance_table2, 
                    lp,
                    sumry_m.ate,
                    m.ate,
                    plot(sumry_m.ate),
                    plot(nom.ate))
  return(to_return)
}


#Examine and create weights
###Model 1
PSAnalysis("dum_all_healthy", "sex.f")

###Model 2
PSAnalysis("dum_all_healthy", c("sex.f", "mat_edu_f"))

#-*- -*- -*- OBS NOGET GALT: hvordan lave splines i bal.tab????
###Model 3
PSAnalysis("dum_all_healthy", c("rcs(psyc.smry)", "c.dis", "rcs(tot.carb)",
                                "rcs(tot.fbr)", "rcs(vitD)", "rcs(hb)", 
                                "rcs(ferritin)"))









