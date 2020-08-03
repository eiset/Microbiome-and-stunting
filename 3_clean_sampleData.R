# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Author:   Andreas Halgreen Eiset, eiset@ph.au.dk
# Title:    Clean sample data informed by 2b
# Licence:  GNU GPLv3
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

source("1b_loadTidy_sampleData.R")

tmps <- new.env()
tmps$pckg <- c("tidyverse")
lapply(tmps$pckg, library, character.only = TRUE)
sapply(tmps$pckg, packageVersion)



# Learned from 2b_0_... ------------------------------------
##drop the variable with no observation of outcome and make birth
##mode "assistet" -> "vaginal"
smpl.dta <- smpl.dta %>% 
  select(- c(astma, disblty, brst.fed6m)) %>% 
  filter(cld.gwt %in% c("1", "3")) %>% 
  mutate(b.mode = replace(b.mode, which(b.mode == "assisted"), 
                          "vaginal"))

##Remove outliers in fibres and carbs
###OBS KAN BRUGE?
#mutate(tot.fbr = replace(tot.fbr, 
#                        which(tot.fbr > quantile(tot.fbr, .96)),
#                       "NA"))
tmps$ids.fbr <- as.character(c(10409, 11416, 50917, 51405))
smpl.dta$tot.fbr <- as.numeric(
  ifelse(smpl.dta$smplID %in% tmps$ids.fbr,
         "NA",
         smpl.dta$tot.fbr))

tmps$ids.carb <- as.character(22117)
smpl.dta$tot.carb <- as.numeric(
  ifelse(smpl.dta$smplID %in% tmps$ids.carb,
         "NA",
         smpl.dta$tot.carb))



# factorise relevant vars ----------------------------------

tmps$tmp.dum.vars <- c("c.dis", 
                       "mat.edu",
                       "num.preg",
                       "sibs",
                       "brst.fed.excl",
                       "cld.gwt")

smpl.dta[ ,tmps$tmp.dum.vars] <- lapply(smpl.dta[tmps$tmp.dum.vars], 
                                        as.integer)
smpl.dta$sex.f <- factor(ifelse(as.numeric(smpl.dta$sex) == 1, 1, 2),
                         labels = c("male", "female"))
smpl.dta$b.mode.f <- factor(smpl.dta$b.mode)
###in dplyr use recode() and if_else()


# orthogonalize fibre and carbs -----------------------------
#plot(smpl.dta$tot.fbr, smpl.dta$tot.carb)

smpl.dta$cabfib.avg <- (smpl.dta$tot.fbr + smpl.dta$tot.carb) / 2
smpl.dta$cabfib.ratio <- smpl.dta$tot.fbr / smpl.dta$tot.carb
  
#plot(smpl.dta$cabfib.avg, smpl.dta$cabfib.ratio)


# psyc.soc -------------------------------------------------
###since none of the items seems redundant one way is to simply
###summarise the score. I.e. "higher score = better"
tmps$hugkis <- smpl.dta %>% 
  select(smplID,
         book,
         matches("^explain{1}[a-z]$")) %>% 
  mutate(psyc.smry = rowSums(select(., -smplID))) %>% 
  select(psyc.smry, smplID)

###combine the summary measure with the sample data
smpl.dta <- left_join(smpl.dta, tmps$hugkis, by = "smplID")

###remove all the individual psyc vars (not the summary)
tmps$psycvars <- c("^cvi", "^toy", "^expla", "book")
smpl.dta <- smpl.dta[, !grepl(paste(tmps$psycvars,
                                    collapse = "|"),
                              names(smpl.dta))]

#ggplot(smpl.dta, aes(psyc.smry)) +
 # geom_density()


