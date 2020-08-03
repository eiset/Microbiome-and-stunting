# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Author:   Andreas Halgreen Eiset, eiset@ph.au.dk
# Title:    Load and tidy sample data for UM microbiome project
# Licence:  GNU GPLv3
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

###I have opened and saved each file as rds first using 
###libraries "haven" (because of .dta file) and "readr":
###write_rds(df, path, compress = "xz")

tmps <- new.env()
tmps$pckg <- c("tidyverse", "rms")
lapply(tmps$pckg, library, character.only = TRUE)
sapply(tmps$pckg, packageVersion)

#read in the files containing additional data (nutrition, 
#delivery mode etc) and remove irrelevant psyc.soc vars
tmps$addta <- read_rds("comprsd_additionalDataCorrespDAG010819.rds") %>% 
  mutate(smplID = as.character(id), id = NULL) %>% 
  select(-matches("^toy{1}[0-9]$"),
         -grep(paste(c("^cvilnt{1}[a-z][0-9]$", "^cviolent{1}[a-z]$"),
                     collapse = "|"), colnames(.)))
###grepl kan ikke bruges i select argument



#Add microbiome and more nutrition data
smpl.dta <- read_rds("comprsd_mapfileNutData.rds") %>% 
  rename(smplID = "#SampleID", 
         cld.gwt = childgrowth, 
         b.mode = Description, 
         tot.fbr = Total_fibre, 
         vitD = Vitamin_D) %>% 
  select(1, 4:6, tot.fbr) %>% 
  mutate(smplID = sub(".*Vtm *(.*?) *.fastq.*", "\\1", 
                      .$smplID)) %>% 
  left_join(tmps$addta, by = "smplID")


#rename
smpl.dta <- smpl.dta %>% 
  rename(c.dis = chronicdisease, 
         astma = asthma, 
         disblty = disability, 
         hb = haemoglobin, 
         tot.carb = Total_carbohydrate, 
         mat.edu = maternaleducation, 
         num.preg = numberofpregnancies, 
         sibs = numberoflivekids, 
         brst.fed6m = breastfedat6months, 
         brst.fed.excl = exclusivelybreastfed, 
         week.compl.food = agecomplementfoodintroweeks, 
         vitD = vitamind)


rm(tmps)
 

