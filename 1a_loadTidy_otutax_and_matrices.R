# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Author:   Andreas Halgreen Eiset, eiset@ph.au.dk
# Title:    Data set for propensity score model
# Licence:  GNU GPLv3
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

pckg <- c("tidyverse", "phyloseq", "taxa", "metacoder", "vegan")
lapply(pckg, library, character.only = TRUE)
sapply(pckg, packageVersion)

#using phyloseq
###load otu and taxa data
dta <- import_biom("vietnam16S_otu_table_biomCOPY.biom")
#otu_table(dta)[1:5, 1:5]
###ok. The "VtmXXXXX.fastq" is the ID. The "denovoXXXX" is the OTU ID
#tax_table(dta)[1:5, 1:5]
#str(tax_table(dta))
###anoyinig notation. "RankX" is the taxanomic rank. "denovoXXXX" is the OTU.
###All bacteria names (at all tax levels) begins with "x__". Some data management
###is nessessary.

#-*- -*- -*-TIDY UP-*- -*- -*-
###using base R and dplyr (jf. grundwaldlab)
###phyloseq package format is too anoying. Make it matrix

otuTbl <- matrix(data = otu_table(dta), ncol = 100,
                 dimnames = dimnames(otu_table(dta)))
taxTbl <- matrix(data = tax_table(dta), ncol = 7,
                 dimnames = dimnames(tax_table(dta)))
colnames(taxTbl) <- c("Kingdom", "Phylum", "Class", "Order", "Family",
                      "Genus", "Species")
rm(dta)
detach("package:phyloseq", unload = TRUE)
tmps <- new.env()

###now get rid of anoying prefix after creating variable with whole
###string
tmp.taxmy <- apply(taxTbl, 1, paste, collapse = ";")
taxTbl <- gsub("\\w__", "", taxTbl)
taxTbl <- cbind(taxTbl, taxmy  = tmp.taxmy)

###make into tibbles and combine
tmp.otuTbl <- data.frame(otuTbl) %>%
  rownames_to_column(var = "otuID") %>%
  tbl_df()
tmp.taxTbl <- data.frame(taxTbl) %>%
  select(taxmy) %>%
  rownames_to_column(var = "otuID") %>%
  tbl_df()
otutax.raw <- left_join(tmp.otuTbl, tmp.taxTbl, by = "otuID")
rm(tmp.otuTbl, tmp.taxTbl, tmp.taxmy)
###The otuID is a unique id for each instance of a given species
###(or at another taxon level). That is, the same individual will
###often have several instances of the same species (or other
###level), however they will all have a unique denovoXXXX number.
###When converting to taxa object the created taxon id ("taxon_id")
###will relate to the deepest common level of any finding and assign
###a common taxon id to this i.e. all of the same findings (on whatever
###taxon level) will have the same taxon id. Thus, there are fewer
###unique taxon ids than otu ids. When doing quality assessment it
###must be done on otuID level since this is indicates where the
###PCR may have introduced errors


#using taxa package
###make otutax into taxa object
otutax <- parse_tax_data(otutax.raw,
                         class_cols = "taxmy",
                         class_sep = ";",
                         class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                         class_key = c("tax_rank" = "taxon_rank",
                                       "name" = "taxon_name"))
###if feel like it, drop class_data table: takes up space in the output and
###the same info is available from taxon_tanks()
otutax$data$class_data <- NULL

###"tax_data" is not very informative. Change name to represent OTU counts
###names(otutax$data) <- c("otuCnts", "class_data")
names(otutax$data) <- c("otuCnts")

###remove non-informative taxa ("" and "NA")
otutax <- filter_taxa(otutax, taxon_names != "", taxon_names != "NA")
###the OTUs that were associated with the nameless and "NA" taxa are reassigned
###to the closest supertaxon that passes the filter

###Are there any otu left with no counts after subsetting?
###Need ID from sample data set ("mapfile")
mapfile <- read_csv("mappingfile_vietnam16S_AMU_COPY.csv")
IDs <- mapfile %>%
  select(smplID = "#SampleID")
IDs <- IDs$smplID
ifelse(sum(rowSums(otutax$data$otuCnts[, IDs]) == 0) == 0,
       "no zero counts", "YES: go on to filter out below")
###Filter out the zero count OTUs (BEWARE also deletes taxa with 
###no observed OTUs, if other tables in obj$data data for these 
###taxa are also removed unless drop_obs = FALSE)
otutax <- filter_obs(otutax, "otuCnts",
                     !rowSums(otutax$data$otuCnts[, IDs]) == 0,
                     drop_taxa = TRUE)

#Create cut data set same as Andies (SEE "...alpha_and_beta..." FILE, search for cut off)
mapfile <- mapfile %>%
  select(1:10) %>%
  rename(smplID = "#SampleID",
         cld.gwt = childgrowth,
         b.mode = Description) %>%
  filter(!cld.gwt == 0) %>% #there is one obs. where the outcome has not been recorded (neither has most of the other variables)
  mutate(cld.gwt = factor(cld.gwt))

tmps$depth_dta <- data.frame(csum = colSums(otutax$data$otuCnts[, mapfile$smplID])) %>%
  rownames_to_column("sampleID") %>%
  as_tibble() %>%
  arrange(csum)

tmps$idsAndie <- head(tmps$depth_dta$sampleID, 8)

otutax$data$otuCnts_cut <- otutax$data$otuCnts %>%
  select(- one_of(tmps$idsAndie))


#-*- -*- -*-QUALITY CONTROL-*- -*- -*-
###From grunwaldlab: Sequencing technologies all have some amount
###of error. Some of this error happens during PCR and some happens
###during sequencing. There are lots of ways to filter out errors
###including: Removing sequences/bases that the sequencer indicates
###are low quality. This information is typically in the .fastq
###files returned by the sequencer. There are many programs, such
###as trimmomatic, that can use the quality scores in these files
###to filter out low-quality sequences. Clustering similar
###sequences together. Many erroneous sequences are very similar
###to the true sequences, so they will be incorporated into the
###same OTU as the true sequence. In that sense, clustering into
###OTUs will hide some sequencing error. ###Removal of chimeric
###sequences. “Chimeras” are more significant errors that occur
###during PCR when an incomplete amplicon acts as a primer for a
###different template in a subsequent cycle. Most amplicon
###metagenomic pipelines (e.g. QIIME, mothur, usearch) have
###functions to detect and remove chimeras.

#Removing low-abundance sequences. Most of the time, erroneous sequences will
###occur much less often than the true sequence. Simply removing any unique
###sequence / OTU that appears less than some minimum number of times is an
###effective way to remove these errors. OTUs represented by only a single
###sequence are commonly called singletons. So when you hear things like
###“singletons were removed” and means all OTUs composed of only a single
###sequence were removed. 
###Rarefying read counts. Usually, we try to make each
###sample in a sequencing run have the same number of reads sequenced, but that
###does not always happen. When a sample has many more reads than another
###sample, its apparent diversity can be artificially inflated since rare taxa
###are more likely to be found. To avoid this, the reads are subsampled to a
###fixed number or “rarefied”. This technique can throw out a lot of data, so
###it is not always appropriate (McMurdie and Holmes (2014)). 
###Converting counts to presence/absence. Due to PCR biases and the nature of compositional
###data, many researchers advocate not using read count information at all.
###Instead, they recommend converting counts to simply “present” or “absent”.

#Using "metacoder"
#Remove low-abundance counts
###First set all counts <3 (i.e. counts of 1 or 2) to 0 i.e.
###remove "doubletons"
otutax$data$otuCnts_zerod3 <- zero_low_counts(otutax, "otuCnts",
                                             min_count = 3,
                                             other_cols = TRUE)
###BEAWARE I have set min n = 3. n seems to be rather abritrately
###set. I have seen examples where n is set to 10.

table(rowSums(otutax$data$otuCnts_zerod3[IDs]) == 0)

###now remove OTUs and taxa that are 0 (equvivalent to remove
###non-informative taxa above)
###Check if thare are any
ifelse(sum(rowSums(otutax$data$otuCnts_zerod3[IDs]) == 0) == 0,
       "no zero counts", "YES: go on to filter out below")

###Filter out the zero count OTUs (BEWARE also deletes taxa with
###no observed OTUs, if other tables in obj$data data for these
###taxa are also removed unless drop_obs = FALSE)
otutax_trimd3 <- filter_obs(otutax, c("otuCnts", "otuCnts_cut"),
                           !rowSums(otutax$data$otuCnts_zerod3[IDs]) == 0,
                           drop_taxa = TRUE)

otutax_trimd3$data$otuCnts_zerod3 <- NULL #no use for this after this point

###same with trim for n = 10
otutax$data$otuCnts_zerod10 <- zero_low_counts(otutax, "otuCnts",
                                               min_count = 10,
                                               other_cols = TRUE)

otutax_trimd10 <- filter_obs(otutax, c("otuCnts", "otuCnts_cut"),
                            !rowSums(otutax$data$otuCnts_zerod10[IDs]) == 0,
                            drop_taxa = TRUE)

otutax_trimd10$data$otuCnts_zerod3 <- NULL
otutax_trimd10$data$otuCnts_zerod10 <- NULL

###in this data set the otuCnts will be a trimmed down version of
###the full data set according to the above criteria for n

table(rowSums(otutax_trimd3$data$otuCnts[IDs]) == 0)
table(rowSums(otutax_trimd10$data$otuCnts[IDs]) == 0)

#Rarefaction
###Rarefaction is used to simulate even numbers of reads per sample.
###Even sampling is important for at least two reasons:
###1. When comparing diversity of samples, more samples make it more
###likely to observe rare species. This will have a larger effect
###on some diversity indexes than others, depending on how they
###weigh rare species. 2.When comparing the similarity of samples,
###the presence of rare species due to higher sampling depth in one
###sample but not another can make the two samples appear more
###different than they actually are. Therefore, when comparing the
###diversity or similarity of samples, it is important to take into
###account differences in sampling depth. Rarefying can waste a
###lot of data and is not needed in many cases. See McMurdie and
###Holmes (2014) for an in-depth discussion on alternatives to
###rarefaction. We cover it here because it is a popular technique.
###Typically, the rarefaction depth chosen is the minimum sample
###depth. If the minimum depth is very small, the samples with the
###smallest depth can be removed and the minimum depth of the
###remaining samples can be used.

###OBS jeg har slettet mapfile data i funktionen
#hist(colSums(otutax$data$otuCnts[, mapfile$smplID]))
#summary(colSums(otutax$data$otuCnts[, mapfile$smplID]))

###Jeg venter med dette skridt af ovennævnte grund

#RAREFACTION FROM PHYLOSEQ:

#Standardise abundances to the median sequencing depth
###(from phyloseq tutorial: joey711.github.io/phyloseq/preprocess).
###This means for every sample (every entry): divide by the sum of
###counts for the same individual and multiply with the overall
###median seq depth

F_stdize <- function(df, df_nbr) {
  df <- df[["data"]]
  x <- case_when(df_nbr == "otuCnts" ~ 1,
                 df_nbr == "otuCnts_cut" ~ 2,
                 TRUE ~ NA_real_)
  y <- tibble(df)[[1]][[x]]
  m <- keep(y, is.numeric) %>%
    map(~ sum(.)) %>%
    unlist() %>%
    median()
  d <- y %>% modify_if(is.numeric, ~ round(m * (. / sum(.))))
  return(d)
}

otutax_trimd3$data$otuCnts_stdizd <- F_stdize(otutax_trimd3,
                                              "otuCnts")
otutax_trimd3$data$otuCnts_cut_stdizd <- F_stdize(otutax_trimd3,
                                                  "otuCnts_cut")

otutax_trimd10$data$otuCnts_stdizd <- F_stdize(otutax_trimd10,
                                               "otuCnts")
otutax_trimd10$data$otuCnts_cut_stdizd <- F_stdize(otutax_trimd10,
                                                   "otuCnts_cut")

#Filter taxa using cutoff of 3.0 for the coef. of variation

F_raref <- function(df, df_nbr){
  df <- df[["data"]]
  x <- case_when(df_nbr == "otuCnts_stdizd" ~ 3,
                 df_nbr == "otuCnts_cut_stdizd" ~ 4,
                 TRUE ~ NA_real_)
  y <- tibble(df)[[1]][[x]]
  y <- gather(y, key = smplID,
               value = obs, -c("taxon_id", "otuID", "taxmy")) %>%
    group_by(otuID) %>%
    mutate_at("obs",
              .funs = list(obs_cov = ~ sd(.) / mean(., na.rm = TRUE))) %>%
    ungroup() %>%
    filter(obs_cov > 3.0) %>%
    spread(key = smplID, value = obs)
  return(y)
}
###mutate_at og filter trinnet kan slås sammen til 1 filter trin
###men jeg har gjort det som ovenstående for at kunne se hvad der
###sker

otutax_trimd3$data$otuCnts_rared <- F_raref(otutax_trimd3,
                                            "otuCnts_stdizd")
otutax_trimd3$data$otuCnts_cut_rared <- F_raref(otutax_trimd3,
                                                "otuCnts_cut_stdizd")

otutax_trimd10$data$otuCnts_rared <- F_raref(otutax_trimd10,
                                             "otuCnts_stdizd")
otutax_trimd10$data$otuCnts_cut_rared <- F_raref(otutax_trimd10,
                                                 "otuCnts_cut_stdizd")
###From wikipedia: In probability theory and statistics, the
###coefficient of variation (CV), also known as relative standard
###deviation (RSD), is a standardized measure of dispersion of a
###probability distribution or frequency distribution. It is often
###expressed as a percentage, and is defined as the ratio of the
###standard deviation "sigma" to the mean "mu" (or its absolute
###value, | "mu" |. It shows the extent of variability in relation
###to the mean of the population. The coefficient of variation
###should be computed only for data measured on a ratio scale, as
###these are the measurements that allow the division operation.
###widely used in analytical chemistry to express the precision and
###repeatability of an assay.


#Proportions
###An easy alternative to rarefaction is to convert read counts to
###proportions of reads per-sample. This makes samples directly comparable,
###but does not solve the problem with inflated diversity estimates discussed
###in the rarefaction section.

###On the raw data set
otutax$data$otuProps <- calc_obs_props(otutax, "otuCnts",
                                       other_cols = TRUE)

###On the trimmed and cut data set (not the rarefyed)
F_propi <- function(df, dtaset){
  y <- calc_obs_props(df, dtaset, other_cols = TRUE)
  return(y)
}

otutax_trimd3$data$otuProps <- F_propi(otutax_trimd3,
                                       "otuCnts")
otutax_trimd3$data$otuProps_cut <- F_propi(otutax_trimd3,
                                           "otuCnts_cut")
otutax_trimd3$data$otuProps_stdizd <- F_propi(otutax_trimd3,
                                              "otuCnts_stdizd")
otutax_trimd3$data$otuProps_cut_stdizd <- F_propi(otutax_trimd3,
                                                  "otuCnts_cut_stdizd")

otutax_trimd10$data$otuProps <- F_propi(otutax_trimd10,
                                        "otuCnts")
otutax_trimd10$data$otuProps_cut <- F_propi(otutax_trimd10,
                                            "otuCnts_cut")
otutax_trimd10$data$otuProps_stdizd <- F_propi(otutax_trimd10,
                                               "otuCnts_stdizd")
otutax_trimd10$data$otuProps_cut_stdizd <- F_propi(otutax_trimd10,
                                                   "otuCnts_cut_stdizd")


#Presence/abscence reads (i.e. "Qualitative")
###Many researchers believe that read depth is not informative due
###to PCR and sequencing biases. Therefore, instead of comparing
###read counts, the counts can be converted to presence or absence
###of an OTU in a given sample. This can be done like so
otutax$data$otuPresAbsc <- counts_to_presence(otutax,
                                              "otuCnts",
                                              other_cols = TRUE)

###On the trimmed and cut data
F_propi <- function(df, dtaset){
  y <- counts_to_presence(df, dtaset, other_cols = TRUE)
  return(y)
}

otutax_trimd3$data$otuPresAbsc <- F_propi(otutax_trimd3,
                                          "otuCnts")
otutax_trimd3$data$otuPresAbsc_cut <- F_propi(otutax_trimd3,
                                              "otuCnts_cut")
otutax_trimd3$data$otuPresAbsc_stdizd <- F_propi(otutax_trimd3,
                                                 "otuCnts_stdizd")
otutax_trimd3$data$otuPresAbsc_cut_stdizd <- F_propi(otutax_trimd3,
                                                     "otuCnts_cut_stdizd")

otutax_trimd10$data$otuPresAbsc <- F_propi(otutax_trimd10,
                                           "otuCnts")
otutax_trimd10$data$otuPresAbsc_cut <- F_propi(otutax_trimd10,
                                               "otuCnts_cut")
otutax_trimd10$data$otuPresAbsc_stdizd <- F_propi(otutax_trimd10,
                                                  "otuCnts_stdizd")
otutax_trimd10$data$otuPresAbsc_cut_stdizd <- F_propi(otutax_trimd10,
                                                      "otuCnts_cut_stdizd")

