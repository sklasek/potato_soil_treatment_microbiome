phyloseq data cleanup
================
Scott Klasek
2023-12-11

### Purpose:

Clean up and add metadata to 16S and ITS phyloseq objects for all states
on the platform, for the second round of sequencing (2021 and 2022
samples). This document is analogous to 03_phyloseq_data_cleanup.Rmd.
Revisited once we resequenced a lot of second round samples that didn’t
have enough sequences (for that, see doc \#33, “Not enough sequences”.).

### Startup

#### Load libraries

``` r
packages <- c("tidyverse", "phyloseq")
invisible(lapply(packages, require, character.only = TRUE))
```

    ## Loading required package: tidyverse

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.1     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
    ## Loading required package: phyloseq

#### Load data

``` r
treatment.decoder <- read.csv(file="/Users/klas0061/Desktop/UMN/phyloseqs/treatment_decoder.csv") # the csv to match treatment #s with plot #s
treatment.decoder$plot <- as.numeric(treatment.decoder$plot)
sk_cats <- read.csv(file="/Users/klas0061/Desktop/UMN/phyloseqs/sk_categories.csv") # the csv to add treatment-specific metadata
larkin.decoder <- read.csv(file="/Users/klas0061/Desktop/UMN/phyloseqs/larkin_decoder.csv") # larkin-specific design needs a larkin-specific decoder file
```

#### Define functions

add_treatment_info uses the files treatment.decoder.csv and
sk_categories.csv to import treatment numbers to the phyloseq object
(corrected treatment numbers, so they are consistent with 1-6
corresponding to three-year and 7-12 to two-year rotations from ALL
states). Then from the treatment numbers, it finds the treatment
categories I assigned and imports those too. ITS database assigns all
taxonomy character strings as beginning with “\[a-z\]\_\_“, and removes
these three annoying characters.

``` r
add_treatment_info <- function(ps, state_abbr){
  see <- left_join(data.frame(sample_data(ps)), 
                 treatment.decoder %>% filter(state==state_abbr) %>% select(rotation, plot, treatment_new, national_control), 
                 by=c("rotation" = "rotation", "plot" = "plot")) # merges the dataframes by plot and rotation, add treatment_new and national_control columns
  see2 <- left_join(see, sk_cats %>% filter(state==state_abbr) %>% select(4:ncol(sk_cats)), by="treatment_new") # add in the treatment metadata
  rownames(see2) <- sample_names(ps) # add sample names to the dataframe
  sample_data(ps) <- see2 # write the new dataframe into the ps object
  return(ps)
}

fix_its_taxa <- function(ps){
  tax_table(ps)[, colnames(tax_table(ps))] <- gsub(tax_table(ps)[, colnames(tax_table(ps))], pattern = "[a-z]__", replacement = "")
  return(ps) # thanks! https://microbiome.github.io/tutorials/cleaning_taxonomy_table.html
}
```

### Cleaning up phyloseqs!

Organized in somewhat chronological order. Showing Idaho 16S and ITS
first, because none of their samples needed to be resequenced.

#### Idaho Obj1 samples are all complete in the first round. Processed 8-4-23.

Recall some of the 2019 samples were missing, and these have been added
to the 2021 phyloseq objects.  
16S.

``` r
# 2021 has some leftover 2019 samples in it
ID_16S_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/ID_16S_2021_processed.ps")

# add season as character- all are Fall in this case
sample_data(ID_16S_2021_processed.ps)$season <- "Fall" 

# add block
sample_data(ID_16S_2021_processed.ps)$block <- as.numeric(substring(sample_data(ID_16S_2021_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
ID_16S_2021_m.ps <- add_treatment_info(ID_16S_2021_processed.ps, "ID") # m signifies metadata-added

# final check and write out
view(sample_data(ID_16S_2021_m.ps))
saveRDS(ID_16S_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/ID_16S_2021_with_fall_2019_m.ps")

# 2022 
ID_16S_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/ID_16S_2022_processed.ps")

# take a look at the metadata
view(ID_16S_2022_processed.ps@sam_data)

# gotta omit this one sample because there were replicates and this was the one with fewer reads
sample_sums(ID_16S_2022_processed.ps)["ID-1-3-306-22-4_S1672"]
ID_16S_2022_processed.ps@sam_data["ID-1-3-306-22-4-REP2_S1675","month"] <- 4 # fix month
ID_16S_2022_processed.ps <- subset_samples(ID_16S_2022_processed.ps, 
                                           sample_names(ID_16S_2022_processed.ps) != "ID-1-3-306-22-4_S1672") # omit the sample

# add season as character- Spring or Summer in this case
# add season
sample_data(ID_16S_2022_processed.ps)$season <- ifelse(sample_data(ID_16S_2022_processed.ps)$month == 4, "Spring",
                                                       ifelse(sample_data(ID_16S_2022_processed.ps)$month == 6, "Summer", NA))

# add block
sample_data(ID_16S_2022_processed.ps)$block <- as.numeric(substring(sample_data(ID_16S_2022_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
ID_16S_2022_m.ps <- add_treatment_info(ID_16S_2022_processed.ps, "ID") # m signifies metadata-added

# final check and write out
view(sample_data(ID_16S_2022_m.ps))
saveRDS(ID_16S_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/ID_16S_2022_m.ps")
```

ITS. The only additional step is to fix taxonomy labels

``` r
# 2021 has some leftover 2019 samples in it
ID_ITS_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/ID_ITS_2021_processed.ps")

# add season as character- all are Fall in this case
sample_data(ID_ITS_2021_processed.ps)$season <- "Fall" 

# add block
sample_data(ID_ITS_2021_processed.ps)$block <- as.numeric(substring(sample_data(ID_ITS_2021_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
ID_ITS_2021_m.ps <- add_treatment_info(ID_ITS_2021_processed.ps, "ID") # m signifies metadata-added

# fix taxonomy labels
ID_ITS_2021_m.ps <- fix_its_taxa(ID_ITS_2021_m.ps)

# final check and write out
view(sample_data(ID_ITS_2021_m.ps))
saveRDS(ID_ITS_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/ID_ITS_2021_with_fall_2019_m.ps")

# 2022 
ID_ITS_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/ID_ITS_2022_processed.ps")

# take a look at the metadata
view(ID_ITS_2022_processed.ps@sam_data)

# gotta omit this one sample because there were replicates and this was the one with fewer reads
sample_sums(ID_ITS_2022_processed.ps)["ID-1-3-306-22-4_S1672"]
ID_ITS_2022_processed.ps@sam_data["ID-1-3-306-22-4-REP2_S1675","month"] <- 4 # fix month
ID_ITS_2022_processed.ps <- subset_samples(ID_ITS_2022_processed.ps, 
                                           sample_names(ID_ITS_2022_processed.ps) != "ID-1-3-306-22-4_S1672") # omit the sample

# add season as character- Spring or Summer in this case
# add season
sample_data(ID_ITS_2022_processed.ps)$season <- ifelse(sample_data(ID_ITS_2022_processed.ps)$month == 4, "Spring",
                                                       ifelse(sample_data(ID_ITS_2022_processed.ps)$month == 6, "Summer", NA))

# add block
sample_data(ID_ITS_2022_processed.ps)$block <- as.numeric(substring(sample_data(ID_ITS_2022_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
ID_ITS_2022_m.ps <- add_treatment_info(ID_ITS_2022_processed.ps, "ID") # m signifies metadata-added

# fix taxonomy labels
ID_ITS_2022_m.ps <- fix_its_taxa(ID_ITS_2022_m.ps)

# final check and write out
view(sample_data(ID_ITS_2022_m.ps))
saveRDS(ID_ITS_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/ID_ITS_2022_m.ps")
```

Sent this to Brenda and Gilbert back in August.

#### Colorado: Processed 12-12-23.

16S.

``` r
# 2021
CO_16S_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/CO_16S_2021_processed.ps") # import ps

# convert month to numeric and add season
sample_data(CO_16S_2021_processed.ps)$month <- as.numeric(sample_data(CO_16S_2021_processed.ps)$month) 
sample_data(CO_16S_2021_processed.ps)$season <- "Fall" # add season as character

# add block
sample_data(CO_16S_2021_processed.ps)$block <- NA # easy. CO has no blocks

# add in the treatment numbers and metadata
CO_16S_2021_m.ps <- add_treatment_info(CO_16S_2021_processed.ps, "CO") # m signifies metadata-added

# final check and write out
# sample_data(CO_16S_2021_m.ps)
saveRDS(CO_16S_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/CO_16S_2021_m.ps")


# 2022
CO_16S_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/CO_16S_2022_processed.ps") # import ps

# convert month to numeric and add season
sample_data(CO_16S_2022_processed.ps)$month <- as.numeric(sample_data(CO_16S_2022_processed.ps)$month)
sample_data(CO_16S_2022_processed.ps)$season <- ifelse(sample_data(CO_16S_2022_processed.ps)$month == 4, "Spring",
                                                       ifelse(sample_data(CO_16S_2022_processed.ps)$month == 7, "Summer", NA))

# add block
sample_data(CO_16S_2022_processed.ps)$block <- NA # easy. CO has no blocks

# add in the treatment numbers and metadata
CO_16S_2022_m.ps <- add_treatment_info(CO_16S_2022_processed.ps, "CO") # m signifies metadata-added

# final check and write out
# sample_data(CO_16S_2022_m.ps)
saveRDS(CO_16S_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/CO_16S_2022_m.ps")
```

ITS.

``` r
# 2021
CO_ITS_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/CO_ITS_2021_processed.ps") # import ps

# convert month to numeric and add season
sample_data(CO_ITS_2021_processed.ps)$month <- as.numeric(sample_data(CO_ITS_2021_processed.ps)$month)
sample_data(CO_ITS_2021_processed.ps)$season <- "Fall" 

# add block
sample_data(CO_ITS_2021_processed.ps)$block <- NA # easy. CO has no blocks

# add in the treatment numbers and metadata
CO_ITS_2021_m.ps <- add_treatment_info(CO_ITS_2021_processed.ps, "CO") # m signifies metadata-added

# fix taxonomy labels
CO_ITS_2021_m.ps <- fix_its_taxa(CO_ITS_2021_m.ps)

# final check and write out
# sample_data(CO_ITS_2021_m.ps)
saveRDS(CO_ITS_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/CO_ITS_2021_m.ps")


# 2022
CO_ITS_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/CO_ITS_2022_processed.ps") # import ps

# convert month to numeric and add season
sample_data(CO_ITS_2022_processed.ps)$month <- as.numeric(sample_data(CO_ITS_2022_processed.ps)$month)
sample_data(CO_ITS_2022_processed.ps)$season <- ifelse(sample_data(CO_ITS_2022_processed.ps)$month == 4, "Spring",
                                                       ifelse(sample_data(CO_ITS_2022_processed.ps)$month == 7, "Summer", NA))

# add block
sample_data(CO_ITS_2022_processed.ps)$block <- NA # easy. CO has no blocks

# add in the treatment numbers and metadata
CO_ITS_2022_m.ps <- add_treatment_info(CO_ITS_2022_processed.ps, "CO") # m signifies metadata-added

# fix taxonomy labels
CO_ITS_2022_m.ps <- fix_its_taxa(CO_ITS_2022_m.ps)

# final check and write out
# sample_data(CO_ITS_2022_m.ps)
saveRDS(CO_ITS_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/CO_ITS_2022_m.ps")
```

Uploaded phyloseqs to Box, Jorge and Jane notified 12-12.

#### Maine: Processed 12-12-23.

16S.

``` r
# 2021
ME_16S_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/ME_16S_2021_processed.ps") # import ps

# convert month to numeric and add season
sample_data(ME_16S_2021_processed.ps)$month <- as.numeric(sample_data(ME_16S_2021_processed.ps)$month) 
sample_data(ME_16S_2021_processed.ps)$season <- "Fall" # add season as character

# add block (first # of plot #)
sample_data(ME_16S_2021_processed.ps)$block <- as.numeric(substring(sample_data(ME_16S_2021_processed.ps)$plot, 1,1)) 

# add in the treatment numbers and metadata
ME_16S_2021_m.ps <- add_treatment_info(ME_16S_2021_processed.ps, "ME") # m signifies metadata-added

# final check and write out
# sample_data(ME_16S_2021_m.ps)
saveRDS(ME_16S_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/ME_16S_2021_m.ps")


# 2022
ME_16S_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/ME_16S_2022_processed.ps") # import ps

# convert month to numeric and add season
sample_data(ME_16S_2022_processed.ps)$month <- as.numeric(sample_data(ME_16S_2022_processed.ps)$month)
sample_data(ME_16S_2022_processed.ps)$season <- ifelse(sample_data(ME_16S_2022_processed.ps)$month == 5, "Spring",
                                                       ifelse(sample_data(ME_16S_2022_processed.ps)$month == 7, "Summer", NA))

# add block
sample_data(ME_16S_2022_processed.ps)$block <- as.numeric(substring(sample_data(ME_16S_2022_processed.ps)$plot, 1,1)) 

# add in the treatment numbers and metadata
ME_16S_2022_m.ps <- add_treatment_info(ME_16S_2022_processed.ps, "ME") # m signifies metadata-added

# final check and write out
# sample_data(ME_16S_2022_m.ps)
saveRDS(ME_16S_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/ME_16S_2022_m.ps")
```

ITS.

``` r
# 2021
ME_ITS_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/ME_ITS_2021_processed.ps") # import ps

# convert month to numeric and add season
sample_data(ME_ITS_2021_processed.ps)$month <- as.numeric(sample_data(ME_ITS_2021_processed.ps)$month)
sample_data(ME_ITS_2021_processed.ps)$season <- "Fall" 

# add block
sample_data(ME_ITS_2021_processed.ps)$block <- as.numeric(substring(sample_data(ME_ITS_2021_processed.ps)$plot, 1,1)) 

# add in the treatment numbers and metadata
ME_ITS_2021_m.ps <- add_treatment_info(ME_ITS_2021_processed.ps, "ME") # m signifies metadata-added

# fix taxonomy labels
ME_ITS_2021_m.ps <- fix_its_taxa(ME_ITS_2021_m.ps)

# final check and write out
# sample_data(ME_ITS_2021_m.ps)
saveRDS(ME_ITS_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/ME_ITS_2021_m.ps")


# 2022
ME_ITS_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/ME_ITS_2022_processed.ps") # import ps

# convert month to numeric and add season
sample_data(ME_ITS_2022_processed.ps)$month <- as.numeric(sample_data(ME_ITS_2022_processed.ps)$month)
sample_data(ME_ITS_2022_processed.ps)$season <- ifelse(sample_data(ME_ITS_2022_processed.ps)$month == 5, "Spring",
                                                       ifelse(sample_data(ME_ITS_2022_processed.ps)$month == 7, "Summer", NA))

# add block
sample_data(ME_ITS_2022_processed.ps)$block <- as.numeric(substring(sample_data(ME_ITS_2022_processed.ps)$plot, 1,1)) 

# add in the treatment numbers and metadata
ME_ITS_2022_m.ps <- add_treatment_info(ME_ITS_2022_processed.ps, "ME") # m signifies metadata-added

# fix taxonomy labels
ME_ITS_2022_m.ps <- fix_its_taxa(ME_ITS_2022_m.ps)

# final check and write out
# sample_data(ME_ITS_2022_m.ps)
saveRDS(ME_ITS_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/ME_ITS_2022_m.ps")
```

Uploaded phyloseqs to Box, Katie and Jay notified 12-12.

#### Michigan: Processed 12-15-23.

16S.

``` r
# 2021
MI_16S_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/MI_16S_2021_processed.ps") # import ps

# convert month to numeric and add season
sample_data(MI_16S_2021_processed.ps)$month <- as.numeric(sample_data(MI_16S_2021_processed.ps)$month) 
sample_data(MI_16S_2021_processed.ps)$season <- "Fall" # add season as character

# add block (first # of plot #)
sample_data(MI_16S_2021_processed.ps)$block <- as.numeric(substring(sample_data(MI_16S_2021_processed.ps)$plot, 1,1)) 

# add in the treatment numbers and metadata
MI_16S_2021_m.ps <- add_treatment_info(MI_16S_2021_processed.ps, "MI") # m signifies metadata-added

# final check and write out
# sample_data(MI_16S_2021_m.ps)
saveRDS(MI_16S_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/MI_16S_2021_m.ps")


# 2022 
MI_16S_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/MI_16S_2022_processed.ps") # import ps

# add season, then covert month to numeric
sample_data(MI_16S_2022_processed.ps)$season <- ifelse(sample_data(MI_16S_2022_processed.ps)$month == 5, "Spring",
                                                       ifelse(sample_data(MI_16S_2022_processed.ps)$month == 8, "Summer", "Spring"))
sample_data(MI_16S_2022_processed.ps)[which(sample_data(MI_16S_2022_processed.ps)$month == "5-REP2"),"month"] <- 5
sample_data(MI_16S_2022_processed.ps)$month <- as.numeric(sample_data(MI_16S_2022_processed.ps)$month)

# add block
sample_data(MI_16S_2022_processed.ps)$block <- as.numeric(substring(sample_data(MI_16S_2022_processed.ps)$plot, 1,1)) 

# add in the treatment numbers and metadata
MI_16S_2022_m.ps <- add_treatment_info(MI_16S_2022_processed.ps, "MI") # m signifies metadata-added

# final check and write out
# sample_data(MI_16S_2022_m.ps)
saveRDS(MI_16S_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/MI_16S_2022_m.ps")
```

ITS.

``` r
# 2021
MI_ITS_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/MI_ITS_2021_processed.ps") # import ps

# convert month to numeric and add season
sample_data(MI_ITS_2021_processed.ps)$month <- as.numeric(sample_data(MI_ITS_2021_processed.ps)$month)
sample_data(MI_ITS_2021_processed.ps)$season <- "Fall" 

# add block
sample_data(MI_ITS_2021_processed.ps)$block <- as.numeric(substring(sample_data(MI_ITS_2021_processed.ps)$plot, 1,1)) 

# add in the treatment numbers and metadata
MI_ITS_2021_m.ps <- add_treatment_info(MI_ITS_2021_processed.ps, "MI") # m signifies metadata-added

# fix taxonomy labels
MI_ITS_2021_m.ps <- fix_its_taxa(MI_ITS_2021_m.ps)

# final check and write out
sample_data(MI_ITS_2021_m.ps)
saveRDS(MI_ITS_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/MI_ITS_2021_m.ps")


# 2022
MI_ITS_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/MI_ITS_2022_processed.ps") # import ps

# convert month to numeric and add season
sample_data(MI_ITS_2022_processed.ps)$month <- as.numeric(sample_data(MI_ITS_2022_processed.ps)$month)
sample_data(MI_ITS_2022_processed.ps)$season <- ifelse(sample_data(MI_ITS_2022_processed.ps)$month == 5, "Spring",
                                                       ifelse(sample_data(MI_ITS_2022_processed.ps)$month == 8, "Summer", NA))

# add block
sample_data(MI_ITS_2022_processed.ps)$block <- as.numeric(substring(sample_data(MI_ITS_2022_processed.ps)$plot, 1,1)) 

# add in the treatment numbers and metadata
MI_ITS_2022_m.ps <- add_treatment_info(MI_ITS_2022_processed.ps, "MI") # m signifies metadata-added

# fix taxonomy labels
MI_ITS_2022_m.ps <- fix_its_taxa(MI_ITS_2022_m.ps)

# final check and write out
# sample_data(MI_ITS_2022_m.ps)
saveRDS(MI_ITS_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/MI_ITS_2022_m.ps")
```

Uploaded phyloseqs to Box, Kurt notified 12-15.

#### Minnesota: Processed 12-15-23.

16S.

``` r
# 2021 
MN_16S_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/MN_16S_2021_processed.ps") # import ps

# inspect the sample data
view(data.frame(MN_16S_2021_processed.ps@sam_data))

# convert month to numeric and add season
sample_data(MN_16S_2021_processed.ps)$month <- as.numeric(sample_data(MN_16S_2021_processed.ps)$month) 
sample_data(MN_16S_2021_processed.ps)$season <- "Fall" 

# add block
sample_data(MN_16S_2021_processed.ps)$block <- as.numeric(substring(sample_data(MN_16S_2021_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
MN_16S_2021_m.ps <- add_treatment_info(MN_16S_2021_processed.ps, "MN") # m signifies metadata-added

# final check 
view(data.frame(MN_16S_2021_m.ps@sam_data))

# write out
saveRDS(MN_16S_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/MN_16S_2021_m.ps")


# 2022
MN_16S_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/MN_16S_2022_processed.ps") # import ps

# inspect
view(MN_16S_2022_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(MN_16S_2022_processed.ps)$month <- as.numeric(sample_data(MN_16S_2022_processed.ps)$month) 
sample_data(MN_16S_2022_processed.ps)$season <- ifelse(sample_data(MN_16S_2022_processed.ps)$month == 4, "Spring",
                                                ifelse(sample_data(MN_16S_2022_processed.ps)$month == 6, "Summer", NA)) 
# add block
sample_data(MN_16S_2022_processed.ps)$block <- as.numeric(substring(sample_data(MN_16S_2022_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
MN_16S_2022_m.ps <- add_treatment_info(MN_16S_2022_processed.ps, "MN") # m signifies metadata-added

# final check and write out
# view(sample_data(MN_16S_2022_m.ps))
saveRDS(MN_16S_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/MN_16S_2022_m.ps")
```

ITS.

``` r
# 2021
MN_ITS_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/MN_ITS_2021_processed.ps") # import ps

# inspect the sample data
view(data.frame(MN_ITS_2021_processed.ps@sam_data))

# convert month to numeric and add season
sample_data(MN_ITS_2021_processed.ps)$month <- as.numeric(sample_data(MN_ITS_2021_processed.ps)$month) 
sample_data(MN_ITS_2021_processed.ps)$season <- "Fall" 

# add block
sample_data(MN_ITS_2021_processed.ps)$block <- as.numeric(substring(sample_data(MN_ITS_2021_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
MN_ITS_2021_m.ps <- add_treatment_info(MN_ITS_2021_processed.ps, "MN") # m signifies metadata-added

# fix ITS taxa
MN_ITS_2021_m.ps <- fix_its_taxa(MN_ITS_2021_m.ps)

# final check 
# view(data.frame(MN_ITS_2021_m.ps@sam_data))

# write out
saveRDS(MN_ITS_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/MN_ITS_2021_m.ps")


# 2022
MN_ITS_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/MN_ITS_2022_processed.ps") # import ps

# inspect
# view(MN_ITS_2022_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(MN_ITS_2022_processed.ps)$month <- as.numeric(sample_data(MN_ITS_2022_processed.ps)$month) 
sample_data(MN_ITS_2022_processed.ps)$season <- ifelse(sample_data(MN_ITS_2022_processed.ps)$month == 4, "Spring",
                                                ifelse(sample_data(MN_ITS_2022_processed.ps)$month == 6, "Summer", NA)) 
# add block
sample_data(MN_ITS_2022_processed.ps)$block <- as.numeric(substring(sample_data(MN_ITS_2022_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
MN_ITS_2022_m.ps <- add_treatment_info(MN_ITS_2022_processed.ps, "MN") # m signifies metadata-added

# fix ITS taxa
MN_ITS_2022_m.ps <- fix_its_taxa(MN_ITS_2022_m.ps)

# final check and write out
# view(sample_data(MN_ITS_2022_m.ps))
saveRDS(MN_ITS_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/MN_ITS_2022_m.ps")
```

Uploaded phyloseqs to Box on 12-15.

#### North Dakota: Processed 12-14-23.

16S.

``` r
# 2021 
ND_16S_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/ND_16S_2021_processed.ps") # import ps

# inspect
view(ND_16S_2021_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(ND_16S_2021_processed.ps)$month <- as.numeric(sample_data(ND_16S_2021_processed.ps)$month) 
sample_data(ND_16S_2021_processed.ps)$season <- "Fall" # add season as character

# add block
sample_data(ND_16S_2021_processed.ps)$block <- as.numeric(substring(sample_data(ND_16S_2021_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
ND_16S_2021_m.ps <- add_treatment_info(ND_16S_2021_processed.ps, "ND") # m signifies metadata-added

# final check and write out
view(ND_16S_2021_m.ps@sam_data)
saveRDS(ND_16S_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/ND_16S_2021_m.ps")


# 2022 
ND_16S_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/ND_16S_2022_processed.ps")

# take a look at the metadata
view(ND_16S_2022_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(ND_16S_2022_processed.ps)$month <- as.numeric(sample_data(ND_16S_2022_processed.ps)$month) 
sample_data(ND_16S_2022_processed.ps)$season <- ifelse(sample_data(ND_16S_2022_processed.ps)$month == 5, "Spring",
                                                       ifelse(sample_data(ND_16S_2022_processed.ps)$month == 7, "Summer", NA))

# add block
sample_data(ND_16S_2022_processed.ps)$block <- as.numeric(substring(sample_data(ND_16S_2022_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
ND_16S_2022_m.ps <- add_treatment_info(ND_16S_2022_processed.ps, "ND") # m signifies metadata-added

# final check and write out
view(sample_data(ND_16S_2022_m.ps))
saveRDS(ND_16S_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/ND_16S_2022_m.ps")
```

ITS.

``` r
# 2021
ND_ITS_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/ND_ITS_2021_processed.ps") # import ps

# inspect
view(ND_ITS_2021_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(ND_ITS_2021_processed.ps)$month <- as.numeric(sample_data(ND_ITS_2021_processed.ps)$month) 
sample_data(ND_ITS_2021_processed.ps)$season <- "Fall" # add season as character

# add block
sample_data(ND_ITS_2021_processed.ps)$block <- as.numeric(substring(sample_data(ND_ITS_2021_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
ND_ITS_2021_m.ps <- add_treatment_info(ND_ITS_2021_processed.ps, "ND") # m signifies metadata-added

# fix taxonomy labels
ND_ITS_2021_m.ps <- fix_its_taxa(ND_ITS_2021_m.ps)

# final check and write out
view(ND_ITS_2021_m.ps@sam_data)
saveRDS(ND_ITS_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/ND_ITS_2021_m.ps")


# 2022 
ND_ITS_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/ND_ITS_2022_processed.ps")

# take a look at the metadata
view(ND_ITS_2022_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(ND_ITS_2022_processed.ps)$month <- as.numeric(sample_data(ND_ITS_2022_processed.ps)$month) 
sample_data(ND_ITS_2022_processed.ps)$season <- ifelse(sample_data(ND_ITS_2022_processed.ps)$month == 5, "Spring",
                                                       ifelse(sample_data(ND_ITS_2022_processed.ps)$month == 7, "Summer", NA))

# add block
sample_data(ND_ITS_2022_processed.ps)$block <- as.numeric(substring(sample_data(ND_ITS_2022_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
ND_ITS_2022_m.ps <- add_treatment_info(ND_ITS_2022_processed.ps, "ND") # m signifies metadata-added

# fix taxonomy labels
ND_ITS_2022_m.ps <- fix_its_taxa(ND_ITS_2022_m.ps)

# final check and write out
view(sample_data(ND_ITS_2022_m.ps))
saveRDS(ND_ITS_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/ND_ITS_2022_m.ps")
```

Uploaded phyloseqs to Box & Kim, Egla, and Julie notified 12-14.

#### Oregon: Processed 12-16-23.

Recall some of the 2019 samples were missing, and these have been added
to the 2021 phyloseq objects.  
16S.

``` r
# 2021 
OR_16S_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/OR_16S_2021_processed.ps") # import ps

# inspect
# sample_data(OR_16S_2021_processed.ps)

# add NAs for the blanks and compost sample metadata
sample_data(OR_16S_2021_processed.ps)[c(25,98:100),1:6] <- NA

# convert month to numeric and add season
sample_data(OR_16S_2021_processed.ps)$month <- as.numeric(sample_data(OR_16S_2021_processed.ps)$month) 
sample_data(OR_16S_2021_processed.ps)$season <- ifelse(sample_data(OR_16S_2021_processed.ps)$month == 4, "Spring",
                                                       ifelse(sample_data(OR_16S_2021_processed.ps)$month == 6, "Summer", 
                                                             ifelse(sample_data(OR_16S_2021_processed.ps)$month == 8, "Fall", NA))) 

# add block info (Oregon has plots numbered 100s to 800s for different rotations, an anomaly)
ordf <- data.frame(sample_data(OR_16S_2021_processed.ps))
ordf$block <- as.numeric(substring(sample_data(OR_16S_2021_processed.ps)$plot, 1,1)) # 3-yr rotation blocks are numbered as usual, the first digit of the plot #
ordf[which(ordf$rotation==2),"block"] <- (ordf[which(ordf$rotation==2),"block"]-4) # 2-yr rotation blocks begin with 5-8, so subtract 4
sample_data(OR_16S_2021_processed.ps) <- ordf

# add in the treatment numbers and metadata
OR_16S_2021_m.ps <- add_treatment_info(OR_16S_2021_processed.ps, "OR") # m signifies metadata-added

# final check and write out
# sample_data(OR_16S_2021_m.ps)
saveRDS(OR_16S_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/OR_16S_2021_with_some_2019_m.ps")


# 2022 
OR_16S_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/OR_16S_2022_processed.ps") # import ps
view(OR_16S_2022_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(OR_16S_2022_processed.ps)$month <- as.numeric(sample_data(OR_16S_2022_processed.ps)$month) 
sample_data(OR_16S_2022_processed.ps)$season <- ifelse(sample_data(OR_16S_2022_processed.ps)$month == 3, "Spring",
                                                ifelse(sample_data(OR_16S_2022_processed.ps)$month == 6, "Summer", NA)) 

# add block info (Oregon has plots numbered 100s to 800s for different rotations, an anomaly)
ordf <- data.frame(sample_data(OR_16S_2022_processed.ps))
ordf$block <- as.numeric(substring(sample_data(OR_16S_2022_processed.ps)$plot, 1,1)) # 3-yr rotation blocks are numbered as usual, the first digit of the plot #
ordf[which(ordf$rotation==2),"block"] <- (ordf[which(ordf$rotation==2),"block"]-4) # 2-yr rotation blocks begin with 5-8, so subtract 4
sample_data(OR_16S_2022_processed.ps) <- ordf

# add in the treatment numbers and metadata
OR_16S_2022_m.ps <- add_treatment_info(OR_16S_2022_processed.ps, "OR") # m signifies metadata-added

# final check and write out
# sample_data(OR_16S_2022_m.ps)
saveRDS(OR_16S_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/OR_16S_2022_m.ps")
```

ITS.

``` r
# 2021
OR_ITS_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/OR_ITS_2021_processed.ps") # import ps

# inspect
# sample_data(OR_ITS_2021_processed.ps)

# convert month to numeric and add season
sample_data(OR_ITS_2021_processed.ps)$month <- as.numeric(sample_data(OR_ITS_2021_processed.ps)$month) 
sample_data(OR_ITS_2021_processed.ps)$season <- ifelse(sample_data(OR_ITS_2021_processed.ps)$month == 4, "Spring",
                                                       ifelse(sample_data(OR_ITS_2021_processed.ps)$month == 6, "Summer", 
                                                             ifelse(sample_data(OR_ITS_2021_processed.ps)$month == 8, "Fall", NA))) 

# add block info (Oregon has plots numbered 100s to 800s for different rotations, an anomaly)
ordf <- data.frame(sample_data(OR_ITS_2021_processed.ps))
ordf$block <- as.numeric(substring(sample_data(OR_ITS_2021_processed.ps)$plot, 1,1)) # 3-yr rotation blocks are numbered as usual, the first digit of the plot #
ordf[which(ordf$rotation==2),"block"] <- (ordf[which(ordf$rotation==2),"block"]-4) # 2-yr rotation blocks begin with 5-8, so subtract 4
sample_data(OR_ITS_2021_processed.ps) <- ordf

# add in the treatment numbers and metadata
OR_ITS_2021_m.ps <- add_treatment_info(OR_ITS_2021_processed.ps, "OR") # m signifies metadata-added

# fix taxonomy labels
OR_ITS_2021_m.ps <- fix_its_taxa(OR_ITS_2021_m.ps)

# final check and write out
# sample_data(OR_ITS_2021_m.ps)
saveRDS(OR_ITS_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/OR_ITS_2021_with_some_2019_m.ps")


# 2022
OR_ITS_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/OR_ITS_2022_processed.ps") # import ps

# inspect
# sample_data(OR_ITS_2022_processed.ps)

# convert month to numeric and add season
sample_data(OR_ITS_2022_processed.ps)$month <- as.numeric(sample_data(OR_ITS_2022_processed.ps)$month) 
sample_data(OR_ITS_2022_processed.ps)$season <- ifelse(sample_data(OR_ITS_2022_processed.ps)$month == 3, "Spring",
                                                       ifelse(sample_data(OR_ITS_2022_processed.ps)$month == 6, "Summer", NA)) 

# add block info (Oregon has plots numbered 100s to 800s for different rotations, an anomaly)
ordf <- data.frame(sample_data(OR_ITS_2022_processed.ps))
ordf$block <- as.numeric(substring(sample_data(OR_ITS_2022_processed.ps)$plot, 1,1)) # 3-yr rotation blocks are numbered as usual, the first digit of the plot #
ordf[which(ordf$rotation==2),"block"] <- (ordf[which(ordf$rotation==2),"block"]-4) # 2-yr rotation blocks begin with 5-8, so subtract 4
sample_data(OR_ITS_2022_processed.ps) <- ordf

# add in the treatment numbers and metadata
OR_ITS_2022_m.ps <- add_treatment_info(OR_ITS_2022_processed.ps, "OR") # m signifies metadata-added

# fix taxonomy labels
OR_ITS_2022_m.ps <- fix_its_taxa(OR_ITS_2022_m.ps)

# final check and write out
# sample_data(OR_ITS_2022_m.ps)
saveRDS(OR_ITS_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/OR_ITS_2022_m.ps")
```

Uploaded phyloseqs to Box, Ken notified 12-16.

#### Wisconsin: Processed 12-16-23.

16S.

``` r
# 2021 
WI_16S_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/WI_16S_2021_processed.ps") # import ps

# inspect
# view(WI_16S_2021_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(WI_16S_2021_processed.ps)$month <- as.numeric(sample_data(WI_16S_2021_processed.ps)$month) 
sample_data(WI_16S_2021_processed.ps)$season <- "Fall" # add season as character

# add block
sample_data(WI_16S_2021_processed.ps)$block <- as.numeric(substring(sample_data(WI_16S_2021_processed.ps)$plot, 1,1)) 

# add in the treatment numbers and metadata
WI_16S_2021_m.ps <- add_treatment_info(WI_16S_2021_processed.ps, "WI") # m signifies metadata-added

# final check and write out
view(WI_16S_2021_m.ps@sam_data)
saveRDS(WI_16S_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/WI_16S_2021_m.ps")


# 2022 
WI_16S_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/WI_16S_2022_processed.ps")

# take a look at the metadata
view(WI_16S_2022_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(WI_16S_2022_processed.ps)$month <- as.numeric(sample_data(WI_16S_2022_processed.ps)$month) 
sample_data(WI_16S_2022_processed.ps)$season <- ifelse(sample_data(WI_16S_2022_processed.ps)$month == 4, "Spring",
                                                       ifelse(sample_data(WI_16S_2022_processed.ps)$month == 6, "Summer", NA))

# add block
sample_data(WI_16S_2022_processed.ps)$block <- as.numeric(substring(sample_data(WI_16S_2022_processed.ps)$plot, 1,1))

# add in the treatment numbers and metadata
WI_16S_2022_m.ps <- add_treatment_info(WI_16S_2022_processed.ps, "WI") # m signifies metadata-added

# final check aWI write out
view(sample_data(WI_16S_2022_m.ps))
saveRDS(WI_16S_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/WI_16S_2022_m.ps")
```

ITS.

``` r
# 2021
WI_ITS_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/WI_ITS_2021_processed.ps") # import ps

# inspect
# view(WI_ITS_2021_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(WI_ITS_2021_processed.ps)$month <- as.numeric(sample_data(WI_ITS_2021_processed.ps)$month) 
sample_data(WI_ITS_2021_processed.ps)$season <- "Fall" # add season as character

# add block
sample_data(WI_ITS_2021_processed.ps)$block <- as.numeric(substring(sample_data(WI_ITS_2021_processed.ps)$plot, 1,1)) 

# add in the treatment numbers and metadata
WI_ITS_2021_m.ps <- add_treatment_info(WI_ITS_2021_processed.ps, "WI") # m signifies metadata-added

# fix taxonomy labels
WI_ITS_2021_m.ps <- fix_its_taxa(WI_ITS_2021_m.ps)

# final check and write out
view(WI_ITS_2021_m.ps@sam_data)
saveRDS(WI_ITS_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/WI_ITS_2021_m.ps")


# 2022 
WI_ITS_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/WI_ITS_2022_processed.ps")

# take a look at the metadata
# view(WI_ITS_2022_processed.ps@sam_data)

# convert month to numeric and add season
sample_data(WI_ITS_2022_processed.ps)$month <- as.numeric(sample_data(WI_ITS_2022_processed.ps)$month) 
sample_data(WI_ITS_2022_processed.ps)$season <- ifelse(sample_data(WI_ITS_2022_processed.ps)$month == 4, "Spring",
                                                       ifelse(sample_data(WI_ITS_2022_processed.ps)$month == 6, "Summer", NA))

# add block
sample_data(WI_ITS_2022_processed.ps)$block <- as.numeric(substring(sample_data(WI_ITS_2022_processed.ps)$plot, 1,1)) # Block is first number of plot #

# add in the treatment numbers and metadata
WI_ITS_2022_m.ps <- add_treatment_info(WI_ITS_2022_processed.ps, "WI") # m signifies metadata-added

# fix taxonomy labels
WI_ITS_2022_m.ps <- fix_its_taxa(WI_ITS_2022_m.ps)

# final check and write out
view(sample_data(WI_ITS_2022_m.ps))
saveRDS(WI_ITS_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/WI_ITS_2022_m.ps")
```

Uploaded phyloseqs to Box, Matt R. and Rick L. notified 12-16.

#### USDA: Processed 12-16-23.

16S.

``` r
# 2021
US_16S_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/larkin_16S_2021_processed.ps") # import ps

# inspect
# sample_data(US_16S_2021_processed.ps)

# convert month to numeric and add season
sample_data(US_16S_2021_processed.ps)$month <- as.numeric(sample_data(US_16S_2021_processed.ps)$month) 
sample_data(US_16S_2021_processed.ps)$season <- ifelse(sample_data(US_16S_2021_processed.ps)$month == 5, "Spring",
                                                       ifelse(sample_data(US_16S_2021_processed.ps)$month == 7, "Summer",
                                                              ifelse(sample_data(US_16S_2021_processed.ps)$month == 9, "Fall", NA))) 
# import everything including block
see <- left_join(data.frame(sample_data(US_16S_2021_processed.ps)), 
                 larkin.decoder %>% select(-state), 
                 by="plot") # merges the dataframes by plot, adding treatment info as well

# rewrite the rotation in column two (rotation.x) with the latter one (rotation.y) and correct column names
see$rotation.x <- see$rotation.y
see <- see %>% select(-rotation.y)
colnames(see)[3] <- "rotation"
colnames(see)[15] <- "general_category"

rownames(see) <- sample_names(US_16S_2021_processed.ps) # add sample names to the dataframe
all.equal(colnames(see), colnames(sample_data(MN_16S_2021_m.ps))) # check to see the column names are the same as other ps objects
sample_data(US_16S_2021_processed.ps) <- see # write the new dataframe into the ps object
US_16S_2021_m.ps <- US_16S_2021_processed.ps # new ps object

# final check and write out
# sample_data(US_16S_2021_m.ps)
saveRDS(US_16S_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/US_16S_2021_m.ps")


# 2022 
US_16S_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/16S_processed/larkin_16S_2022_processed.ps") # import ps

# inspect
# sample_data(US_16S_2022_processed.ps)

# convert month to numeric and add season
sample_data(US_16S_2022_processed.ps)$month <- str_sub(sample_data(US_16S_2022_processed.ps)$month, 1, 1)
sample_data(US_16S_2022_processed.ps)$month <- as.numeric(sample_data(US_16S_2022_processed.ps)$month) 
sample_data(US_16S_2022_processed.ps)$season <- ifelse(sample_data(US_16S_2022_processed.ps)$month == 5, "Spring",
                                                ifelse(sample_data(US_16S_2022_processed.ps)$month == 7, "Summer", NA))

# import everything including block
see <- left_join(data.frame(sample_data(US_16S_2022_processed.ps)), 
                 larkin.decoder %>% select(-state), 
                 by="plot") # merges the dataframes by plot, adding treatment info as well

# rewrite the rotation in column two (rotation.x) with the latter one (rotation.y) and correct column names
see$rotation.x <- see$rotation.y
see <- see %>% select(-rotation.y)
colnames(see)[3] <- "rotation"
colnames(see)[15] <- "general_category"

rownames(see) <- sample_names(US_16S_2022_processed.ps) # add sample names to the dataframe
all.equal(colnames(see), colnames(sample_data(MN_16S_2022_m.ps))) # check to see the column names are the same as other ps objects
sample_data(US_16S_2022_processed.ps) <- see # write the new dataframe into the ps object
US_16S_2022_m.ps <- US_16S_2022_processed.ps # new ps object

# final check and write out
# sample_data(US_16S_2022_m.ps)
saveRDS(US_16S_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/16S_metadata_added/US_16S_2022_m.ps")
```

ITS.

``` r
# 2021
US_ITS_2021_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/larkin_ITS_2021_processed.ps") # import ps

# inspect
# sample_data(US_ITS_2021_processed.ps)

# convert month to numeric and add season
sample_data(US_ITS_2021_processed.ps)$month <- as.numeric(sample_data(US_ITS_2021_processed.ps)$month) 
sample_data(US_ITS_2021_processed.ps)$season <- ifelse(sample_data(US_ITS_2021_processed.ps)$month == 5, "Spring",
                                                       ifelse(sample_data(US_ITS_2021_processed.ps)$month == 7, "Summer",
                                                              ifelse(sample_data(US_ITS_2021_processed.ps)$month == 9, "Fall", NA))) 
# import everything including block
see <- left_join(data.frame(sample_data(US_ITS_2021_processed.ps)), 
                 larkin.decoder %>% select(-state), 
                 by="plot") # merges the dataframes by plot, adding treatment info as well

# rewrite the rotation in column two (rotation.x) with the latter one (rotation.y) and correct column names
see$rotation.x <- see$rotation.y
see <- see %>% select(-rotation.y)
colnames(see)[3] <- "rotation"
colnames(see)[15] <- "general_category"

rownames(see) <- sample_names(US_ITS_2021_processed.ps) # add sample names to the dataframe
all.equal(colnames(see), colnames(sample_data(WI_ITS_2021_m.ps))) # check to see the column names are the same as other ps objects
sample_data(US_ITS_2021_processed.ps) <- see # write the new dataframe into the ps object
US_ITS_2021_m.ps <- US_ITS_2021_processed.ps # new ps object

# fix taxonomy labels
US_ITS_2021_m.ps <- fix_its_taxa(US_ITS_2021_m.ps)

# final check and write out
# sample_data(US_ITS_2021_m.ps)
saveRDS(US_ITS_2021_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/US_ITS_2021_m.ps")


# 2022
US_ITS_2022_processed.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/ITS_processed/larkin_ITS_2022_processed.ps") # import ps

# inspect
# sample_data(US_ITS_2022_processed.ps)

# convert month to numeric and add season
sample_data(US_ITS_2022_processed.ps)$month <- as.numeric(str_sub(rownames(sample_data(US_ITS_2022_processed.ps)), 15, 15)) # anomalous fix
sample_data(US_ITS_2022_processed.ps)$season <- ifelse(sample_data(US_ITS_2022_processed.ps)$month == 5, "Spring",
                                                       ifelse(sample_data(US_ITS_2022_processed.ps)$month == 7, "Summer", NA)) 
# import everything including block
see <- left_join(data.frame(sample_data(US_ITS_2022_processed.ps)), 
                 larkin.decoder %>% select(-state), 
                 by="plot") # merges the dataframes by plot, adding treatment info as well

# rewrite the rotation in column two (rotation.x) with the latter one (rotation.y) and correct column names
see$rotation.x <- see$rotation.y
see <- see %>% select(-rotation.y)
colnames(see)[3] <- "rotation"
colnames(see)[15] <- "general_category"

rownames(see) <- sample_names(US_ITS_2022_processed.ps) # add sample names to the dataframe
all.equal(colnames(see), colnames(sample_data(WI_ITS_2022_m.ps))) # check to see the column names are the same as other ps objects
sample_data(US_ITS_2022_processed.ps) <- see # write the new dataframe into the ps object
US_ITS_2022_m.ps <- US_ITS_2022_processed.ps # new ps object

# fix taxonomy labels
US_ITS_2022_m.ps <- fix_its_taxa(US_ITS_2022_m.ps)

# final check and write out
# sample_data(US_ITS_2022_m.ps)
saveRDS(US_ITS_2022_m.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/ITS_metadata_added/US_ITS_2022_m.ps")
```

Uploaded phyloseqs to Box, Bob notified 12-16.

### Additional notes, conclusions:

In conclusion, I added all the metadata corresponding to treatment to
all the phyloseq objects. Could I have written another function to add
and curate sample data and apply it to all input phyloseqs so that I
didn’t end up with nearly a thousand lines of code? Yeah, probably,
but 1) there were a few anomalies for certain sites, and 2) more
importantly, I wanted to make sure I laid eyes on all of these to verify
accuracy. I checked at least one plot number for each input phyloseq
object to verify that it corresponded to the correct treatment
description.

### Session info

``` r
sessionInfo()
```

    ## R version 4.2.3 (2023-03-15)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] phyloseq_1.42.0 lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0  
    ##  [5] dplyr_1.1.1     purrr_1.0.1     readr_2.1.4     tidyr_1.3.0    
    ##  [9] tibble_3.2.1    ggplot2_3.4.2   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.10            ape_5.7-1              lattice_0.20-45       
    ##  [4] Biostrings_2.66.0      digest_0.6.31          foreach_1.5.2         
    ##  [7] utf8_1.2.3             R6_2.5.1               GenomeInfoDb_1.34.9   
    ## [10] plyr_1.8.8             stats4_4.2.3           evaluate_0.20         
    ## [13] pillar_1.9.0           zlibbioc_1.44.0        rlang_1.1.0           
    ## [16] data.table_1.14.8      rstudioapi_0.14        vegan_2.6-4           
    ## [19] S4Vectors_0.36.2       Matrix_1.5-4.1         rmarkdown_2.21        
    ## [22] splines_4.2.3          igraph_1.4.1           RCurl_1.98-1.12       
    ## [25] munsell_0.5.0          compiler_4.2.3         xfun_0.38             
    ## [28] pkgconfig_2.0.3        BiocGenerics_0.44.0    multtest_2.54.0       
    ## [31] mgcv_1.8-42            htmltools_0.5.5        biomformat_1.26.0     
    ## [34] tidyselect_1.2.0       GenomeInfoDbData_1.2.9 IRanges_2.32.0        
    ## [37] codetools_0.2-19       permute_0.9-7          fansi_1.0.4           
    ## [40] crayon_1.5.2           tzdb_0.4.0             withr_2.5.0           
    ## [43] rhdf5filters_1.10.1    MASS_7.3-58.2          bitops_1.0-7          
    ## [46] grid_4.2.3             nlme_3.1-162           jsonlite_1.8.4        
    ## [49] gtable_0.3.3           lifecycle_1.0.3        magrittr_2.0.3        
    ## [52] scales_1.2.1           cli_3.6.1              stringi_1.7.12        
    ## [55] XVector_0.38.0         reshape2_1.4.4         generics_0.1.3        
    ## [58] vctrs_0.6.1            Rhdf5lib_1.20.0        iterators_1.0.14      
    ## [61] tools_4.2.3            ade4_1.7-22            Biobase_2.58.0        
    ## [64] glue_1.6.2             hms_1.1.3              survival_3.5-3        
    ## [67] parallel_4.2.3         fastmap_1.1.1          yaml_2.3.7            
    ## [70] rhdf5_2.42.0           timechange_0.2.0       colorspace_2.1-0      
    ## [73] cluster_2.1.4          knitr_1.42
