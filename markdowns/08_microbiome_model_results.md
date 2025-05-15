model_results
================
Scott Klasek
2024-09-10

## Purpose

Summarize model information here because the document 45_SEM_stuff is
too big and scary.

## Setup

#### load libraries

``` r
packages <- c("tidyverse", "phyloseq", "speedyseq", "patchwork", "compositions", "piecewiseSEM", "lme4")
invisible(lapply(packages, require, character.only = TRUE))
```

    ## Loading required package: tidyverse

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
    ## Loading required package: phyloseq
    ## 
    ## Loading required package: speedyseq
    ## 
    ## 
    ## Attaching package: 'speedyseq'
    ## 
    ## 
    ## The following objects are masked from 'package:phyloseq':
    ## 
    ##     filter_taxa, plot_bar, plot_heatmap, plot_tree, psmelt, tax_glom,
    ##     tip_glom, transform_sample_counts
    ## 
    ## 
    ## Loading required package: patchwork
    ## 
    ## Loading required package: compositions
    ## 
    ## Welcome to compositions, a package for compositional data analysis.
    ## Find an intro with "? compositions"
    ## 
    ## 
    ## 
    ## Attaching package: 'compositions'
    ## 
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     anova, cor, cov, dist, var
    ## 
    ## 
    ## The following object is masked from 'package:graphics':
    ## 
    ##     segments
    ## 
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     %*%, norm, scale, scale.default
    ## 
    ## 
    ## Loading required package: piecewiseSEM
    ## 
    ## 
    ##   This is piecewiseSEM version 2.3.0.1.
    ## 
    ## 
    ##   Questions or bugs can be addressed to <jslefche@gmail.com>.
    ## 
    ## Loading required package: lme4
    ## 
    ## Loading required package: Matrix
    ## 
    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

#### load data

``` r
# lists of phyloseq objects
# these include metadata updated in doc 39.Rmd
# they are ALL years, ASVs and samples NOT subsetted or transformed
its.ps.list <- readRDS(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/its.ps.list")
its.ps.list <- its.ps.list[c(1:7,9)] # remove Larkin samples
bact.ps.list <- readRDS(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/bact.ps.list")
bact.ps.list <- bact.ps.list[c(1:7,9)] # remove Larkin samples

# jim's spreadsheet of soil chemical metadata
jim.info <- read.csv(file="/Users/klas0061/Desktop/UMN/jim_info/PSHP_ALL_obj_1_all_data_2024_02_29.csv") # import
```

#### reformat soil chemical metadata so it can be added to phyloseq objects using merge_with_jim_data

``` r
# basic clean-up stuff
colnames(jim.info) <- jim.info[1,] # fix column names (there are two header cols)
jim.info <- jim.info[2:nrow(jim.info),] # omit redundant column name 

for (i in c(2:3,5:6)){jim.info[,i] <- as.numeric(jim.info[,i])} # convert columns that should be numeric into numeric

# remove extraneous year columns that correspond to yield data only
jim.info <- jim.info[,c(1:101,103:129,131:ncol(jim.info))]

# select data columns to merge by, and to add (an example, can always select more later)
jim.info.s <- jim.info %>% dplyr::select(State, Objective, Rotation, Plot, Year, Month, 
                                  pH, `Nitrate-N (ppm)`, `NH4-N (ppm)`, `P-Olsen (ppm)`, `P-Bray (ppm)`,
                                  `OM (%)`, `Solvita (ppm)`, `Total C (%)`, `Total organic C (%)`, `POX-C (ppm)`) 

for (i in 7:16){jim.info.s[,i] <- as.numeric(jim.info.s[,i])} # again, convert columns that should be numeric into numeric
```

    ## Warning: NAs introduced by coercion
    ## Warning: NAs introduced by coercion
    ## Warning: NAs introduced by coercion
    ## Warning: NAs introduced by coercion
    ## Warning: NAs introduced by coercion
    ## Warning: NAs introduced by coercion
    ## Warning: NAs introduced by coercion
    ## Warning: NAs introduced by coercion
    ## Warning: NAs introduced by coercion
    ## Warning: NAs introduced by coercion

``` r
# change OR Spring 2020 month to 4 so Cultivar (whether field is in potato or not) is not NA
jim.info.s[which(jim.info.s$State=="OR" & jim.info.s$Year==20 & jim.info.s$Month=="3"),"Month"] <- "4" 

jim.info.s$Month <- as.numeric(jim.info.s$Month) # more character to numeric
```

#### define functions

``` r
# drop_ghost_asvs drops ASVs that are not present in any samples- but are still kept because you've 
# subsetted samples within a larger phyloseq object. This returns the phyloseq object with the subsetted count table.
drop_ghost_asvs <- function(ps){
  
  # extract the count table from the phyloseq, with ASVs as columns
    if(taxa_are_rows(ps)){counts <- t(ps@otu_table)} else{ # write count table with ASVs as rows
    counts <- ps@otu_table
    }
  # remove empty rows corresponding to ASVs that are present in 0 samples
  counts <- counts[,colSums(counts)>0]
  
  # write the count table back in
  otu_table(ps) <- counts # drops all zero-count ASVs from the tax table and refseq as well, in contrast to 'ps@otu_table <- counts'
  return(ps)
}

# subset_occupancy drops ASVs detected in less than a certain proportion of the samples
# all it takes is a phyloseq object. it drops the ASVs below the occupancy cutoff, and moves their counts to 
# a "summing" column. Tax table now has a summing column as well. Output is a phyloseq object. 
subset_occupancy <- function(ps, occ_cutoff){
  
  # first, calculate occupancy of each ASV (proportion of samples detected in, from 0 to 1)
  occ <- vector("numeric")
    for (i in 1:ncol(ps@otu_table)) {
      occ[i] <- sum(ps@otu_table[,i] != 0)/nrow(ps@otu_table[,i])
    }
  # write ASV names
  names(occ) <- colnames(ps@otu_table)
    
  # select ASVs to keep based on occupancy 
  keep.asvs <- occ[which(occ >= (occ_cutoff))] # keep only the ASVs above the cutoff
  
  # print some useful information
  cat((length(keep.asvs)), "features are kept. \n")
  
  # subset the count table to maintain compositionality
  keep.table <- ps@otu_table[,names(keep.asvs)] # the kept ASVs, which does not include a summing column
  drop.table <- ps@otu_table[,setdiff(colnames(ps@otu_table), names(keep.asvs))] # the dropped ASVs
  summing <- rowSums(drop.table) # add up counts for dropped features to make a summing column
  keep.new <- cbind(keep.table, summing) # add new summing to the keep table
  
  # make a subsetted tax_table
  taxa <- tax_table(ps)[names(keep.asvs),] # with only the kept ASVs over the occupancy cutoff
  summing <- c("Summing", rep(NA, times = ncol(taxa)-1))
  taxa2 <- rbind(taxa, summing)
    
  # put the new count table and tax table into a new ps object
  # we must say goodbye to our refseq here because "summing" doesn't have a sequence associated with it.
  new.ps <- phyloseq(sample_data(ps),
                     otu_table(keep.new, taxa_are_rows = FALSE),
                     tax_table(taxa2))
  return(new.ps)
}

# transform_clr takes a phyloseq object and transforms the count table using a centered-log-ratio
# it spits out a phyloseq object with transformed counts.
transform_clr <- function(ps){
  counts.clr <- as.matrix(clr(ps@otu_table))
  ps.new <- phyloseq(otu_table(counts.clr, taxa_are_rows = FALSE),
                     tax_table(ps),
                     sample_data(ps))
  cat("Counts transformed with CLR. \n")
  return(ps.new)
}

# merge_with_jim_data takes as input a phyloseq object and a dataframe of objective 1 data that Jim has curated, which can be subsetted 
# it writes data from Jim's spreadsheet into the corresponding phyloseq object, based on combinations
# make sure jim's data is cleaned up according to the steps above!
merge_with_jim_data <- function(ps, jim.info.subset){
  df <- left_join(data.frame(ps@sam_data), jim.info.subset, 
                by = c("state" = "State", "objective" = "Objective", "rotation" = "Rotation",
                       "plot" = "Plot", "year" = "Year", "month" = "Month")) # merges 
  #the unique combo of columns from phyloseq sample data and jim's data
  rownames(df) <- rownames(ps@sam_data) # ok as long as the rows are in the same order, which they are
  
  # put df into the sample_data.
  sample_data(ps) <- df # this also automatically drops the ghost ASVs and refsequences from the ps object
  return(ps)
}
```

## Microbiome-informed models

### Scenario 1: OR 2-yr bacterial consortia in response to amendment

From the 2-yr rotation only

``` r
# begin by subsetting the phyloseq object by site, rotation, and year
ps6 <- bact.ps.list[["OR.16S.ps"]] %>% 
  subset_samples(year == 22 & rotation == 2 & !is.na(Total.yield) & general_category != "Fumigated") 

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples 
sub.df <- jim.info.s %>% filter(State == "OR" & Year == 22)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps6 <- merge_with_jim_data(ps6, sub.df)

# define ASVs to model
or.amd.bact.summer.asvs <- c("ASV1254","ASV21630", "ASV2399",
                          "ASV3632",  "ASV4948", "ASV5076",  "ASV7919", "ASV792")

# keep only the summer samples
ps6 <- ps6 %>% subset_samples(season == "Summer") %>% 
  drop_ghost_asvs() %>%
  subset_occupancy(0.5) %>% 
  transform_clr() 
```

    ## 5246 features are kept. 
    ## Counts transformed with CLR.

``` r
# add in the abundance of the 2 yr rotation SUMMER associated ASVs to the sample data 
df1 <- cbind(data.frame(ps6@sam_data), ps6@otu_table[,or.amd.bact.summer.asvs])

# make a consortium column of the sum of the ASVs
df1$consortium <- rowSums(df1[,startsWith(colnames(df1), "ASV")])

### tidying up the dataframe 
# remove NA cols
df1 <- df1 %>% dplyr::select(-cover_2020, -green_cov_pct_pre_kill, -Total.C...., -Total.organic.C....)

# there was one outlier OM value at 2.8, next highest OM recorded in Oregon any time was 1.8
# also keep only Norkotah cultivar
df1.no <- df1 %>% filter(`OM....` < 2 & cultivar == "Norkotah")

# lm where amendment directly impacts yield as well (this cannot be modeled in SEM)
lm4 <- lm(Total.yield ~ Amended + consortium, df1.no)
summary(lm4)
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ Amended + consortium, data = df1.no)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -64.277 -19.503  -8.684  29.540  66.776 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   531.14      22.26  23.860 1.01e-08 ***
    ## AmendedTRUE  -224.50     118.84  -1.889   0.0956 .  
    ## consortium     32.02      11.29   2.836   0.0220 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 46.47 on 8 degrees of freedom
    ## Multiple R-squared:  0.7281, Adjusted R-squared:  0.6601 
    ## F-statistic: 10.71 on 2 and 8 DF,  p-value: 0.005467

``` r
extractAIC(lm4) # gives AIC of 87, which is lower than the PSEM
```

    ## [1]  3.00000 86.95183

``` r
# calculate 
est_amd <- summary(lm(consortium ~ Amended, df1.no))$coefficients[2, 1]
est_cons <- summary(lm(Total.yield ~ consortium, df1.no))$coefficients[2, 1]

# Yield increase as an effect of treatment on microbiome (cwt/ac) is the product of the individual coefficients 
yield_inc <- est_amd * est_cons

# percent increase in total yield from treatment-related microbiome changes, relative to control  
pct_yield_inc <- est_amd * est_cons / mean(df1.no %>% filter(general_category == "Control") %>% pull(Total.yield)) * 100

# best model
summary(lm(Total.yield ~ consortium + Amended, df1.no))
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ consortium + Amended, data = df1.no)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -64.277 -19.503  -8.684  29.540  66.776 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   531.14      22.26  23.860 1.01e-08 ***
    ## consortium     32.02      11.29   2.836   0.0220 *  
    ## AmendedTRUE  -224.50     118.84  -1.889   0.0956 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 46.47 on 8 degrees of freedom
    ## Multiple R-squared:  0.7281, Adjusted R-squared:  0.6601 
    ## F-statistic: 10.71 on 2 and 8 DF,  p-value: 0.005467

``` r
# extract p- and R2 values
overall_p <- function(my_model) {
    f <- summary(my_model)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

pval <- overall_p(lm(Total.yield ~ consortium + Amended, df1.no))
rsq <- summary(lm(Total.yield ~ consortium + Amended, df1.no))[[8]]
```

### Scenario 2: ID 3-yr bacterial consortia in response to amendment

From the 3-yr rotation only

``` r
# begin by subsetting the phyloseq object by site, rotation, and year
ps7 <- bact.ps.list[["ID.16S.ps"]] %>% 
  subset_samples(year == 22 & rotation == 3 & !is.na(Total.yield) & Fumigated == FALSE) 

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples 
sub.df <- jim.info.s %>% filter(State == "ID" & Year == 22 & Rotation == 3)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps7 <- merge_with_jim_data(ps7, sub.df)

# define ASVs
id.amd.bact.summer.asvs <- c("ASV3781", "ASV12474", "ASV17047", "ASV2340",  "ASV2345",  
                             "ASV2589", "ASV4416",  "ASV4825",  "ASV4830")

# keep only the summer samples, remove any below 50% occupancy, and transform abundances
ps7 <- ps7 %>% subset_samples(season == "Summer") %>% 
  drop_ghost_asvs() %>%
  subset_occupancy(0.5) %>% 
  transform_clr() 
```

    ## 5180 features are kept. 
    ## Counts transformed with CLR.

``` r
# add in the abundance of the summer associated ASVs to the sample data 
df2 <- cbind(data.frame(ps7@sam_data), ps7@otu_table[,id.amd.bact.summer.asvs])

# make a consortium column of the sum of the ASVs
df2$consortium <- rowSums(df2[,startsWith(colnames(df2), "ASV")])

# remember that Green Manure treatments are now classified as controls
df2[which(df2$amendment == "Green Manure"),"general_category"] <- "Control"

# control for cultivar: keep only Norkotahs
df2 <- df2 %>% filter(cultivar == "Norkotah")

# drop NA cols
df2 <- df2[, colSums(is.na(df2)) != nrow(df2)]

# convert amended treatment to 0/1 (NECESSARY FOR EXTRACTING MODEL INFO FROM PSEMs)
df2$Amended <- as.integer(as.logical(df2$Amended))

# SEM 2: amendment increases OM, OM and amendment increase consortium, and consortium increases total yields
psem2 <- psem(lm(`OM....` ~ Amended, df2),
               lm(consortium ~ `OM....` + Amended, df2),
               lm(Total.yield ~ consortium, df2))  

fisherC(psem2) # valid
```

    ##   Fisher.C df P.Value
    ## 1    5.122  4   0.275

``` r
coefs(psem2)
```

    ##      Response  Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ## 1      OM....    Amended   0.1125    0.0409 14     2.7495  0.0157       0.5922
    ## 2  consortium     OM....  -4.0992    6.0953 13    -0.6725  0.5130      -0.0649
    ## 3  consortium    Amended  11.9669    1.1580 13    10.3340  0.0000       0.9969
    ## 4 Total.yield consortium   4.5492    0.9639 14     4.7195  0.0003       0.7836
    ##      
    ## 1   *
    ## 2    
    ## 3 ***
    ## 4 ***

``` r
# extract p and R2 values 
psemrsq <- rsquared(psem2)[which(rsquared(psem2)$Response == "Total.yield"),"R.squared"]
rsq <- append(rsq, psemrsq)
pval <- append(pval, fisherC(psem2, conserve = T)[1,3])

# calculate 
est_amd <- summary(lm(consortium ~ Amended, df2))$coefficients[2, 1]
est_cons <- summary(lm(Total.yield ~ consortium, df2))$coefficients[2, 1]

# Yield increase as an effect of treatment on microbiome (cwt/ac) is the product of the individual coefficients 
yield_inc <- append(yield_inc, (est_amd * est_cons))

# percent increase in total yield from treatment-related microbiome changes, relative to control  
pct <- est_amd * est_cons / mean(df2 %>% filter(general_category == "Control") %>% pull(Total.yield)) * 100
pct_yield_inc <- append(pct_yield_inc, pct)
```

### Scenario 3: ME 2-yr ITS in response to amendment

The 3-yr rotation is the one that had the yield increase, but the 2-yr
one had the two amendment-influenced, yield-positive fungal ASVs.

``` r
# begin by subsetting the phyloseq object by site, rotation, and year
ps8 <- its.ps.list[["ME.ITS.ps"]] %>% 
  subset_samples(year == 22 & rotation == 2 & !is.na(Total.yield) & general_category != "Fumigated") 

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples
sub.df <- jim.info.s %>% filter(State == "ME" & Year == 22 & Rotation == 2)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps8 <- merge_with_jim_data(ps8, sub.df)

# define ASVs
me.amd.its.spring.asvs <- c("ASV501", "ASV1906")

# keep only the spring samples, drop, and transform
ps8 <- ps8 %>% subset_samples(season == "Spring") %>% 
  drop_ghost_asvs() %>%
  subset_occupancy(0.5) %>% 
  transform_clr() 
```

    ## 803 features are kept. 
    ## Counts transformed with CLR.

``` r
# add in the abundance of the spring associated ASVs to the sample data 
df3 <- cbind(data.frame(ps8@sam_data), ps8@otu_table[,me.amd.its.spring.asvs])

# make a consortium column of the sum of the ASVs
df3$consortium <- rowSums(df3[,startsWith(colnames(df3), "ASV")])

# since we have two cultivars of nearly equal abundance, control for cultivar as a factor
df3$cultivar <- ifelse(df3$cultivar == "Burbank", 0, 
                       ifelse(df3$cultivar == "Caribou", 1, NA))

# model
summary(lm(Total.yield ~ Amended + consortium + cultivar, df3))
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ Amended + consortium + cultivar, data = df3)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -57.365 -13.953   2.008  12.530  49.537 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  165.251     10.536  15.684 4.54e-13 ***
    ## AmendedTRUE   -8.876     17.752  -0.500   0.6223    
    ## consortium     9.501      4.426   2.147   0.0436 *  
    ## cultivar     -15.018     15.276  -0.983   0.3367    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 26.18 on 21 degrees of freedom
    ## Multiple R-squared:  0.5301, Adjusted R-squared:  0.4629 
    ## F-statistic: 7.896 on 3 and 21 DF,  p-value: 0.001029

``` r
pval <- append(pval, overall_p(lm(Total.yield ~ Amended + consortium + cultivar, df3)))
rsq <- append(rsq, summary(lm(Total.yield ~ Amended + consortium + cultivar, df3))[[8]])

# calculate coefficients
# cultivar left out of consortium -> yield model because insignificant
est_amd <- summary(lm(consortium ~ Amended + cultivar, df3))$coefficients[2, 1]
est_cons <- summary(lm(Total.yield ~ consortium, df3))$coefficients[2, 1]

# Yield increase as an effect of treatment on microbiome (cwt/ac) is the product of the individual coefficients 
yield_inc <- append(yield_inc, (est_amd * est_cons))

# percent increase in total yield from treatment-related microbiome changes, relative to control  
pct <- est_amd * est_cons / mean(df3 %>% filter(general_category == "Control") %>% pull(Total.yield)) * 100
pct_yield_inc <- append(pct_yield_inc, pct)
```

### Scenario 4: MI 3-yr ITS in response to amendment

``` r
# begin by subsetting the phyloseq object by site, rotation, and year
ps9 <- its.ps.list[["MI.ITS.ps"]] %>% 
  subset_samples(year == 22 & rotation == 3 & !is.na(Total.yield) & general_category != "Control") 

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples
sub.df <- jim.info.s %>% filter(State == "MI" & Year == 22 & Rotation == 3)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps9 <- merge_with_jim_data(ps9, sub.df)

# define ASVs to model
mi.amd.its.summer.asvs <- c("ASV1099", "ASV1194")

# keep only the spring samples, drop, and transform
ps9 <- ps9 %>% subset_samples(season == "Summer") %>% 
  drop_ghost_asvs() %>%
  subset_occupancy(0.5) %>% 
  transform_clr() 
```

    ## 793 features are kept. 
    ## Counts transformed with CLR.

``` r
# add in the abundance of the spring associated ASVs to the sample data 
df4 <- cbind(data.frame(ps9@sam_data), ps9@otu_table[,mi.amd.its.summer.asvs])

# make a consortium column of the sum of the ASVs
df4$consortium <- rowSums(df4[,startsWith(colnames(df4), "ASV")])

# convert amended treatment to 0/1
df4$Amended <- as.integer(as.logical(df4$Amended))

### tidying up the dataframe 
# remove NA cols
df4 <- df4 %>% dplyr::select(-cover_2019, -cover_2020, -green_cov_pct_pre_kill, -Total.C...., 
                             -Total.organic.C...., -VPPG, -Root.Lesion.100.cc, -rotation_comparison)

# encode cultivar as a dummy variable
df4$cultivar <- ifelse(df4$cultivar == "Superior", 0,
                       ifelse(df4$cultivar == "Burbank", 1, 2))

# SEM: amendment directly stimulates consortium, which increases yields
psem1 <- psem(lm(consortium ~ Amended + cultivar, df4),
              lm(Total.yield ~ consortium + cultivar, df4)) 

fisherC(psem1) # valid
```

    ##   Fisher.C df P.Value
    ## 1    0.156  2   0.925

``` r
coefs(psem1)
```

    ##      Response  Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ## 1  consortium    Amended   6.7092    0.7709 20     8.7034  0.0000       1.1131
    ## 2  consortium   cultivar   1.5034    0.4658 20     3.2273  0.0042       0.4128
    ## 3 Total.yield consortium   9.4768    3.4580 20     2.7405  0.0126       0.5498
    ## 4 Total.yield   cultivar  12.0465   12.5952 20     0.9564  0.3503       0.1919
    ##      
    ## 1 ***
    ## 2  **
    ## 3   *
    ## 4

``` r
# extract p and R2 values
psemrsq <- rsquared(psem1)[which(rsquared(psem1)$Response == "Total.yield"),"R.squared"]
rsq <- append(rsq, psemrsq)
pval <- append(pval, fisherC(psem1, conserve = T)[1,3])

# calculate coefficients
# cultivar left out of consortium -> yield model because insignificant
est_amd <- summary(lm(consortium ~ Amended + cultivar, df4))$coefficients[2, 1]
est_cons <- summary(lm(Total.yield ~ consortium, df4))$coefficients[2, 1]

# percent increase in total yield from treatment-related microbiome changes, relative to control  
pct <- est_amd * est_cons / mean(df4 %>% filter(general_category == "Fumigated") %>% pull(Total.yield)) * 100

# Yield increase as an effect of treatment on microbiome (cwt/ac) is the product of the individual coefficients 
yield_inc <- append(yield_inc, (est_amd * est_cons))

# percent increase in total yield from treatment-related microbiome changes, relative to control  
pct_yield_inc <- append(pct_yield_inc, pct)
```

### Scenario 5: OR 3-yr ITS in response to fumigation

``` r
# begin by subsetting the phyloseq object by site, rotation, and year
ps1 <- its.ps.list[["OR.ITS.ps"]] %>% 
  subset_samples(year == 22 & rotation == 3 & !is.na(Total.yield) & general_category != "Amended") 

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples 
sub.df <- jim.info.s %>% filter(State == "OR" & Year == 22)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps1 <- merge_with_jim_data(ps1, sub.df)

# keep only the spring samples, drop, and transform
ps1 <- ps1 %>% subset_samples(season == "Summer") %>% 
  drop_ghost_asvs() %>%
  subset_occupancy(0.5) %>% 
  transform_clr() 
```

    ## 703 features are kept. 
    ## Counts transformed with CLR.

``` r
# only one eukaryotic ASV associated with fumigation/yield, and it is 60 DAP
df5 <- cbind(data.frame(ps1@sam_data), ps1@otu_table[,"ASV2495"])

# encode cultivar as a dummy variable
df5$cultivar <- ifelse(df5$cultivar == "Burbank", 0,
                       ifelse(df5$cultivar == "Norkotah", 1, 2))

# Total.yield ~ 60 DAP abundance of ASV2495 + Fumigated + cultivar. Cultivars coded as dummy variables.
summary(lm(Total.yield ~ ASV2495 + Fumigated + cultivar, df5))
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ ASV2495 + Fumigated + cultivar, data = df5)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -93.80 -50.03  18.04  44.05 103.17 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     527.93      39.23  13.456 3.55e-08 ***
    ## ASV2495          74.12      25.39   2.919   0.0140 *  
    ## FumigatedTRUE   -88.03      65.19  -1.350   0.2040    
    ## cultivar        115.45      42.40   2.723   0.0198 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 68.56 on 11 degrees of freedom
    ## Multiple R-squared:  0.6088, Adjusted R-squared:  0.5022 
    ## F-statistic: 5.707 on 3 and 11 DF,  p-value: 0.01321

``` r
pval <- append(pval, overall_p(lm(Total.yield ~ ASV2495 + Fumigated + cultivar, df5)))
rsq <- append(rsq, summary(lm(Total.yield ~ ASV2495 + Fumigated + cultivar, df5))[[8]])

# calculate coefficients
# cultivar left out of consortium -> yield model because insignificant
est_amd <- summary(lm(ASV2495 ~ Amended + cultivar, df5))$coefficients[2, 1]
est_cons <- summary(lm(Total.yield ~ ASV2495, df5))$coefficients[2, 1]

# Yield increase as an effect of treatment on microbiome (cwt/ac) is the product of the individual coefficients 
yield_inc <- append(yield_inc, (est_amd * est_cons))

# percent increase in total yield from treatment-related microbiome changes, relative to control  
pct <- est_amd * est_cons / mean(df5 %>% filter(general_category == "Control") %>% pull(Total.yield)) * 100
pct_yield_inc <- append(pct_yield_inc, pct)
```

### Scenario 6: OR 3-yr bacteria in response to mustard incorporation

``` r
# begin by subsetting the phyloseq object by site, rotation, and year
ps <- bact.ps.list[["OR.16S.ps"]] %>% 
  subset_samples(year == 22 & rotation == 3 & !is.na(Total.yield) & general_category != "Fumigated") 

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples
sub.df <- jim.info.s %>% filter(State == "OR" & Year == 22 & Rotation == 3)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps <- merge_with_jim_data(ps, sub.df)

# they're both in spring
or.must.bact.spring.asvs <- c("ASV270", "ASV2762", "ASV3009", "ASV3745", "ASV574")

# keep only the spring samples, drop, and transform
ps <- ps %>% subset_samples(season == "Spring") %>% 
  drop_ghost_asvs() %>%
  subset_occupancy(0.5) %>% 
  transform_clr() 
```

    ## 5378 features are kept. 
    ## Counts transformed with CLR.

``` r
# add in the abundance of the spring associated ASVs to the sample data 
df6 <- cbind(data.frame(ps@sam_data), ps@otu_table[,or.must.bact.spring.asvs])

# make a consortium column of the sum of the ASVs
df6$consortium <- rowSums(df6[,startsWith(colnames(df6), "ASV")])

# drop non Norkotah cultivars
df6.no <- df6 %>% filter(cultivar == "Norkotah")

# drop NA cols
df6.no <- df6.no[, colSums(is.na(df6.no)) != nrow(df6.no)]

# convert Brassica to integer
df6.no$Brassica <- as.integer(as.logical(df6.no$Brassica))

# SEM: mustard increases nitrate, decreases pH, and increases cumulative target ASV abundances at planting, 
# abundances decrease yields. Best SEM included pH, but not nitrate. Omitted non-Norkotah cultivars.
psem3 <- psem(lm(pH ~ Brassica, df6.no),
              lm(consortium ~ Brassica, df6.no),
              lm(Total.yield ~ consortium + pH, df6.no))

fisherC(psem3) # valid
```

    ##   Fisher.C df P.Value
    ## 1    6.893  4   0.142

``` r
coefs(psem3)
```

    ##      Response  Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ## 1          pH   Brassica  -1.0167    0.4013  7    -2.5332  0.0391      -0.6916
    ## 2  consortium   Brassica   9.3067    1.5021  7     6.1959  0.0004       0.9197
    ## 3 Total.yield consortium -17.0278    4.9085  6    -3.4691  0.0133      -0.9390
    ## 4 Total.yield         pH -62.6043   33.7888  6    -1.8528  0.1133      -0.5015
    ##      
    ## 1   *
    ## 2 ***
    ## 3   *
    ## 4

``` r
# extract p and R2 values
psemrsq <- rsquared(psem3)[which(rsquared(psem3)$Response == "Total.yield"),"R.squared"]
rsq <- append(rsq, psemrsq)
pval <- append(pval, fisherC(psem3, conserve = T)[1,3])

# calculate coefficients
# cultivar left out of consortium -> yield model because insignificant
est_amd <- summary(lm(consortium ~ Brassica + pH, df6.no))$coefficients[2, 1]
est_cons <- summary(lm(Total.yield ~ consortium, df6.no))$coefficients[2, 1]

# Yield increase as an effect of treatment on microbiome (cwt/ac) is the product of the individual coefficients 
yield_inc <- append(yield_inc, (est_amd * est_cons))

# percent increase in total yield from treatment-related microbiome changes, relative to control  
pct <- est_amd * est_cons / mean(df6.no %>% filter(general_category == "Control") %>% pull(Total.yield)) * 100
pct_yield_inc <- append(pct_yield_inc, pct)
```

## Summary table /figure

``` r
# add site, soil treatment, rotation length, and model type 
site <- c("OR", "ID", "ME1", "MI", "OR", "OR")
treatment <- c(rep("Amendment", times = 4), "Fumigation", "Mustard")
rotation <- c("2-yr", "3-yr", "2-yr", "3-yr", "3-yr", "3-yr")
modeltype <- c("lm", "SEM", "lm", "SEM", "lm", "SEM")

# bind these together with the model output data
model.table <- data.frame(cbind(site, treatment, rotation, modeltype, 
                                "p_value" = round(pval, 3), "mult_R2" = round(rsq, 2), 
                                "yield_inc" = round(yield_inc, 1), "pct_yield_inc" = round(pct_yield_inc, 1)))
# coerce back to numeric
for (i in 5:8){
  model.table[,i] <- as.numeric(model.table[,i])
}

# and here we are
model.table
```

    ##   site  treatment rotation modeltype p_value mult_R2 yield_inc pct_yield_inc
    ## 1   OR  Amendment     2-yr        lm   0.005    0.73     115.5          23.2
    ## 2   ID  Amendment     3-yr       SEM   0.275    0.61      52.3          10.4
    ## 3  ME1  Amendment     2-yr        lm   0.001    0.53      31.8          21.3
    ## 4   MI  Amendment     3-yr       SEM   0.925    0.27      56.7          23.4
    ## 5   OR Fumigation     3-yr        lm   0.013    0.61      14.2           2.6
    ## 6   OR    Mustard     3-yr       SEM   0.142    0.67    -140.5         -26.4

``` r
write_tsv(model.table, file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/figures_and_tables/supplemental_tables/microbiome_model_info_table.tsv")
```

## Conclusions

#### session info

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Chicago
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lme4_1.1-35.5        Matrix_1.7-0         piecewiseSEM_2.3.0.1
    ##  [4] compositions_2.0-8   patchwork_1.2.0.9000 speedyseq_0.5.3.9021
    ##  [7] phyloseq_1.48.0      lubridate_1.9.3      forcats_1.0.0       
    ## [10] stringr_1.5.1        dplyr_1.1.4          purrr_1.0.2         
    ## [13] readr_2.1.5          tidyr_1.3.1          tibble_3.2.1        
    ## [16] ggplot2_3.5.1        tidyverse_2.0.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] sandwich_3.1-0          permute_0.9-7           rlang_1.1.4            
    ##   [4] magrittr_2.0.3          multcomp_1.4-26         ade4_1.7-22            
    ##   [7] compiler_4.4.1          mgcv_1.9-1              vctrs_0.6.5            
    ##  [10] reshape2_1.4.4          pkgconfig_2.0.3         crayon_1.5.3           
    ##  [13] fastmap_1.2.0           XVector_0.44.0          utf8_1.2.4             
    ##  [16] rmarkdown_2.28          tzdb_0.4.0              nloptr_2.1.1           
    ##  [19] UCSC.utils_1.0.0        bit_4.0.5               xfun_0.47              
    ##  [22] zlibbioc_1.50.0         GenomeInfoDb_1.40.1     jsonlite_1.8.8         
    ##  [25] biomformat_1.32.0       rhdf5filters_1.16.0     Rhdf5lib_1.26.0        
    ##  [28] parallel_4.4.1          cluster_2.1.6           R6_2.5.1               
    ##  [31] stringi_1.8.4           RColorBrewer_1.1-3      boot_1.3-30            
    ##  [34] car_3.1-2               estimability_1.5.1      Rcpp_1.0.13            
    ##  [37] iterators_1.0.14        knitr_1.48              zoo_1.8-12             
    ##  [40] IRanges_2.38.1          splines_4.4.1           igraph_2.0.3           
    ##  [43] timechange_0.3.0        tidyselect_1.2.1        rstudioapi_0.16.0      
    ##  [46] abind_1.4-5             yaml_2.3.10             MuMIn_1.48.4           
    ##  [49] vegan_2.6-8             codetools_0.2-20        lattice_0.22-6         
    ##  [52] plyr_1.8.9              Biobase_2.64.0          withr_3.0.1            
    ##  [55] evaluate_0.24.0         survival_3.6-4          bayesm_3.1-6           
    ##  [58] Biostrings_2.72.1       pillar_1.9.0            carData_3.0-5          
    ##  [61] DiagrammeR_1.0.11       tensorA_0.36.2.1        foreach_1.5.2          
    ##  [64] stats4_4.4.1            insight_0.20.4          generics_0.1.3         
    ##  [67] vroom_1.6.5             S4Vectors_0.42.1        hms_1.1.3              
    ##  [70] munsell_0.5.1           scales_1.3.0            minqa_1.2.8            
    ##  [73] xtable_1.8-4            glue_1.7.0              emmeans_1.10.4         
    ##  [76] tools_4.4.1             robustbase_0.99-4       data.table_1.16.0      
    ##  [79] visNetwork_2.1.2        mvtnorm_1.3-1           rhdf5_2.48.0           
    ##  [82] grid_4.4.1              ape_5.8                 colorspace_2.1-1       
    ##  [85] nlme_3.1-164            GenomeInfoDbData_1.2.12 performance_0.12.3     
    ##  [88] cli_3.6.3               fansi_1.0.6             gtable_0.3.5           
    ##  [91] DEoptimR_1.1-3          digest_0.6.37           BiocGenerics_0.50.0    
    ##  [94] TH.data_1.1-2           htmlwidgets_1.6.4       farver_2.1.2           
    ##  [97] htmltools_0.5.8.1       multtest_2.60.0         lifecycle_1.0.4        
    ## [100] httr_1.4.7              bit64_4.0.5             MASS_7.3-60.2
