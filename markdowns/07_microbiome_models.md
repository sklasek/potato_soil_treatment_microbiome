SEM
================
Scott Klasek
25 February, 2025

## Purpose

Incorporate ASV abundances into treatment –\> yield models using SEM.
Also, provide some background about directionality about the
treatment-responsive ASVs: did they all increase in the particular
treatment, etc.

Make generalizations about when/where/which treatments affected yields
THROUGH changing the microbiome.

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

### treatment and yield associated ASVs
its.asvs <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/its.asvs.treatment.and.yield.associated.csv")
bact.asvs <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/bact.asvs.treatment.and.yield.associated.csv")
all.trt.yld.asvs <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/all.asvs.treatment.and.yield.associated.filtered.csv")

### just yield-associated ASVs with rotation/season info for context
its.yield.df <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/its.yield.lm.results.txt")
bact.yield.df <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/bact.yield.lm.results.txt")

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
# plot.biomarkers takes a phyloseq object (of counts, not already transformed to percent abundances etc), 
# a list of ASVs, and a taxonomic level to plot them by
# it returns a bar plot showing the percent abundances of those ASVs colored by the specified taxonomic level,
# across all the samples in the ps object. 
plot.biomarkers <- function(ps, biomarkerlist, taxlevel){
  
  # transform counts to percent abundances
  ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU) * 100)
  
  # calculate mean percent abundances for all biomarker ASVs
  num <- vector("numeric") # define a numeric vector
  for (i in biomarkerlist) {num[i] <- mean(otu_table(ps)[,i])}
  
  # obtains the row numbers corresponding to biomarker ASVs
  asvs.to.subset <- vector("integer") # define an integer vector
  for (i in names(num)) {asvs.to.subset[i] <- which(rownames(tax_table(ps))==i)} 
  
  # subset taxonomy and ASV tables
  bmtt <- data.frame(tax_table(ps)[asvs.to.subset,])
  bmtt$ASV <- rownames(bmtt)
  bmtt <- as.matrix(bmtt)
  bmasvt <- otu_table(ps)[,names(num)] 
  
  # make a new phyloseq object
  ps.bm <- phyloseq(tax_table(bmtt), 
                 sample_data(ps),
                 otu_table(bmasvt, taxa_are_rows = FALSE))
  
  # make a barplot of % abundances across all samples
  bm.barplot <- plot_bar(ps.bm, fill=taxlevel)+
    geom_bar(stat="identity", position="stack")+
    scale_y_continuous("% Abundance")+
    theme(axis.text.x = element_blank())+
    facet_grid(~general_category, scales = "free", space = "free")

  return(bm.barplot) # return the plot
}

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

## Inspect direction of treatment influence on treatment-yield ASVs

This narrowed down the ASVs from its.asvs and bact.asvs into the ones
found in all.trt.yld.asvs. See also doc 40.Rmd, which plotted ITS and
16S ASVs separately (before filtering out the ones that decreased in the
particular treatment) and then together.

#### ITS

``` r
# see which rotation they correspond to 
its.yield.df %>% filter(site == "OR" & ASV %in% c("ASV416", "ASV2495"))

# get the phyloseq object
ps1 <- its.ps.list[["ME.ITS.ps"]] %>% 
  subset_samples(year == 22 & !is.na(Total.yield))

# plot the treatment-influenced ASVs
plot.biomarkers(ps1, c("ASV1642"), "ASV")+
  facet_grid(~general_category, scales = "free", space = "free")+
  scale_fill_discrete()
```

#### Inspect Sordariomycete Class

``` r
ps2 <- its.ps.list[["MN.ITS.ps"]] %>% 
  subset_samples(year == 22 & !is.na(Total.yield)) %>% 
  tax_glom(taxrank = "Class")

ps2 <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU) * 100) %>% 
  subset_taxa(Class == "Sordariomycetes")

plot_bar(ps2, fill= "Class")+
    geom_bar(stat="identity", position="stack")+
    scale_y_continuous("% Abundance")+
    theme(axis.text.x = element_blank())+
    facet_grid(~general_category, scales = "free", space = "free")
```

#### 16S by amendment, where we saw the most ASVs.

``` r
# Idaho
asvs.to.plot <- bact.asvs %>% filter(site == "ID" & treatment == "Amended") %>% pull(ASV)

# see which rotation they correspond to 
asvs.to.plot <- bact.yield.df %>% filter(ASV %in% asvs.to.plot & site == "ID" & rotation == 3) %>% pull(ASV)

asvs.to.plot <- c("ASV6")

ps4 <- bact.ps.list[["ID.16S.ps"]] %>% 
  subset_samples(year == 22 & !is.na(Total.yield))

plot.biomarkers(ps4, asvs.to.plot, "ASV")+
  facet_grid(~general_category, scales = "free", space = "free")+
  scale_fill_discrete()

# Oregon 
asvs.to.plot <- bact.asvs %>% filter(site == "OR" & treatment == "Amended") %>% pull(ASV)

# see which rotation they correspond to 
asvs.to.plot <- bact.yield.df %>% filter(ASV %in% asvs.to.plot & site == "OR" & rotation == 2) %>% pull(ASV)

ps4 <- bact.ps.list[["OR.16S.ps"]] %>% 
  subset_samples(year == 22 & rotation == 2 & !is.na(Total.yield))

plot.biomarkers(ps4, asvs.to.plot[18:19], "ASV")+
  facet_grid(~general_category, scales = "free", space = "free")+
  scale_fill_discrete()
```

#### 16S one-offs

``` r
# specify the ASVs to investigate
asvs.to.plot <- c("ASV331")

# see which rotation they correspond to 
bact.yield.df %>% filter(ASV %in% asvs.to.plot & site == "ND")

# subset a phyloseq to include just the samples from that year & rotation
ps4 <- bact.ps.list[["MN.16S.ps"]] %>% 
  subset_samples(year == 22 & !is.na(Total.yield))
  
plot.biomarkers(ps4, asvs.to.plot, "ASV")+
  facet_grid(~general_category, scales = "free", space = "free")+
  scale_fill_discrete()
```

From filtering these, it appeared that there were several bacteria in OR
2-yr and ID 3-yr rotations that increased in amended treatments and were
positively associated with yields. Fumigations also increased yields
here, but without influencing as many ASVs. Only one
fumigation-influenced yield-associated ASV decreased in the fumigated
treatment, which was ASV1428 (Chryseolinea sp, 16S) in ID.

## Tracking consortia over time and space

#### Which season/rotation are these treatment/yield ASVs associated with?

I looked for all combinations of bacterial and eukaryotic ASVs across
all sites, treatments, and rotations, separating ASV by which season
(spring or summer) they were yield-associated with. Results are in the
excel spreadsheet treatment_asv_yield_modeling. One example given here:

``` r
# filter trt-yld-ASVs by site and treatment... ITS
lookup <- all.trt.yld.asvs %>% filter(site == "MI" & treatment == "Amended") %>% pull(ASV)
its.yield.df %>% filter(ASV %in% lookup & lm_mult_r_sq >= 0.2 & site == "MI")
```

    ##       ASV    loglik dfree lm_mult_r_sq lm_adj_r_sq  estimate statistic
    ## 1 ASV1099 -148.2095    27    0.2412358   0.2131334 16.604508  3.673240
    ## 2 ASV1099 -145.0761    26    0.2180398   0.1879644  9.837806  2.407935
    ## 3 ASV1194 -149.0213    27    0.2220518   0.1932389 32.470728  6.070325
    ## 4 ASV1194 -144.3887    26    0.2527223   0.2239809 23.695247  2.587787
    ## 5  ASV161 -148.6650    27    0.2140044   0.1848935 18.944759  4.054420
    ## 6   ASV80 -148.7268    27    0.2104671   0.1812251 16.750033  4.414474
    ##              p        padj site amplicon season rotation      Kingdom
    ## 1 1.043481e-03 0.166746463   MI      ITS spring        3        Fungi
    ## 2 2.343300e-02 0.958545626   MI      ITS summer        3        Fungi
    ## 3 1.758140e-06 0.001872419   MI      ITS spring        3 Heterolobosa
    ## 4 1.560104e-02 0.958545626   MI      ITS summer        3 Heterolobosa
    ## 5 3.831722e-04 0.081615687   MI      ITS spring        3        Fungi
    ## 6 1.467432e-04 0.052093847   MI      ITS spring        3        Fungi
    ##                            Phylum           Class          Order         Family
    ## 1                      Ascomycota Sordariomycetes   Microascales  Microascaceae
    ## 2                      Ascomycota Sordariomycetes   Microascales  Microascaceae
    ## 3 Heterolobosa_phy_Incertae_sedis   Heterolobosea Schizopyrenida Vahlkampfiidae
    ## 4 Heterolobosa_phy_Incertae_sedis   Heterolobosea Schizopyrenida Vahlkampfiidae
    ## 5                      Ascomycota Sordariomycetes   Microascales  Microascaceae
    ## 6                      Ascomycota Sordariomycetes    Sordariales  Chaetomiaceae
    ##             Genus      Species                closest_tax
    ## 1          Kernia   columnaris          Kernia columnaris
    ## 2          Kernia   columnaris          Kernia columnaris
    ## 3 Allovahlkampfia         <NA>        Allovahlkampfia sp.
    ## 4 Allovahlkampfia         <NA>        Allovahlkampfia sp.
    ## 5          Kernia   columnaris          Kernia columnaris
    ## 6   Botryotrichum spirotrichum Botryotrichum spirotrichum

``` r
# same for 16S
bact.yield.df %>% filter(ASV %in% lookup & lm_mult_r_sq >= 0.2 & site == "MI")
```

    ##       ASV    loglik dfree lm_mult_r_sq lm_adj_r_sq estimate statistic
    ## 1 ASV1426 -147.0370    27    0.2943043   0.2681674 21.23211  3.197438
    ## 2  ASV534 -138.2291    26    0.2228267   0.1929355 23.77813  2.405470
    ##             p      padj site amplicon season rotation  Kingdom           Phylum
    ## 1 0.003522083 0.5950911   MI      16S spring        3 Bacteria       Firmicutes
    ## 2 0.023561945 0.4235070   MI      16S spring        2 Bacteria Actinobacteriota
    ##            Class                               Order                Family
    ## 1     Clostridia Peptostreptococcales-Tissierellales Peptostreptococcaceae
    ## 2 Actinobacteria                       Micrococcales Promicromonosporaceae
    ##                Genus Species            closest_tax
    ## 1         Romboutsia    <NA>         Romboutsia sp.
    ## 2 Cellulosimicrobium    <NA> Cellulosimicrobium sp.

#### OR and ID amendment-associated bacterial consortia

``` r
# get vectors of OR or ID amendment/yield-associated bacterial ASVs
or.amd.bact.asvs <- all.trt.yld.asvs %>% filter(site == "OR" & treatment == "Amended") %>% pull(ASV)
id.amd.bact.asvs <- all.trt.yld.asvs %>% filter(site == "ID" & treatment == "Amended") %>% pull(ASV)

# separate them by rotation
or.amd.bact.2yr.asvs <- bact.yield.df %>% filter(ASV %in% or.amd.bact.asvs & site == "OR" & rotation == 2) %>% pull(ASV)
or.amd.bact.3yr.asvs <- bact.yield.df %>% filter(ASV %in% or.amd.bact.asvs & site == "OR" & rotation == 3) %>% pull(ASV)
id.amd.bact.3yr.asvs <- bact.yield.df %>% filter(ASV %in% id.amd.bact.asvs & site == "ID" & rotation == 3) %>% pull(ASV) 

### oregon bacterial consortium across BOTH rotations
# subset phyloseq
ps1 <- bact.ps.list[["OR.16S.ps"]] %>% 
  subset_samples(!is.na(Total.yield))

# add a year/treatment label
ps1@sam_data$year_cat <- paste(ps1@sam_data$year, ps1@sam_data$general_category, sep = "_")
  
# plot over treatment and time
plot.biomarkers(ps1, union(or.amd.bact.2yr.asvs, or.amd.bact.3yr.asvs), "ASV")+
  facet_grid(~year_cat, scales = "free", space = "free")+
  scale_fill_discrete()+
  theme_bw()+ggtitle("OR amendment/yield bacteria")+
  theme(axis.text.x = element_blank())
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
### idaho 3-yr bacterial consortium
# subset phyloseq
ps2 <- bact.ps.list[["ID.16S.ps"]] %>% 
  subset_samples(rotation == 3 & !is.na(Total.yield))

# remember that Green Manure treatments are now classified as controls
ps2@sam_data[which(ps2@sam_data$amendment == "Green Manure"),"general_category"] <- "Control"

# add a year/treatment label
ps2@sam_data$year_cat <- paste(ps2@sam_data$year, ps2@sam_data$general_category, sep = "_")
  
# plot over treatment and time
id.3yr.amd.bact.gg <- plot.biomarkers(ps2, unique(id.amd.bact.3yr.asvs), "ASV")+
  facet_grid(~year_cat, scales = "free", space = "free")+
  scale_fill_discrete()+
  theme_bw()+ggtitle("ID 3-yr amendment/yield bacteria")+
  theme(axis.text.x = element_blank())
id.3yr.amd.bact.gg
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-10-2.png)<!-- --> In
both OR and ID, amendments stimulate these ASVs in the second potato
year and increase their combined relative abundances from ~0.1% to
~0.3-0.4% of the community.

In MN1, ASV792 also increased in amendment treatments by the second
potato year. This and three other amendment-yield-positive ASVs were
found across 6 sites total, so it seems they only appear to increase
under certain treatments/circumstances. See what I mean here:

``` r
all.trt.yld.asvs %>% filter(ASV %in% or.amd.bact.2yr.asvs) %>% dplyr::select(site, ASV, num_sites)
```

    ##    site      ASV num_sites
    ## 1   MN1   ASV792         7
    ## 2    OR   ASV792         7
    ## 3    OR  ASV1254         7
    ## 4    OR  ASV2399         6
    ## 5    OR  ASV3632         4
    ## 6    OR  ASV4948         2
    ## 7    OR  ASV5076         4
    ## 8    OR  ASV6264         3
    ## 9    OR  ASV7919         2
    ## 10   OR ASV21630         2

#### OR brassica-yield-associated bacterial consortia

``` r
# get vectors of OR brassica/yield-associated bacterial ASVs
or.brs.bact.asvs <- all.trt.yld.asvs %>% filter(site == "OR" & treatment == "Must.") %>% pull(ASV)

# all but one are 3-yr
or.brs.bact.3yr.asvs <- bact.yield.df %>% filter(ASV %in% or.brs.bact.asvs & site == "OR" & rotation == 3) %>% pull(ASV)

# subset phyloseq
ps3 <- bact.ps.list[["OR.16S.ps"]] %>% 
  subset_samples(rotation == 3 & !is.na(Total.yield))

# add a year/treatment label
ps3@sam_data$year_cat <- paste(ps3@sam_data$year, ps3@sam_data$treatment_description, sep = "_")
  
# plot over treatment and time
plot.biomarkers(ps3, or.brs.bact.3yr.asvs, "ASV")+
  facet_grid(~year_cat, scales = "free", space = "free")+
  scale_fill_discrete()+
  theme_bw()+ggtitle("OR 3-yr brassica/yield bacteria")+
  theme(axis.text.x = element_blank())
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

These brassica-ASVs are all negatively-associated with yield, and
clearly only in the second potato year.

#### ME amendment-yield-associated fungi

``` r
# get the asvs
me.its.amd.asvs <- its.asvs %>% filter(site == "ME1" & treatment == "Amended") %>% pull(ASV)

# all are 2-yr
me.its.amd.asvs <- unique(its.yield.df %>% filter(ASV %in% me.its.amd.asvs & site == "ME" & rotation == 2) %>% pull(ASV))

# subset phyloseq
ps4 <- its.ps.list[["ME.ITS.ps"]] %>% 
  subset_samples(rotation == 2 & !is.na(Total.yield))

# add a year/treatment label
ps4@sam_data$year_cat <- paste(ps4@sam_data$year, ps4@sam_data$general_category, sep = "_")
  
# plot over treatment and time
plot.biomarkers(ps4, me.its.amd.asvs, "ASV")+
  facet_grid(~year_cat, scales = "free", space = "free")+
  scale_fill_discrete()+
  theme_bw()+ggtitle("ME 2-yr amended/yield fungi")+
  theme(axis.text.x = element_blank())
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

These Ascomycete fungi increased over the two potato growing years, but
were already higher in the first year. Notice though that their
occupancy increases in the final year.

#### MI amendment-yield-associated euks

``` r
# get asvs
mi.its.amd.asvs <- its.asvs %>% filter(site == "MI" & treatment == "Amended") %>% pull(ASV)

# all are 3-yr
mi.its.amd.asvs <- unique(its.yield.df %>% filter(ASV %in% mi.its.amd.asvs & site == "MI" & rotation == 3) %>% pull(ASV))

# subset phyloseq
ps5 <- its.ps.list[["MI.ITS.ps"]] %>% 
  subset_samples(rotation == 3 & !is.na(Total.yield))

# add a year/treatment label
ps5@sam_data$year_cat <- paste(ps5@sam_data$year, ps5@sam_data$general_category, sep = "_")
  
# plot over treatment and time
plot.biomarkers(ps5, mi.its.amd.asvs, "ASV")+
  facet_grid(~year_cat, scales = "free", space = "free")+
  scale_fill_discrete()+
  theme_bw()+ggtitle("MI 3-yr amended/yield euks")+
  theme(axis.text.x = element_blank())
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

These euks in MI increased with treatment and time, but seemed to have
had a bit of a head start in 2019.

### Pause for a mini-conclusion

Now we have a sense of WHICH yield-associated ASVs increased WHERE, in
response to WHAT TREATMENTS. Next we should establish a causal model,
particularly for the Amendment-associated treatments. Does amendment
increase organic matter or other nutrients in the soil, which then
increase abundances of these taxa? Or does amendment inoculate the soil
with these taxa? (We might not be able to answer that). Or rather, how
much does OM or nutrients explain the yield increase, and how much is it
the taxa themselves? And are the members of these consortia part of a
co-occurring module?

For the brassica treatment, I suppose OM or C or pH or something else
could be used, but I still don’t know how brassica relates to yields.

For fumigated treatments, this is tricky because few soil taxa actually
varied with this treatment. The only yield-associated ASV that decreased
in fumigated treatments does not seem to be a pathogen (on the contrary,
it’s implicated as a beneficial bacterium in combating fusarium wilt).
So how would fumigation change yields through the microbiome? By
changing something like alpha diversity or bacteria:fungi ratios?

## Structural equation modeling (SEM)

### Scenario 1: Oregon bacterial consortia in response to amendment

From the 2-yr rotation only

#### Make a dataframe with the ASV, yield, and soil chemical data across samples.

``` r
# begin by subsetting the phyloseq object by site, rotation, and year
ps6 <- bact.ps.list[["OR.16S.ps"]] %>% 
  subset_samples(year == 22 & rotation == 2 & !is.na(Total.yield) & general_category != "Fumigated") 

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples 
sub.df <- jim.info.s %>% filter(State == "OR" & Year == 22)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps6 <- merge_with_jim_data(ps6, sub.df)

# see which season the ASVs whose abundances we want to include were yield-associated with
bact.yield.df %>% filter(ASV %in% or.amd.bact.2yr.asvs & site == "OR") %>% dplyr::select(ASV, season) 
```

    ##        ASV season
    ## 1  ASV1254 summer
    ## 2 ASV21630 summer
    ## 3  ASV2399 summer
    ## 4  ASV3632 summer
    ## 5  ASV4948 summer
    ## 6  ASV5076 summer
    ## 7  ASV6264 spring
    ## 8  ASV7919 summer
    ## 9   ASV792 summer

``` r
# subset the ASVs by season, then subset ASVs below 50% occupancy and transform counts to CLR
or.amd.bact.summer.asvs <- bact.yield.df %>% 
  filter(ASV %in% or.amd.bact.2yr.asvs & site == "OR" & season == "summer") %>% 
  pull(ASV)

# keep only the summer samples
ps6 <- ps6 %>% subset_samples(season == "Summer")%>% 
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
```

#### Check our premises

``` r
# is yield higher in amended?
ggplot(df1, aes(Amended, Total.yield))+geom_jitter(width = 0.1) # yes
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
# is consortium higher in amended? 
ggplot(df1, aes(Amended, consortium))+geom_jitter(width = 0.2) # yes
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
# is OM higher in amended?
ggplot(df1, aes(Amended, OM....))+geom_jitter(width = 0.1) # yes (one outlier)
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

``` r
# also Solvita, pH, P were noticeably higher in amended
ggplot(df1, aes(Amended, Solvita..ppm.))+geom_jitter(width = 0.1) # yes
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-16-4.png)<!-- -->

``` r
# convert amended treatment to 0/1
df1$Amended <- as.integer(as.logical(df1$Amended))
```

#### Modeling

``` r
### tidying up the dataframe 
# remove NA cols
df1 <- df1 %>% dplyr::select(-cover_2020, -green_cov_pct_pre_kill, -Total.C...., -Total.organic.C....)

# there was one outlier OM value at 2.8, next highest OM recorded in Oregon any time was 1.8
# also keep only Norkotah cultivar
df1.no <- df1 %>% filter(`OM....` < 2 & cultivar == "Norkotah")

### evaluate linear models to include in the SEM
# first equation: soil organic amendment increases OM%
# since amendment is encoded as a 0/1, can use linear model
amend <- lm(`OM....` ~ Amended, df1.no) 
summary(amend) # it's NOT a bad model
```

    ## 
    ## Call:
    ## lm(formula = OM.... ~ Amended, data = df1.no)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.31667 -0.06667  0.02000  0.10167  0.18333 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.31667    0.06616  19.903 9.48e-09 ***
    ## Amended      0.36333    0.09812   3.703   0.0049 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.162 on 9 degrees of freedom
    ## Multiple R-squared:  0.6037, Adjusted R-squared:  0.5597 
    ## F-statistic: 13.71 on 1 and 9 DF,  p-value: 0.004899

``` r
summary(lm(Total.yield ~ Amended, df1.no)) # Amendment itself DOES explain total yield
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ Amended, data = df1.no)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -86.672 -39.358   3.651  36.552  93.013 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   498.12      25.33   19.67 1.05e-08 ***
    ## Amended       102.93      37.57    2.74   0.0229 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 62.04 on 9 degrees of freedom
    ## Multiple R-squared:  0.4547, Adjusted R-squared:  0.3942 
    ## F-statistic: 7.506 on 1 and 9 DF,  p-value: 0.02286

``` r
# next equation: OM increases soil respiration
ompct <- lm(`Solvita..ppm.` ~ `OM....`, df1.no)
summary(ompct) # it's a very good model
```

    ## 
    ## Call:
    ## lm(formula = Solvita..ppm. ~ OM...., data = df1.no)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -22.218 -12.369   8.081  10.733  16.881 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -63.99      29.74  -2.152 0.059847 .  
    ## OM....        102.01      19.82   5.146 0.000607 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 15.31 on 9 degrees of freedom
    ## Multiple R-squared:  0.7463, Adjusted R-squared:  0.7181 
    ## F-statistic: 26.48 on 1 and 9 DF,  p-value: 0.0006068

``` r
summary(lm(Total.yield ~ `OM....`, df1.no)) # OM does not explain total yield
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ OM...., data = df1.no)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -115.858  -54.984    0.893   48.553  132.676 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)    401.5      155.9   2.576   0.0299 *
    ## OM....          96.8      103.9   0.932   0.3759  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 80.24 on 9 degrees of freedom
    ## Multiple R-squared:  0.08795,    Adjusted R-squared:  -0.01339 
    ## F-statistic: 0.8678 on 1 and 9 DF,  p-value: 0.3759

``` r
summary(lm(Total.yield ~ `Solvita..ppm.`, df1.no)) # respiration does not explain total yield
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ Solvita..ppm., data = df1.no)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -102.88  -48.93  -13.99   39.60  111.19 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   420.1399    72.0488   5.831  0.00025 ***
    ## Solvita..ppm.   1.4314     0.7883   1.816  0.10278    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 71.88 on 9 degrees of freedom
    ## Multiple R-squared:  0.2681, Adjusted R-squared:  0.1868 
    ## F-statistic: 3.297 on 1 and 9 DF,  p-value: 0.1028

``` r
# soil respiration increases consortium abundance
resp <- lm(consortium ~ `Solvita..ppm.`, df1.no)
summary(resp) # it's a good model
```

    ## 
    ## Call:
    ## lm(formula = consortium ~ Solvita..ppm., data = df1.no)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7913 -2.2278 -0.1764  1.2971  4.1348 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -11.13128    2.67356  -4.163 0.002435 ** 
    ## Solvita..ppm.   0.16920    0.02925   5.784 0.000265 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.667 on 9 degrees of freedom
    ## Multiple R-squared:  0.788,  Adjusted R-squared:  0.7645 
    ## F-statistic: 33.46 on 1 and 9 DF,  p-value: 0.0002647

``` r
# next, consortium abundance increases total yield
cons <- lm(Total.yield ~ consortium, df1.no)
summary(cons) 
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ consortium, data = df1.no)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -72.859 -38.915  -2.148  37.332  75.540 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  504.050     19.302  26.114 8.55e-10 ***
    ## consortium    11.298      3.032   3.727  0.00472 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 52.69 on 9 degrees of freedom
    ## Multiple R-squared:  0.6068, Adjusted R-squared:  0.5631 
    ## F-statistic: 13.89 on 1 and 9 DF,  p-value: 0.004721

``` r
# BACK UP: does amendment explain consortium?
trt.on.cons <- lm(consortium ~ Amended, df1.no) # yes indeedy
summary(trt.on.cons)
```

    ## 
    ## Call:
    ## lm(formula = consortium ~ Amended, data = df1.no)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7468 -1.0002 -0.3676  1.0246  2.0530 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -1.0312     0.5600  -1.841   0.0987 .  
    ## Amended      10.2247     0.8306  12.309  6.2e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.372 on 9 degrees of freedom
    ## Multiple R-squared:  0.9439, Adjusted R-squared:  0.9377 
    ## F-statistic: 151.5 on 1 and 9 DF,  p-value: 6.197e-07

``` r
# add amendment back into the consortium-yield model
cons2 <- lm(Total.yield ~ consortium + Amended, df1.no)
summary(cons2)
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
    ## Amended      -224.50     118.84  -1.889   0.0956 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 46.47 on 8 degrees of freedom
    ## Multiple R-squared:  0.7281, Adjusted R-squared:  0.6601 
    ## F-statistic: 10.71 on 2 and 8 DF,  p-value: 0.005467

``` r
### Build SEMs
# SEM 1: amendment increases OM, which increases respiration, which increases consortium, then total yields
psem1 <- psem(amend, ompct, resp, cons)  
coefs(psem1)
```

    ##        Response     Predictor Estimate Std.Error DF Crit.Value P.Value
    ## 1        OM....       Amended   0.3633    0.0981  9     3.7028  0.0049
    ## 2 Solvita..ppm.        OM.... 102.0061   19.8242  9     5.1455  0.0006
    ## 3    consortium Solvita..ppm.   0.1692    0.0293  9     5.7840  0.0003
    ## 4   Total.yield    consortium  11.2978    3.0316  9     3.7267  0.0047
    ##   Std.Estimate    
    ## 1       0.7770  **
    ## 2       0.8639 ***
    ## 3       0.8877 ***
    ## 4       0.7790  **

``` r
fisherC(psem1) # rats it's bad. perhaps its overfitted?
```

    ##   Fisher.C df P.Value
    ## 1   43.078 12       0

``` r
# plot(psem1)

# SEM 2: amendment increases consortium, which increases total yields
psem2 <- psem(trt.on.cons, cons)  
coefs(psem2)
```

    ##      Response  Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ## 1  consortium    Amended  10.2247    0.8306  9    12.3094  0.0000       0.9716
    ## 2 Total.yield consortium  11.2978    3.0316  9     3.7267  0.0047       0.7790
    ##      
    ## 1 ***
    ## 2  **

``` r
fisherC(psem2) # bad
```

    ##   Fisher.C df P.Value
    ## 1    4.696  2   0.096

``` r
# plot(psem2)

# SEM 3: amendment increases consortium and soil chem, soil chem increases(?) consortium, 
# and amendment and consortium increase yield
# didnt get anywhere because none of the soil chem parameters increases consortium independently of Amendment
psem3 <- psem(lm(OM.... ~ Amended, df1.no),
              lm(consortium ~ Amended + OM...., df1.no),
              lm(Total.yield ~ consortium + Amended, df1.no))
fisherC(psem3)
```

    ##   Fisher.C df P.Value
    ## 1    0.936  2   0.626

``` r
coefs(psem3)
```

    ##      Response  Predictor  Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ## 1      OM....    Amended    0.3633    0.0981  9     3.7028  0.0049       0.7770
    ## 2  consortium    Amended   11.8586    1.1857  8    10.0012  0.0000       1.1268
    ## 3  consortium     OM....   -4.4970    2.5357  8    -1.7735  0.1141      -0.1998
    ## 4 Total.yield consortium   32.0237   11.2926  8     2.8358  0.0220       2.2080
    ## 5 Total.yield    Amended -224.5035  118.8437  8    -1.8891  0.0956      -1.4708
    ##      
    ## 1  **
    ## 2 ***
    ## 3    
    ## 4   *
    ## 5

``` r
summary(psem3)$AIC # 155 for OM
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%

    ##       AIC    AICc  K  n
    ## 1 155.457 172.219 11 11

``` r
# scenario 4 is not an SEM, rather a lm where amendment directly impacts yield as well (this cannot be modeled in SEM)
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
    ## Amended      -224.50     118.84  -1.889   0.0956 .  
    ## consortium     32.02      11.29   2.836   0.0220 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 46.47 on 8 degrees of freedom
    ## Multiple R-squared:  0.7281, Adjusted R-squared:  0.6601 
    ## F-statistic: 10.71 on 2 and 8 DF,  p-value: 0.005467

``` r
extractAIC(lm4) # gives AIC of 86, which is lower than the PSEM
```

    ## [1]  3.00000 86.95183

The initial approach I took was to examine OR 2-yr rotation samples in
2022 only, using the “consortium” of the bacterial ASVs that were
positively associated with organic amendments and yields in the 2-yr
rotation. I build two SEMs of lms, one where amendment –\> soil OM% –\>
respiration –\> consortium abundance –\> Total yield, and another where
amendment –\> consortium –\> Total yield. I could see the first one
being considered as overfitted, but both SEMs I made had very low (ie
BAD) p-values, so maybe their power was too low anyway.

Next I went back and plotted OR samples in 2022, from both rotations,
using the combination of bacterial ASVs that were positively associated
with organic amendments and yields across both rotations. These ASVs
were still much higher in the 2022 amended samples regardless of
rotation, so I thought it fair to combine them into one analysis.
However, by including samples (and ASVs) from both rotations, the
association between consortium abundance and total yield became
insignificant. I tweaked this to remove the one negative ASV, and then
the model became just-barely significant.

I also checked to see whether Amendment also directly affected Total
yield, but this yielded NA’s for p-value and Fisher C. THIS CAN BE
BETTER MODELED WITH LM: Total.yield ~ Amended + consortium. Curiously,
Amendment had a slightly negative relationship with yield, but at the
same time increased the consortium which increased yields. Either way, I
probably shouldn’t bother justifying the addition of the 3-yr rotation
where amendment did nothing for yields: we just don’t have enough
observations to make a good model in the 2-yr rotation where amendment
actually increased yields.

The reason I saw direct relationships between treatments and yields (in
figure 1) was that I was comparing each treatment accordingly to its
control. Previously I had not been doing this in the SEM section (I was
comparing Amendments to non-Amendments, and the latter included
Fumigated samples).

My conclusion here is that amendment increases the consortium, which
increases the yield MORESO than the amendment’s direct effects on
yields, which were actually negative. Amendment increased OM, but OM did
not increase consortium abundance, nor did solvita, N or P. Something
else that is NOT related to these soil chem parameters is stimulating
this yield-positive consortium.

### Scenario 2: Idaho 3-yr rotation bacterial consortia in response to amendment

``` r
# begin by subsetting the phyloseq object by site, rotation, and year
ps7 <- bact.ps.list[["ID.16S.ps"]] %>% 
  subset_samples(year == 22 & rotation == 3 & !is.na(Total.yield) & Fumigated == FALSE) 

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples 
sub.df <- jim.info.s %>% filter(State == "ID" & Year == 22 & Rotation == 3)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps7 <- merge_with_jim_data(ps7, sub.df)

# All the ASVs were summer-associated with yield
bact.yield.df %>% filter(ASV %in% id.amd.bact.3yr.asvs & site == "ID") %>% dplyr::select(ASV, season) 
```

    ##         ASV season
    ## 1  ASV12474 summer
    ## 2  ASV17047 summer
    ## 3   ASV2340 summer
    ## 4   ASV2345 summer
    ## 5   ASV2589 summer
    ## 6   ASV3781 spring
    ## 7   ASV3781 summer
    ## 8   ASV4416 summer
    ## 9   ASV4825 summer
    ## 10  ASV4830 summer

``` r
# subset the ASVs by season
id.amd.bact.summer.asvs <- bact.yield.df %>% 
  filter(ASV %in% id.amd.bact.3yr.asvs & site == "ID" & season == "summer") %>% 
  pull(ASV)

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

### check our premises
# is yield higher in amended?
ggplot(df2, aes(Amended, Total.yield))+geom_jitter(width = 0.1) # yeah
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
# is consortium higher in amended? 
ggplot(df2, aes(Amended, consortium))+geom_jitter(width = 0.2) # yes
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
# is OM higher in amended?
ggplot(df2, aes(Amended, OM....))+geom_jitter(width = 0.1) # yes 
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->

``` r
# convert amended treatment to 0/1
df2$Amended <- as.integer(as.logical(df2$Amended))


### tidying up the dataframe 
# remove NA cols
df2 <- df2 %>% dplyr::select(-cover_2019, -green_cov_pct_pre_kill, -Total.C...., 
                             -Total.organic.C...., -VPPG, -Root.Lesion.100.cc)

### evaluate linear models to include in the SEM
# first equation: soil organic amendment increases OM % (but not yield)
# since amendment is encoded as a 0/1, can use linear model
amend <- lm(`OM....` ~ Amended, df2) 
summary(amend) # amendment increases OM
```

    ## 
    ## Call:
    ## lm(formula = OM.... ~ Amended, data = df2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.07500 -0.06563 -0.01875  0.03750  0.13750 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.26250    0.02893   43.64 2.32e-16 ***
    ## Amended      0.11250    0.04092    2.75   0.0157 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08183 on 14 degrees of freedom
    ## Multiple R-squared:  0.3506, Adjusted R-squared:  0.3043 
    ## F-statistic:  7.56 on 1 and 14 DF,  p-value: 0.01566

``` r
summary(lm(Total.yield ~ Amended, df2)) # Amendment itself explains total yield
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ Amended, data = df2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -37.812 -14.137   0.838  15.575  32.188 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  503.312      7.482   67.27  < 2e-16 ***
    ## Amended       57.350     10.581    5.42 9.03e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 21.16 on 14 degrees of freedom
    ## Multiple R-squared:  0.6772, Adjusted R-squared:  0.6542 
    ## F-statistic: 29.37 on 1 and 14 DF,  p-value: 9.032e-05

``` r
# OM increases consortium abundance and yield
org <- lm(consortium ~ `OM....`, df2)
summary(org) # it's a good model
```

    ## 
    ## Call:
    ## lm(formula = consortium ~ OM...., data = df2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -9.8101 -3.3047 -0.6798  3.7169  9.8511 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   -38.50      19.00  -2.027   0.0622 .
    ## OM....         33.20      14.37   2.311   0.0366 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 5.459 on 14 degrees of freedom
    ## Multiple R-squared:  0.2761, Adjusted R-squared:  0.2244 
    ## F-statistic:  5.34 on 1 and 14 DF,  p-value: 0.0366

``` r
summary(lm(Total.yield ~ `OM....`, df2)) 
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ OM...., data = df2)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -54.11 -26.13  10.34  22.75  53.88 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   330.36     117.82   2.804   0.0141 *
    ## OM....        152.89      89.11   1.716   0.1083  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 33.86 on 14 degrees of freedom
    ## Multiple R-squared:  0.1737, Adjusted R-squared:  0.1147 
    ## F-statistic: 2.944 on 1 and 14 DF,  p-value: 0.1083

``` r
# next, consortium abundance increases total yield 
cons <- lm(Total.yield ~ consortium, df2)
summary(cons) # 
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ consortium, data = df2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -36.171 -18.103   7.922  17.676  27.559 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 507.9411     7.7092  65.887  < 2e-16 ***
    ## consortium    4.5492     0.9639   4.719 0.000329 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 23.14 on 14 degrees of freedom
    ## Multiple R-squared:  0.614,  Adjusted R-squared:  0.5865 
    ## F-statistic: 22.27 on 1 and 14 DF,  p-value: 0.0003289

``` r
amd.to.cons <- lm(consortium ~ Amended, df2)
summary(amd.to.cons) # amendment DOES influence consortium
```

    ## 
    ## Call:
    ## lm(formula = consortium ~ Amended, data = df2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8389 -1.2785  0.0027  1.0224  3.4758 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -0.4670     0.6468  -0.722    0.482    
    ## Amended      11.5057     0.9147  12.578  5.1e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.829 on 14 degrees of freedom
    ## Multiple R-squared:  0.9187, Adjusted R-squared:  0.9129 
    ## F-statistic: 158.2 on 1 and 14 DF,  p-value: 5.099e-09

``` r
# so make another model where amendment and OM both increase consortium
om.amd <- lm(consortium ~ `OM....` + Amended, df2)
summary(om.amd) 
```

    ## 
    ## Call:
    ## lm(formula = consortium ~ OM.... + Amended, data = df2)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -3.146 -0.966 -0.096  0.897  3.168 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    4.708      7.724   0.610    0.553    
    ## OM....        -4.099      6.095  -0.673    0.513    
    ## Amended       11.967      1.158  10.334 1.23e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.866 on 13 degrees of freedom
    ## Multiple R-squared:  0.9214, Adjusted R-squared:  0.9094 
    ## F-statistic: 76.24 on 2 and 13 DF,  p-value: 6.588e-08

``` r
# amendment increases yield through consortium
summary(lm(Total.yield ~ consortium + OM...., df2))
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ consortium + OM...., data = df2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -35.926 -18.239   8.094  17.664  27.751 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  504.672     95.032   5.311 0.000141 ***
    ## consortium     4.528      1.176   3.851 0.002002 ** 
    ## OM....         2.565     74.283   0.035 0.972982    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 24.01 on 13 degrees of freedom
    ## Multiple R-squared:  0.6141, Adjusted R-squared:  0.5547 
    ## F-statistic: 10.34 on 2 and 13 DF,  p-value: 0.002052

``` r
### Build SEMs
# SEM 1: amendment increases OM, which increases consortium, then total yields
psem1 <- psem(amend, org, cons)  
coefs(psem1)
```

    ##      Response  Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ## 1      OM....    Amended   0.1125    0.0409 14     2.7495  0.0157       0.5922
    ## 2  consortium     OM....  33.2001   14.3676 14     2.3108  0.0366       0.5254
    ## 3 Total.yield consortium   4.5492    0.9639 14     4.7195  0.0003       0.7836
    ##      
    ## 1   *
    ## 2   *
    ## 3 ***

``` r
fisherC(psem1) # model is bad
```

    ##   Fisher.C df P.Value
    ## 1   36.948  6       0

``` r
# plot(psem1)

# SEM 2: amendment increases OM, OM and amendment increase consortium, and consortium increases total yields
psem2 <- psem(amend, om.amd, cons)  
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
fisherC(psem2) # it's good, though links to OM are not significant
```

    ##   Fisher.C df P.Value
    ## 1    5.122  4   0.275

``` r
# plot(psem2)

# SEM 3: amendment and OM increase consortium, consortium increases yields
psem3 <- psem(om.amd, cons)
coefs(psem3)
```

    ##      Response  Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ## 1  consortium     OM....  -4.0992    6.0953 13    -0.6725  0.5130      -0.0649
    ## 2  consortium    Amended  11.9669    1.1580 13    10.3340  0.0000       0.9969
    ## 3 Total.yield consortium   4.5492    0.9639 14     4.7195  0.0003       0.7836
    ##      
    ## 1    
    ## 2 ***
    ## 3 ***

``` r
fisherC(psem3) # not as good as SEM2
```

    ##   Fisher.C df P.Value
    ## 1    4.082  4   0.395

``` r
# SEM 4: amendment increases consortium, consortium increases yields
psem4 <- psem(amd.to.cons, cons)  
coefs(psem4)
```

    ##      Response  Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ## 1  consortium    Amended  11.5057    0.9147 14    12.5785   0e+00       0.9585
    ## 2 Total.yield consortium   4.5492    0.9639 14     4.7195   3e-04       0.7836
    ##      
    ## 1 ***
    ## 2 ***

``` r
fisherC(psem4) 
```

    ##   Fisher.C df P.Value
    ## 1    4.027  2   0.134

``` r
# plot(psem4) 

# which SEM is best? 
summary(psem2)$AIC # lower AIC: better model 
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

    ##       AIC   AICc  K  n
    ## 1 189.023 196.66 10 16

``` r
summary(psem3)$AIC 
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

    ##       AIC    AICc K  n
    ## 1 219.853 225.489 7 16

``` r
# summarize our best model
summary(psem2)
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

    ## 
    ## Structural Equation Model of psem2 
    ## 
    ## Call:
    ##   OM.... ~ Amended
    ##   consortium ~ OM.... + Amended
    ##   Total.yield ~ consortium
    ## 
    ##     AIC
    ##  189.023
    ## 
    ## ---
    ## Tests of directed separation:
    ## 
    ##                Independ.Claim Test.Type DF Crit.Value P.Value 
    ##   Total.yield ~ Amended + ...      coef 13     1.6003  0.1335 
    ##    Total.yield ~ OM.... + ...      coef 12    -0.5712  0.5784 
    ## 
    ## --
    ## Global goodness-of-fit:
    ## 
    ## Chi-Squared = 3.306 with P-value = 0.191 and on 2 degrees of freedom
    ## Fisher's C = 5.122 with P-value = 0.275 and on 4 degrees of freedom
    ## 
    ## ---
    ## Coefficients:
    ## 
    ##      Response  Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ##        OM....    Amended   0.1125    0.0409 14     2.7495  0.0157       0.5922
    ##    consortium     OM....  -4.0992    6.0953 13    -0.6725  0.5130      -0.0649
    ##    consortium    Amended  11.9669    1.1580 13    10.3340  0.0000       0.9969
    ##   Total.yield consortium   4.5492    0.9639 14     4.7195  0.0003       0.7836
    ##      
    ##     *
    ##      
    ##   ***
    ##   ***
    ## 
    ##   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05
    ## 
    ## ---
    ## Individual R-squared:
    ## 
    ##      Response method R.squared
    ##        OM....   none      0.35
    ##    consortium   none      0.92
    ##   Total.yield   none      0.61

``` r
# test directed separation for psem2
summary(lm(Total.yield ~ `OM....` + consortium, df2))
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ OM.... + consortium, data = df2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -35.926 -18.239   8.094  17.664  27.751 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  504.672     95.032   5.311 0.000141 ***
    ## OM....         2.565     74.283   0.035 0.972982    
    ## consortium     4.528      1.176   3.851 0.002002 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 24.01 on 13 degrees of freedom
    ## Multiple R-squared:  0.6141, Adjusted R-squared:  0.5547 
    ## F-statistic: 10.34 on 2 and 13 DF,  p-value: 0.002052

``` r
summary(lm(Total.yield ~ Amended + consortium, df2))
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ Amended + consortium, data = df2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -38.275 -14.310   0.864  15.532  32.360 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 503.1399     7.9039  63.657   <2e-16 ***
    ## Amended      61.6025    38.4939   1.600    0.134    
    ## consortium   -0.3696     3.2068  -0.115    0.910    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 21.95 on 13 degrees of freedom
    ## Multiple R-squared:  0.6776, Adjusted R-squared:  0.628 
    ## F-statistic: 13.66 on 2 and 13 DF,  p-value: 0.0006381

``` r
# they suggest there should NOT be paths directly from OM or Amendment directly to Total yield. 

# now that green manure is a "control", can we improve this model?
psem2.1 <- psem(lm(P.Bray..ppm. ~ Amended, df2),
                lm(pH ~ Amended, df2),
                lm(consortium ~ P.Bray..ppm. + pH, df2),
                lm(Total.yield ~ consortium, df2))  
fisherC(psem2.1) 
```

    ##   Fisher.C df P.Value
    ## 1   30.084 10   0.001

``` r
summary(psem2.1)$AIC
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%

    ##       AIC    AICc  K  n
    ## 1 337.373 347.009 13 16

Here in Idaho in the 3-yr rotation, Amendment increases consortium
abundance, and consortium abundance increases yield. Now that I
understand we want a LOWER AIC (lol), keeping OM gives the best model
even though the paths between amendment and OM, and OM and consortium
abundance are not themselves statistically significant at p \< 0.05. The
interpretation would be that the amendment increased the abundance of
consortia members and that the increase in consortia stimulated tuber
yields. OM may have a minimal/marginal influence here because keeping
these paths to/from it led to a better model. If we had higehr n we
might find OM significant.

Still, in the discussion, we should mention something like “further
research should focus on determining the mechanisms of causality behind
how adding organic matter could lead to the promotion of
yield-associated microbiome members”.

Might be worth it to develop a good model from another site (like MI or
ME) that might show individual taxa or fungi or both. This way we could
say that while predicting how soil treatments affect microbiomes and
yields remains a challenge, these patterns are not anomalous or
restricted to certain field sites, or bacteria vs fungi.

### Scenario 3: ME 2-yr rotation ITS in response to amendment.

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

# see which season the ASVs whose abundances we want to include were yield-associated with
its.yield.df %>% filter(ASV %in% c("ASV501", "ASV1906") & site == "ME") %>% dplyr::select(ASV, season) 
```

    ##       ASV season
    ## 1 ASV1906 spring
    ## 2  ASV501 spring
    ## 3  ASV501 summer

``` r
# they're both in spring
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

### check our premises
# is yield higher in amended?
ggplot(df3, aes(Amended, Total.yield))+geom_jitter(width = 0.1) # yeah
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
# is consortium higher in amended? 
ggplot(df3, aes(Amended, consortium))+geom_jitter(width = 0.2) # yes
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

``` r
# is OM higher in amended?
ggplot(df3, aes(Amended, OM....))+geom_jitter(width = 0.1) # maybe a little bit
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->

``` r
# convert amended treatment to 0/1
df3$Amended <- as.integer(as.logical(df3$Amended))


### tidying up the dataframe 
# remove NA cols
df3 <- df3 %>% dplyr::select(-cover_2020, -green_cov_pct_pre_kill, -Total.C...., 
                             -Total.organic.C...., -VPPG, -Root.Lesion.100.cc, -rotation_comparison)

### try a lm first
lm5 <- lm(Total.yield ~ Amended + consortium + cultivar, df3)
summary(lm5) # it's legit
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
    ## Amended       -8.876     17.752  -0.500   0.6223    
    ## consortium     9.501      4.426   2.147   0.0436 *  
    ## cultivar     -15.018     15.276  -0.983   0.3367    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 26.18 on 21 degrees of freedom
    ## Multiple R-squared:  0.5301, Adjusted R-squared:  0.4629 
    ## F-statistic: 7.896 on 3 and 21 DF,  p-value: 0.001029

``` r
extractAIC(lm5) # AIC of 166
```

    ## [1]   4.0000 166.8833

``` r
### Build SEM
# amendment did NOT increase OM or Nitrate, so we can't really work with anything here
psem1 <- psem(lm(Nitrate.N..ppm. ~ Amended, df3),
              lm(consortium ~ Nitrate.N..ppm., df3),
              lm(Total.yield ~ consortium, df3))  
coefs(psem1)
```

    ##          Response       Predictor Estimate Std.Error DF Crit.Value P.Value
    ## 1 Nitrate.N..ppm.         Amended   1.6667    1.1495 23     1.4500  0.1606
    ## 2      consortium Nitrate.N..ppm.   0.3275    0.1654 23     1.9808  0.0597
    ## 3     Total.yield      consortium  10.0205    2.1728 23     4.6117  0.0001
    ##   Std.Estimate    
    ## 1       0.2894    
    ## 2       0.3817    
    ## 3       0.6931 ***

``` r
fisherC(psem1) # model is good, individual trends are borderline
```

    ##   Fisher.C df P.Value
    ## 1   21.741  6   0.001

``` r
# what if amendment directly stimulated consortium, which increased yields?
psem2 <- psem(lm(consortium ~ Amended, df3),
              lm(Total.yield ~ consortium, df3)) 
fisherC(psem2) # valid
```

    ##   Fisher.C df P.Value
    ## 1    2.582  2   0.275

``` r
coefs(psem2) # valid and good
```

    ##      Response  Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ## 1  consortium    Amended   3.5824    0.7097 23     5.0479   0e+00       0.7250
    ## 2 Total.yield consortium  10.0205    2.1728 23     4.6117   1e-04       0.6931
    ##      
    ## 1 ***
    ## 2 ***

In Maine 2-yr rotation, Amendment increases consortium, which increases
yield. No direct path between Amended and total yield, or path through
OM or anything.

### Scenario 4: MI 3-yr rotation ITS in response to amendment.

``` r
# inspect MI amendment-yield-associated ITS ASVs
view(its.yield.df %>% filter(ASV %in% c("ASV80", "ASV161", "ASV1099", "ASV1194") & site == "MI"))
# all in spring 3-yr

# subset the phyloseq object by site, rotation, and year, and leave out controls- must compare amended/fumigated to fumigated!
ps9 <- its.ps.list[["MI.ITS.ps"]] %>% 
  subset_samples(year == 22 & rotation == 3 & !is.na(Total.yield) & general_category != "Control") 

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples
sub.df <- jim.info.s %>% filter(State == "MI" & Year == 22 & Rotation == 3)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps9 <- merge_with_jim_data(ps9, sub.df)

# define ASVs to model
mi.amd.its.spring.asvs <- c("ASV80", "ASV161", "ASV1099", "ASV1194")
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


### check our premises
# Yield is higher in amended, consortium as well
# What about any soil characteristics? not OM, solvita, N, P, or anything. maybe pH but still not significant

# convert amended treatment to 0/1
df4$Amended <- as.integer(as.logical(df4$Amended))

### tidying up the dataframe 
# remove NA cols
df4 <- df4 %>% dplyr::select(-cover_2019, -cover_2020, -green_cov_pct_pre_kill, -Total.C...., 
                             -Total.organic.C...., -VPPG, -Root.Lesion.100.cc, -rotation_comparison)

# encode cultivar as a dummy variable
df4$cultivar <- ifelse(df4$cultivar == "Superior", 0,
                       ifelse(df4$cultivar == "Burbank", 1, 2))

# which soil chemical constituents are different across treatments?
df4 %>% select_if(is.numeric) %>% 
  select(pH, Nitrate.N..ppm., NH4.N..ppm., P.Olsen..ppm., P.Bray..ppm., OM...., Solvita..ppm., POX.C..ppm.) %>% 
  map_df(~ broom::tidy(t.test(. ~df4$general_category)), .id = "var") %>% 
  filter(p.value < 0.05) %>% 
  pull(var) # none
```

    ## character(0)

``` r
# linear modeling
lm6 <- lm(Total.yield ~ consortium + Amended + cultivar, df4)
summary(lm6)
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ consortium + Amended + cultivar, data = df4)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -64.918 -34.903  -5.779  25.577  94.334 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  234.149     28.323   8.267 1.03e-07 ***
    ## consortium     8.816      7.761   1.136    0.270    
    ## Amended        5.603     58.542   0.096    0.925    
    ## cultivar      13.500     19.938   0.677    0.507    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 48.3 on 19 degrees of freedom
    ## Multiple R-squared:  0.2737, Adjusted R-squared:  0.159 
    ## F-statistic: 2.386 on 3 and 19 DF,  p-value: 0.1011

``` r
### Build SEM
# amendment directly stimulates consortium, which increases yields
psem1 <- psem(lm(consortium ~ Amended + cultivar, df4),
              lm(Total.yield ~ consortium + cultivar, df4)) 
fisherC(psem1) # valid
```

    ##   Fisher.C df P.Value
    ## 1    0.156  2   0.925

``` r
coefs(psem1) # amendment does increase consortium in summer, but not spring
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
# but amendment itself DOES stimulate yields. 
summary(lm(consortium ~ Amended, df4))
```

    ## 
    ## Call:
    ## lm(formula = consortium ~ Amended, data = df4)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0635 -1.1664 -0.4388  1.1599  3.5627 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -0.6844     0.4645  -1.474    0.155    
    ## Amended       5.0902     0.7044   7.226 4.04e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.675 on 21 degrees of freedom
    ## Multiple R-squared:  0.7132, Adjusted R-squared:  0.6995 
    ## F-statistic: 52.22 on 1 and 21 DF,  p-value: 4.041e-07

``` r
# is this consortium even enriched in amendment here? 
# plot over treatment and time

ps10 <- its.ps.list[["MI.ITS.ps"]] %>% 
  subset_samples(year == 22 & rotation == 3 & !is.na(Total.yield) & general_category == "Amended/Fumigated")

ps10@sam_data$year_cat <- paste(ps10@sam_data$year, ps10@sam_data$general_category, sep = "_")

plot.biomarkers(ps10, mi.amd.its.spring.asvs, "ASV")+
  facet_grid(~season, scales = "free", space = "free")+
  scale_fill_discrete()+
  theme_bw()+ggtitle("MI 3-yr amended/yield euks")+
  theme(axis.text.x = element_blank())
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

In MI 3-yr ITS spring, the four ASVs together didn’t make a valid model,
as amendment didn’t even increase the consortium. Plotting their %
abundances shows they all increased from spring to summer. Two of them
also were summer-associated, so I went back to try and model those.

#### Other scenarios

Here is where I’m checking for trends more exhaustively across sites,
treatments, and rotations.  
Not interested in putting any soil data here, I just want to know, in
either spring or summer of a growing season, was there evidence for a
treatment increasing any members of the microbiome that then influenced
yield? And if so, is there a direct link between treatment and yield, or
is it only through the microbiome???

``` r
# subset the phyloseq object by site, rotation, year, and season, and drop ASVs and transform counts
ps <- bact.ps.list[["OR.16S.ps"]] %>% 
  subset_samples(year == 22 & rotation == 3 & !is.na(Total.yield) &
                   season == "Spring" & general_category != "Fumigated") %>%  
  drop_ghost_asvs() %>%
  subset_occupancy(0.5) %>% 
  transform_clr() 
```

    ## 5378 features are kept. 
    ## Counts transformed with CLR.

``` r
# define ASVs to model
asvs.to.model <- c("ASV270", "ASV2762", "ASV3009", "ASV3745", "ASV574")
# asvs.to.model <- "ASV1705"

# add in the abundance of the spring associated ASVs to the sample data 
df <- cbind(data.frame(ps@sam_data), ps@otu_table[,asvs.to.model])
# df$consortium <- rowSums(df[,startsWith(colnames(df), "ASV")])

# make a consortium column of the sum of the ASVs
df$consortium <- rowSums(df[,startsWith(colnames(df), "ASV")])

# convert amended treatment to 0/1
df$Fumigated <- as.integer(as.logical(df$Fumigated))
df$Amended <- as.integer(as.logical(df$Amended))
df$Brassica <- as.integer(as.logical(df$Brassica))

df <- df %>% 
  filter(cultivar == "Norkotah") %>% 
  dplyr::select(-cover_2019, -cover_2020, -VPPG, -Root.Lesion.100.cc, 
                           -rotation_comparison, -green_cov_pct_pre_kill)

### build lm
summary(lm(Total.yield ~ consortium + Brassica, df))
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ consortium + Brassica, data = df)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -70.71 -22.24 -12.86  35.09  62.52 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  558.033     20.393  27.364 1.57e-07 ***
    ## consortium   -35.970      8.797  -4.089  0.00644 ** 
    ## Brassica     257.890     89.018   2.897  0.02744 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 49.44 on 6 degrees of freedom
    ## Multiple R-squared:  0.7823, Adjusted R-squared:  0.7097 
    ## F-statistic: 10.78 on 2 and 6 DF,  p-value: 0.01032

#### Other scenarios but I do want to examine for soil data too

**OR, mustard treatment, 3-yr rotation, bacterial ASVs**

``` r
# begin by subsetting the phyloseq object by site, rotation, and year
ps <- bact.ps.list[["OR.16S.ps"]] %>% 
  subset_samples(year == 22 & rotation == 3 & !is.na(Total.yield) & general_category != "Fumigated") 

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples
sub.df <- jim.info.s %>% filter(State == "OR" & Year == 22 & Rotation == 3)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps <- merge_with_jim_data(ps, sub.df)

# see which season the ASVs whose abundances we want to include were yield-associated with
bact.yield.df %>% filter(ASV %in% c("ASV270", "ASV2762", "ASV3009", "ASV3745", "ASV574") & site == "OR") %>% dplyr::select(ASV, season) 
```

    ##       ASV season
    ## 1  ASV270 spring
    ## 2 ASV2762 spring
    ## 3 ASV3009 spring
    ## 4 ASV3745 spring
    ## 5  ASV574 spring
    ## 6  ASV574 summer

``` r
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
df5 <- cbind(data.frame(ps@sam_data), ps@otu_table[,or.must.bact.spring.asvs])

# make a consortium column of the sum of the ASVs
df5$consortium <- rowSums(df5[,startsWith(colnames(df5), "ASV")])

# Cultivar alert: most are Norkotahs but some Burbank controls

# which chem differs?
df5 %>% select_if(is.numeric) %>% 
  select(pH, Nitrate.N..ppm., NH4.N..ppm., P.Olsen..ppm., P.Bray..ppm., OM...., Solvita..ppm., POX.C..ppm.) %>% 
  map_df(~ broom::tidy(t.test(. ~df5$Brassica)), .id = "var") %>% 
  filter(p.value < 0.05) %>% 
  pull(var) # pH, Nitrate N 
```

    ## [1] "pH"              "Nitrate.N..ppm."

``` r
# drop all-NA cols  
df5 <- df5[, colSums(is.na(df5)) != nrow(df5)]

# drop burbank cultivars and convert Brassica column from logical to numeric
df5.nork <- df5 %>% filter(cultivar == "Norkotah")
df5.nork$Brassica <- as.numeric(as.logical(df5.nork$Brassica))

summary(lm(consortium ~ Brassica, df5.nork)) # brassica strongly inc. consortium, dec. pH, and inc. nitrate
```

    ## 
    ## Call:
    ## lm(formula = consortium ~ Brassica, data = df5.nork)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4572 -1.7521 -0.0642  0.9722  3.9818 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.3317     0.8672   0.382 0.713443    
    ## Brassica      9.3067     1.5021   6.196 0.000447 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.124 on 7 degrees of freedom
    ## Multiple R-squared:  0.8458, Adjusted R-squared:  0.8237 
    ## F-statistic: 38.39 on 1 and 7 DF,  p-value: 0.000447

``` r
summary(lm(consortium ~ Nitrate.N..ppm., df5.nork)) # nitrate, but not pH, increases consortium
```

    ## 
    ## Call:
    ## lm(formula = consortium ~ Nitrate.N..ppm., data = df5.nork)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.8661 -1.8568 -0.3769  0.7936  4.9449 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)     -10.3111     4.0977  -2.516   0.0400 *
    ## Nitrate.N..ppm.   0.9130     0.2623   3.480   0.0103 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.274 on 7 degrees of freedom
    ## Multiple R-squared:  0.6337, Adjusted R-squared:  0.5814 
    ## F-statistic: 12.11 on 1 and 7 DF,  p-value: 0.01027

``` r
summary(lm(Total.yield ~ consortium, df5.nork)) # consortium, but NOT Brassica, nitrate, or pH, affects yields (negatively)
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ consortium, data = df5.nork)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -77.80 -65.62 -11.36  50.77 109.41 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  563.516     29.116   19.35 2.45e-07 ***
    ## consortium   -12.533      4.953   -2.53   0.0392 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 70.89 on 7 degrees of freedom
    ## Multiple R-squared:  0.4777, Adjusted R-squared:  0.4031 
    ## F-statistic: 6.402 on 1 and 7 DF,  p-value: 0.03922

``` r
psem5 <- psem(lm(consortium ~ Brassica + pH, df5.nork),
     lm(pH ~ Brassica, df5.nork),
     lm(Total.yield ~ consortium, df5.nork))
coefs(psem5)
```

    ##      Response  Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate
    ## 1  consortium   Brassica  11.2088    1.9454  6     5.7618  0.0012       1.1076
    ## 2  consortium         pH   1.8710    1.3233  6     1.4139  0.2071       0.2718
    ## 3          pH   Brassica  -1.0167    0.4013  7    -2.5332  0.0391      -0.6916
    ## 4 Total.yield consortium -12.5332    4.9535  7    -2.5302  0.0392      -0.6911
    ##     
    ## 1 **
    ## 2   
    ## 3  *
    ## 4  *

``` r
fisherC(psem5) # valid
```

    ##   Fisher.C df P.Value
    ## 1    7.803  4   0.099

``` r
# plot(psem5)
summary(psem5)$AIC
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

    ##       AIC    AICc  K n
    ## 1 167.316 186.916 10 9

``` r
# Mustard effects on yields
or.3yr.df <- jim.info %>% 
  filter(State == "OR" & Rotation == 3 & Year == 22 & Cultivar == "Russet Norkotah") %>% 
  dplyr::select(State, Objective, Rotation, Plot, Year, Month, `Treatment #`, `Total yield.1`)
for (i in 6:8){or.3yr.df[,i] <- as.numeric(or.3yr.df[,i])} # convert columns that should be numeric into numeric

# treatments 
or.3yr.df <- or.3yr.df %>% 
  mutate(trt = case_when(`Treatment #` == 2 ~ "control",
                         `Treatment #` == 4 ~ "mustard",
                         `Treatment #` == 5 ~ "compost",
                         `Treatment #` == 6 ~ "compost_mustard"))

# pull yields of mustard and non-mustard treatments
mustard.ylds <- or.3yr.df %>% filter(Month == 6 & str_detect(trt, "mustard")) %>% pull(`Total yield.1`)
non.mustard.ylds <- or.3yr.df %>% filter(Month == 6 & !str_detect(trt, "mustard")) %>% pull(`Total yield.1`)

t.test(mustard.ylds, non.mustard.ylds) # they're actually NOT different. 
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  mustard.ylds and non.mustard.ylds
    ## t = -1.4054, df = 11.379, p-value = 0.1866
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -152.49556   33.34613
    ## sample estimates:
    ## mean of x mean of y 
    ##  510.8051  570.3799

``` r
mean(mustard.ylds) / mean(non.mustard.ylds)
```

    ## [1] 0.8955526

The other (non-amendment) scenarios:  
OR, mustard, 3-yr: Brassica increases nitrate, decreases pH, and
increases consortium; consortium decreases yields. Best SEM included pH,
but not nitrate.

#### OR, fumigation, 3-yr.

No soil chemical data changed with treatment, so same as before.

``` r
# begin by subsetting the phyloseq object by site, rotation, and year
ps <- its.ps.list[["OR.ITS.ps"]] %>% 
  subset_samples(year == 22 & rotation == 3 & !is.na(Total.yield) & general_category != "Amended")

# subset soil chemical data, and transfer spring measurements of OM % to their corresponding summer samples
sub.df <- jim.info.s %>% filter(State == "OR" & Year == 22 & Rotation == 3)
sub.df <- sub.df %>% group_by(Plot) %>% fill(`OM (%)`)

# import this soil data into the phyloseq object
ps <- merge_with_jim_data(ps, sub.df)

# keep only the spring samples, drop, and transform
ps <- ps %>% subset_samples(season == "Summer") %>% 
  drop_ghost_asvs() %>%
  subset_occupancy(0.5) %>% 
  transform_clr() 
```

    ## 703 features are kept. 
    ## Counts transformed with CLR.

``` r
# add in the abundance of the spring associated ASVs to the sample data 
df6 <- cbind(data.frame(ps@sam_data), ps@otu_table[,"ASV2495"])

# modeling
summary(lm(ASV2495 ~ Fumigated + cultivar, df6))# ASV is affected by fumigation
```

    ## 
    ## Call:
    ## lm(formula = ASV2495 ~ Fumigated + cultivar, data = df6)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.25582 -0.45990  0.04555  0.48199  0.98498 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)       -0.7516     0.3898  -1.928  0.07783 . 
    ## FumigatedTRUE      1.8031     0.5278   3.416  0.00511 **
    ## cultivarNorkotah  -0.2334     0.4774  -0.489  0.63378   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7796 on 12 degrees of freedom
    ## Multiple R-squared:  0.4999, Adjusted R-squared:  0.4165 
    ## F-statistic: 5.997 on 2 and 12 DF,  p-value: 0.01565

``` r
summary(lm(Total.yield ~ ASV2495 + Fumigated + cultivar, df6)) # ASV (but not fumigation) affects yield
```

    ## 
    ## Call:
    ## lm(formula = Total.yield ~ ASV2495 + Fumigated + cultivar, data = df6)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -93.80 -50.03  18.04  44.05 103.17 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        527.93      39.23  13.456 3.55e-08 ***
    ## ASV2495             74.12      25.39   2.919   0.0140 *  
    ## FumigatedTRUE      -88.03      65.19  -1.350   0.2040    
    ## cultivarNorkotah   115.45      42.40   2.723   0.0198 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 68.56 on 11 degrees of freedom
    ## Multiple R-squared:  0.6088, Adjusted R-squared:  0.5022 
    ## F-statistic: 5.707 on 3 and 11 DF,  p-value: 0.01321

ND, fumigation, 3-yr: Revisited, no model found.

## A figure showing the breakdown of scenarios across site, rotation, and treatment

``` r
# import the manually curated csv file showing the scenarios (treatment x site x rotation) where yields increased/decreased, and whether members of the microbiome could be modeled as influencing this 
sem_summ <- read_csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/SEM_scenario_exhaustive_breakdown.csv")
```

    ## Rows: 33 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (6): site, treatment, trt_yield_notes, model_notes, interpretation, label
    ## dbl (6): rotation, did_trt_inc_yield (relative to control, these values are ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# reorder the labeling
sem_summ$label <- factor(sem_summ$label, levels = c("treatment decreased yields through microbiome", 
                                                     "treatment did not increase yield", 
                                                     "treatment/microbiome had partial influence on yield", 
                                                     "treatment increased yield irrespective of microbiome", 
                                                     "treatment increased yields through microbiome"))

sem_summ$site <- factor(sem_summ$site, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1"))

# barplots showing counts of scenarios, by site and by treatment
by.site.gg <- ggplot(sem_summ, aes(site, fill = label))+
  geom_bar()+
  scale_x_discrete("")+scale_y_continuous("Counts")+
  scale_fill_manual("", values = c('#fc8d62', "gray70", '#8da0cb','#e78ac3', '#66c2a5'))+
  theme_bw()

by.treatment.gg <- ggplot(sem_summ, aes(treatment, fill = label))+
  geom_bar()+
  scale_x_discrete("")+scale_y_continuous("Counts")+
  scale_fill_manual("", values = c('#fc8d62', "gray70", '#8da0cb','#e78ac3', '#66c2a5'))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# both plots side by side
by.site.gg + by.treatment.gg + plot_layout(guides = "collect", widths = c(8,3))
```

![](45_SEM_stuff_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

## Conclusions

Here I made basic SEMs from some of the 2022 field
sites/seasons/rotations where we saw interesting results (generally
these were highest numbers of treatment-yield ASVs, which were
associated with organic amendments).

I can show that soil treatments affect yields by mediating microbiome
members: across sites, rotations/seasons, bacteria or eukaryotes, and
individual ASVs or consortia. Plausible models could be made even when
treatments didn’t statistically affect yields.

Organic amendments increase abundances of more yield-associated ASVs
than fumigation or brassica treatments, although this trend is very
idiosyncratic and varies by site, rotation, and season.

Soil chemistry (OM%, pH) can sometimes provide explanation for how
yield-associated taxa respond to treatments.

Higher n in field sites/rotations of interest would allow us to better
model causal relationships to better understand mechanisms of how
microbiome members respond to soil treatments and affect yields.

Later, exhaustive efforts to model how treatments influenced ASVs that
influenced yields showed a few interesting patterns. Where we observed
treatments influencing ASVs that influenced yields (in western and
eastern sites), amendment-associated ASVs increased yields while
fumigation-associated ASVs decreased them in both Minnesota sites.

#### session info

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.7.3
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
    ##  [13] fastmap_1.2.0           backports_1.5.0         XVector_0.44.0         
    ##  [16] labeling_0.4.3          utf8_1.2.4              rmarkdown_2.28         
    ##  [19] tzdb_0.4.0              nloptr_2.1.1            UCSC.utils_1.0.0       
    ##  [22] bit_4.0.5               xfun_0.47               zlibbioc_1.50.0        
    ##  [25] GenomeInfoDb_1.40.1     jsonlite_1.8.8          biomformat_1.32.0      
    ##  [28] highr_0.11              rhdf5filters_1.16.0     Rhdf5lib_1.26.0        
    ##  [31] broom_1.0.6             parallel_4.4.1          cluster_2.1.6          
    ##  [34] R6_2.5.1                stringi_1.8.4           RColorBrewer_1.1-3     
    ##  [37] boot_1.3-30             car_3.1-2               estimability_1.5.1     
    ##  [40] Rcpp_1.0.13             iterators_1.0.14        knitr_1.48             
    ##  [43] zoo_1.8-12              IRanges_2.38.1          splines_4.4.1          
    ##  [46] igraph_2.0.3            timechange_0.3.0        tidyselect_1.2.1       
    ##  [49] rstudioapi_0.16.0       abind_1.4-5             yaml_2.3.10            
    ##  [52] MuMIn_1.48.4            vegan_2.6-8             codetools_0.2-20       
    ##  [55] lattice_0.22-6          plyr_1.8.9              Biobase_2.64.0         
    ##  [58] withr_3.0.1             evaluate_0.24.0         survival_3.6-4         
    ##  [61] bayesm_3.1-6            Biostrings_2.72.1       pillar_1.9.0           
    ##  [64] carData_3.0-5           DiagrammeR_1.0.11       tensorA_0.36.2.1       
    ##  [67] foreach_1.5.2           stats4_4.4.1            insight_0.20.4         
    ##  [70] generics_0.1.3          vroom_1.6.5             S4Vectors_0.42.1       
    ##  [73] hms_1.1.3               munsell_0.5.1           scales_1.3.0           
    ##  [76] minqa_1.2.8             xtable_1.8-4            glue_1.7.0             
    ##  [79] emmeans_1.10.4          tools_4.4.1             robustbase_0.99-4      
    ##  [82] data.table_1.16.0       visNetwork_2.1.2        mvtnorm_1.3-1          
    ##  [85] rhdf5_2.48.0            grid_4.4.1              ape_5.8                
    ##  [88] colorspace_2.1-1        nlme_3.1-164            GenomeInfoDbData_1.2.12
    ##  [91] performance_0.12.3      cli_3.6.3               fansi_1.0.6            
    ##  [94] gtable_0.3.5            DEoptimR_1.1-3          digest_0.6.37          
    ##  [97] BiocGenerics_0.50.0     TH.data_1.1-2           htmlwidgets_1.6.4      
    ## [100] farver_2.1.2            htmltools_0.5.8.1       multtest_2.60.0        
    ## [103] lifecycle_1.0.4         httr_1.4.7              bit64_4.0.5            
    ## [106] MASS_7.3-60.2
