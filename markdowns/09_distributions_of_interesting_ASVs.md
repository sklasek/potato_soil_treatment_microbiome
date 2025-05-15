distributions of interesting ASVs
================
Scott Klasek
07 May, 2025

## Purpose

Now we have a somewhat-complete understanding of which ASVs in which
sites were enriched due to soil treatments and their affects on yields.
In some cases, we can see how they are influenced by soil chemistry.
Step 1: visualize these treatment scenarios as sankey plots.

Another final question I’d like to resolve is whether these “special”
ASVs (I really need to think of a term for them) are most highly
distributed in the sites in which they were treatment/yield-associated.
This will resolve whether this is 1) a biogeography issue (ie, if only
we had these magic ASVs in Minnesota, or even just higher relative
abundances of them, then we’d see similar yield-increases with
amendments as we saw in, say, Oregon) or if it’s 2) a supremely complex
issue of soil chemistry and microbiomes and a bunch of other variables
all interacting together that make responsiveness to treatments and
affects on yields so variable across taxonomy, space, and time. HINT: I
STRONGLY, STRONGLY suspect it’s \#2, but we need to convince people.

## Setup

#### load libraries

``` r
packages <- c("tidyverse", "phyloseq", "speedyseq", "patchwork", "networkD3", "ggalluvial", "NatParksPalettes")
invisible(lapply(packages, require, character.only = TRUE))
```

    ## Loading required package: tidyverse

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4     
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
    ## Loading required package: networkD3

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'networkD3'

    ## Loading required package: ggalluvial
    ## Loading required package: NatParksPalettes

#### load data

``` r
### phyloseq objects across all sites
its.all.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/all_obj1_by_site/all.ITS.ps")

# trim down bacterial phyloseq to only the 2022 samples
bact.all.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/all_obj1_by_site/all.16S.ps") 
# bact.22.ps <- bact.all.ps %>% subset_samples(year == 22)
# saveRDS(bact.22.ps, file = "/Users/klas0061/Desktop/bact.22.ps")
bact.22.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/all_obj1_by_site/bact.22.ps")

# subset giant bacterial phyloseq by state, year, and treatment
oreida.19.22.ps <- bact.all.ps %>% subset_samples(state %in% c("ID", "OR") &
                                               year %in% c(19, 22) & general_category != "Fumigated")

# separate subsets of OR and ID bacterial phyloseqs, by rotation
id.select.ps <- bact.all.ps %>% subset_samples(state == "ID" & year %in% c(19, 22) & rotation == 3 & general_category != "Fumigated")
or.select.ps <- bact.all.ps %>% subset_samples(state == "OR" & year %in% c(20, 22) & rotation == 2 & general_category != "Fumigated")

rm(bact.all.ps)
gc()
```

    ##              used   (Mb) gc trigger    (Mb) limit (Mb)   max used    (Mb)
    ## Ncells    5831556  311.5    8437255   450.6         NA    6708171   358.3
    ## Vcells 1195208724 9118.8 3512283555 26796.6      32768 3511650201 26791.8

``` r
### treatment and yield associated ASVs
its.asvs <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/its.asvs.treatment.and.yield.associated.csv")
bact.asvs <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/bact.asvs.treatment.and.yield.associated.csv")
all.trt.yld.asvs <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/all.asvs.treatment.and.yield.associated.filtered.csv")

### just yield-associated ASVs with rotation/season info for context
its.yield.df <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/its.yield.lm.results.txt")
bact.yield.df <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/bact.yield.lm.results.txt")
bact.yield.all.df <- read.csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/bact.yield.lm.results.UNFILTERED.txt")

# manually curated csv for alluvial plot
alluv <- read_csv(file = "/Users/klas0061/Desktop/UMN/treatment_variance_modeling/alluvial_modeling_info.csv")
```

    ## Rows: 33 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (4): site, treatment, did treatment increase yield? (relative to control...
    ## dbl (2): rotation, were there treatment- and yield-associated ASVs found in ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

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
    facet_grid(~general_category, scales = "free", space = "free")

  return(bm.barplot) # return the plot
}


plot.biomarkers.2 <- function(ps, biomarkerlist, taxlevel){
  
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
  
  return(ps.bm) # return the phyloseq instead of a plot, so you can tailor it 
}
```

## Sankey plots

### With SankeyNetwork

#### all treatments combined

``` r
# make it from scratch... it's very fiddly
source <- c(0, 0, 0, 1, 1, 2, 2, 3, 4, 4, 5)
target <- c(1, 2, 3, 4, 5, 4, 5, 5, 6, 7, 7)
value <- c(7, 25, 1, 3, 4, 13, 12, 1, 7, 9, 17)
nodes <- data.frame("name" = c("treatment", "yield increase", "no yield change", "yield decrease", "ASVs", 
                               "no ASVs", "microbiome model", "no microbiome model"))

sank <- data.frame(source, target, value)
colors <- 'd3.scaleOrdinal() .domain(["treatment", "yield increase", "no yield change", "yield decrease", "ASVs", 
                               "no ASVs", "microbiome model", "no microbiome model"]) .range(["blue", "black", "purple", "yellow", "red", "yellow", "orange", "brown"])'

sank.plot <- sankeyNetwork(Links = sank, Nodes = nodes,
 Source = "source", Target = "target",
 Value = "value", NodeID = "name",
 colourScale=colors,
 fontSize= 12, nodeWidth = 30)

sank.plot
```

#### amendment treatments only

``` r
source <- c(0, 0, 0, 1, 1, 2, 2, 3, 4, 4, 5)
target <- c(1, 2, 3, 4, 5, 4, 5, 5, 6, 7, 7)
value <- c(4, 9, 1, 3, 1, 5, 4, 1, 4, 4, 6)
nodes <- data.frame("name" = c("treatment", "yield increase", "no yield change", "yield decrease", "ASVs", 
                               "no ASVs", "microbiome model", "no microbiome model"))

sank <- data.frame(source, target, value)

sank.amd <- sankeyNetwork(Links = sank, Nodes = nodes,
 Source = "source", Target = "target",
 Value = "value", NodeID = "name",
 fontSize= 12, nodeWidth = 25)

sank.amd
```

#### fumigated treatments only

``` r
source <- c(0, 0, 1, 2, 2, 3, 3, 4)
target <- c(1, 2, 4, 3, 4, 5, 6, 6)
value <- c(3, 10, 3, 6, 4, 2, 4, 7)
nodes <- data.frame("name" = c("treatment", "yield increase", "no yield change", "ASVs", 
                               "no ASVs", "microbiome model", "no microbiome model"))

sank <- data.frame(source, target, value)

sank.fum <- sankeyNetwork(Links = sank, Nodes = nodes,
 Source = "source", Target = "target",
 Value = "value", NodeID = "name",
 fontSize= 12, nodeWidth = 25)

sank.fum
```

#### mustard treatments only

``` r
source <- c(0, 1, 1, 2, 2, 3)
target <- c(1, 2, 3, 4, 5, 5)
value <- c(6, 2, 4, 1, 1, 4)
nodes <- data.frame("name" = c("treatment", "no yield change", "ASVs", 
                               "no ASVs", "microbiome model", "no microbiome model"))

sank <- data.frame(source, target, value)

sank.must <- sankeyNetwork(Links = sank, Nodes = nodes,
 Source = "source", Target = "target",
 Value = "value", NodeID = "name",
 fontSize= 12, nodeWidth = 25)

sank.must
```

### With ggalluvial package

SankeyNetwork is fine except you can’t knit these plots, and horrible
things happen when you try to specify colors.  
Also, ggalluvial lets you work with data in tabular form without having
to count and summarize each node.

``` r
# make the column names a bit shorter
colnames(alluv) <- c("site", "treatment", "rotation",
                     "trt_inc_yield", "trt_yld_asvs_found", "microbiome_model")

# convert 0s and 1s into True/False
alluv$trt_yld_asvs_found <- as.logical(as.numeric(alluv$trt_yld_asvs_found))

# make site label and color scheme
alluv$site <- factor(alluv$site, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1"))
cols.order <- c("#DC4405", "#B3A369", "#236192", "#7A0019", "#FFC72A", "#C5050C", "#18453B", "#B0D7FF")

# make plots
rot2.alluv <- ggplot(alluv %>% filter(rotation == 2), aes(axis1 = trt_inc_yield, axis2 = trt_yld_asvs_found, axis3 = microbiome_model))+
  geom_alluvium(aes(fill = site), alpha = 0.6)+
  scale_fill_manual("Field Site", values = cols.order)+
  geom_stratum(fill = "gray70", width = 0.5)+
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  facet_grid(treatment~., scales = "free", space = "free", switch = "y")+
  scale_y_continuous("", breaks = seq(0,7,1))+
  scale_x_discrete(limits = c("Treatment effect on yield", "Target ASVs found", "Microbiome-informed model"), expand = c(.05, .05)) +
  ggtitle("2-year rotation")+theme_classic()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

rot3.alluv <- ggplot(alluv %>% filter(rotation == 3), aes(axis1 = trt_inc_yield, axis2 = trt_yld_asvs_found, axis3 = microbiome_model))+
  geom_alluvium(aes(fill = site), alpha = 0.6)+
  scale_fill_manual("Field Site", values = cols.order)+
  geom_stratum(fill = "gray70", width = 0.5)+
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3)+
  facet_grid(treatment~., scales = "free", space = "free", switch = "y")+
  scale_y_continuous("", breaks = seq(0,7,1))+
  scale_x_discrete(limits = c("Treatment effect on yield", "Target ASVs found", "Effect of target ASVs on yields\n in microbiome model"), expand = c(.05, .05)) +
  ggtitle("3-year rotation")+theme_classic()

sep.by.rotation.gg <- rot2.alluv / rot3.alluv + plot_layout(guides = "collect")
sep.by.rotation.gg
```

    ## Warning in to_lodes_form(data = data, axes = axis_ind, discern =
    ## params$discern): Some strata appear at multiple axes.
    ## Warning in to_lodes_form(data = data, axes = axis_ind, discern =
    ## params$discern): Some strata appear at multiple axes.
    ## Warning in to_lodes_form(data = data, axes = axis_ind, discern =
    ## params$discern): Some strata appear at multiple axes.
    ## Warning in to_lodes_form(data = data, axes = axis_ind, discern =
    ## params$discern): Some strata appear at multiple axes.
    ## Warning in to_lodes_form(data = data, axes = axis_ind, discern =
    ## params$discern): Some strata appear at multiple axes.
    ## Warning in to_lodes_form(data = data, axes = axis_ind, discern =
    ## params$discern): Some strata appear at multiple axes.

``` r
grid::grid.draw(grid::textGrob("Scenario count", x = 0.02, rot = 90))
```

![](46_distributions_of_interesting_ASVs_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Trying to make an alluvial plot with NA values omitted… it’s harder than
it should be

``` r
# a reproducible example
data(majors)
majors$curriculum <- as.factor(majors$curriculum)
ggplot(majors,
   aes(x = semester, stratum = curriculum, alluvium = student,
          fill = curriculum, label = student)) + #changed from label = alluvium
   scale_fill_brewer(type = "qual", palette = "Set2") +
   geom_flow(stat = "alluvium", lode.guidance = "frontback",
           color = "darkgray") + #can change lode.guidance parameter here in geom_flow
   geom_stratum() +
   geom_text(stat = "alluvium", size = 3)


### try to do this with my data

# subset df
testdf <- alluv %>% filter(rotation == 3) 

# convert logical to character
testdf$trt_yld_asvs_found <- as.character(as.logical(testdf$trt_yld_asvs_found))

# pivot longer
testdf <- testdf %>% pivot_longer(-c("site", "treatment", "rotation"), values_to = "label", names_to = "step")

# concatenate site, treatment, and rotationk, then specify THAT charcter vector as an alluvium in the aes below? Is this what they mean by lodes form?
testdf$scenario <- paste(testdf$site, testdf$treatment, testdf$rotation, sep = "_")

testdf <- testdf[sample(nrow(testdf)),]

ggplot(testdf,
       aes(x = factor(step, levels = c("trt_inc_yield", "trt_yld_asvs_found", "microbiome_model")), 
           stratum = label, alluvium = scenario, fill = site, label = label))+
  geom_flow(stat = "alluvium")+
  geom_stratum(aes(fill = site), na.rm = TRUE)+
  geom_label(stat = "stratum")+
  scale_fill_manual("Field Site", values = cols.order)+
  facet_grid(treatment~., scales = "free", space = "free", switch = "y")+
  scale_x_discrete("")+
  ggtitle("3-year rotation")+theme_classic()

    
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3)+

  scale_x_discrete(limits = c("Treatment effect on yield", "Target ASVs found", "Effect of target ASVs on yields\n in microbiome model"), expand = c(.05, .05)) +
  ggtitle("3-year rotation")+theme_classic()
```

Most of the microbiome changes that were associated with increased
yields were from amended treatments. Also interesting to note that
treatment effects on yield didn’t need to be statistically significant
to detect meaningful links between treatments –\> ASVs –\> yields.

## Plot distributions of the ASVs (in 2022)

#### Michigan OM and a template for ITS

``` r
# subset phyloseq
its.22.ps <- its.all.ps %>% 
  subset_samples(year == 22)

# if only a few ASVs, just define them manually
mi.its.amd.asvs <- c("ASV80", "ASV161", "ASV1099", "ASV1194")
me.its.amd.asvs <- c("ASV501", "ASV1906", "ASV1642")

# plot over treatment and time
mi.its.gg <- plot.biomarkers(its.22.ps, mi.its.amd.asvs, "ASV")+
  facet_grid(~state, scales = "free", space = "free")+
  scale_fill_discrete()+
  theme_bw()+ggtitle("MI OM-influenced euks across all sites")+
  theme(axis.text.x = element_blank())
# from this plot I could see that MI ASVs were also present in OR and MN2. ME ASVs also only present in MN1.

# re-subset phyloseq and make a label
ps1 <- its.22.ps %>% subset_samples(state %in% c("ND", "OR", "MI"))
ps1@sam_data$site_cat <- paste(ps1@sam_data$state, ps1@sam_data$general_category, sep = "_")

# plot them
plot.biomarkers(ps1, mi.its.amd.asvs, "ASV")+
    facet_grid(~site_cat, scales = "free", space = "free")+
  scale_fill_discrete()+
  theme_bw()+ggtitle("MI Amendment-influenced euks across other sites they are abundant in")+
  theme(axis.text.x = element_blank())
```

![](46_distributions_of_interesting_ASVs_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Check their yield-associations:

``` r
its.yield.df %>% filter(ASV %in% c("ASV80", "ASV161")) %>% 
  dplyr::select(ASV, lm_mult_r_sq, lm_adj_r_sq, estimate, p, site, closest_tax)
```

    ##      ASV lm_mult_r_sq lm_adj_r_sq  estimate            p site
    ## 1 ASV161    0.2140044   0.1848935  18.94476 0.0003831722   MI
    ## 2  ASV80    0.2104671   0.1812251  16.75003 0.0001467432   MI
    ## 3  ASV80    0.1765818   0.1391537 -44.54685 0.0211386766   OR
    ##                  closest_tax
    ## 1          Kernia columnaris
    ## 2 Botryotrichum spirotrichum
    ## 3 Botryotrichum spirotrichum

MI amendment euks (calling them euks, not fungi because one is a
heterolobosa) are also highly abundant in MN2 and OR, and two are
particularly responsive to amendment in MN2 but they are NOT
yield-associated. (Actually if anything, ASV80 is borderline negatively
associated with yields in OR).

Maine: Of the two amendment-yield ITS ASVs, only ASV501 was present in
one other site: MN. It was somewhat enriched in amended plots, but not
strongly yield-associated.

#### ID amend and a template for 16S

``` r
# get the asvs
id.bact.amd.asvs <- bact.asvs %>% filter(site == "ID" & treatment == "Amended") %>% pull(ASV)

# filter by rotation and season
id.bact.amd.asvs <- bact.yield.df %>% filter(ASV %in% id.bact.amd.asvs & 
                                               site == "ID" & rotation == 3 & season == "summer") %>% pull(ASV)

# plot over treatment and time
oreida.bact.ps <- bact.22.ps %>% subset_samples(state %in% c("ID", "OR"))
oreida.bact.ps@sam_data$site_cat <- paste(oreida.bact.ps@sam_data$state, 
                                          oreida.bact.ps@sam_data$general_category, sep = "_")
oreida.bact.gg <- plot.biomarkers(oreida.bact.ps, id.bact.amd.asvs, "ASV")+
  facet_grid(~site_cat, scales = "free", space = "free")+
  scale_fill_discrete()+
  theme_bw()+ggtitle("ID OM-influenced bacteria across Ore-Ida")+
  theme(axis.text.x = element_blank())
oreida.bact.gg
```

![](46_distributions_of_interesting_ASVs_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Check their yield-associations:

``` r
bact.yield.all.df %>% filter(ASV %in% id.bact.amd.asvs & site == "OR") %>% 
  dplyr::select(ASV, lm_mult_r_sq, lm_adj_r_sq, estimate, p, site, closest_tax)
```

    ##        ASV lm_mult_r_sq lm_adj_r_sq    estimate         p site
    ## 1 ASV17047 1.222637e-01  0.04246953 -60.0426100 0.4347977   OR
    ## 2  ASV2340 4.953873e-02 -0.02966638  18.6833525 0.4244175   OR
    ## 3  ASV2340 1.942913e-02 -0.06971368  -6.0431880 0.8330254   OR
    ## 4  ASV2345 7.426491e-05 -0.07684310  -0.6000458 0.9756835   OR
    ## 5  ASV2345 7.585734e-03 -0.08263374   0.1824918 0.9943709   OR
    ## 6  ASV2589 2.550960e-02 -0.05569793  14.0021660 0.5745482   OR
    ## 7  ASV2589 1.567213e-03 -0.08919940   1.4542501 0.9453599   OR
    ## 8  ASV3781 1.637289e-03 -0.08912296  13.4553424 0.7540595   OR
    ## 9  ASV4416 2.178068e-02 -0.06714835  -8.9609731 0.7787209   OR
    ##                 closest_tax
    ## 1         Virgibacillus sp.
    ## 2            Longispora sp.
    ## 3            Longispora sp.
    ## 4 Novibacillus thermophilus
    ## 5 Novibacillus thermophilus
    ## 6      Fam. Geminicoccaceae
    ## 7      Fam. Geminicoccaceae
    ## 8    Fam. Hyphomicrobiaceae
    ## 9        Hyphomicrobium sp.

ID 3-yr amendment-yield bacterial ASVs are present in OR in nearly the
same abundances, and are likewise enriched in amended plots, but they
are not even closely associated with yields.

The ID ASVs were present in CO as well, and enriched in amendment
treatments, but patchy. Same with ME, but only three of the seven were
found. Two were barely present in MI and not associated with amendment,
and none were found in MN1, MN2, or WI.

#### inspect whether these ID ASVs really increased over time in OR.

``` r
or.bact.ps <- oreida.19.22.ps %>% subset_samples(state == "OR")
or.bact.ps@sam_data$year_cat <- paste(or.bact.ps@sam_data$year, 
                                          or.bact.ps@sam_data$general_category, sep = "_")
or.id.asvs.bact.gg <- plot.biomarkers(or.bact.ps, id.bact.amd.asvs, "ASV")+
  facet_grid(~year_cat, scales = "free", space = "free")+
  scale_fill_discrete()+
  theme_bw()+ggtitle("ID OM-influenced bacteria in Oregon across time and treatment")+
  theme(axis.text.x = element_blank())
or.id.asvs.bact.gg
```

![](46_distributions_of_interesting_ASVs_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## Plot distributions of the ASVs (any year)

#### ID 3-yr Amended consortium in ID and OR, across 2019 and 2022, by season and treatment.

Excluding non-Norkotah cultivars, 2-yr rotation samples, and fumigated
treatments.  
Thought about keeping specific amendment treatments separate, but A)
they didn’t differ, and B) don’t put too much extraneous information on
the plot if you’re not going to explain it.

``` r
#### take separately subsetted ps obects from OR years 2020, 2022 (2yr) and ID years 2019, 2022 (3yr) and then recombine them
test.ps <- merge_phyloseq(id.select.ps, or.select.ps)

# reclassify ID green manure as control
test.ps@sam_data[which(test.ps@sam_data$state == "ID" & 
                         test.ps@sam_data$amendment == "Green Manure"),"general_category"] <- "Control"

# keep only Norkotah cultivars 
test.ps <- test.ps %>% subset_samples(cultivar == "Norkotah")

# four digit year
test.ps@sam_data$year <- paste(20, test.ps@sam_data$year, sep = "")

# add label based on state, year, general category, rotation, season, amendment 
test.ps@sam_data$label <- paste(test.ps@sam_data$state, test.ps@sam_data$year, 
                                test.ps@sam_data$general_category, test.ps@sam_data$season, sep = "_")

# average all samples by label. careful here: 
# default behavior is to average every numeric value in the sample data, and to wipe out every categorical value with NA, 
# even if it was repeated for each entry in the group
test.ps <- merge_samples(test.ps, "label")
```

    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion
    ## Warning in asMethod(object): NAs introduced by coercion

``` r
# rewrite state & general_category
test.ps@sam_data$state <- str_sub(rownames(test.ps@sam_data), 1, 2)
test.ps@sam_data$general_category <- str_sub(rownames(test.ps@sam_data), 7, 13)

# concatenate state & year
test.ps@sam_data$facet_label <- paste(test.ps@sam_data$state, test.ps@sam_data$year, sep = "_")

# add taxonomy information
tax.df <- data.frame(test.ps@tax_table)
tax.df$ASV <- rownames(tax.df)
tax.df$closest_tax <- ifelse(!is.na(tax.df$Species), paste(tax.df$Genus, tax.df$Species), 
                    ifelse(!is.na(tax.df$Genus), paste(tax.df$Genus, "sp."),     
                    ifelse(!is.na(tax.df$Family), paste("Fam.", tax.df$Family),
                    ifelse(!is.na(tax.df$Order), paste("Ord.", tax.df$Order),
                    ifelse(!is.na(tax.df$Class), paste("Cl.", tax.df$Class),
                    ifelse(!is.na(tax.df$Phylum), paste("Phy.", tax.df$Phylum),
                    NA))))))
tax.df$label <- paste(tax.df$closest_tax, tax.df$ASV, sep = " ")
tax_table(test.ps) <- as.matrix(as.data.frame(tax.df))

# I gave up on trying to manipulate the row names because they are required to be unique. 
# but can I plot x-axis by something else? 
test.ps@sam_data$new_label <- substr(rownames(test.ps@sam_data), 9, nchar(test.ps@sam_data))

# plot
oreida.gg <- plot.biomarkers.2(test.ps, id.bact.amd.asvs, "label") %>% 
  plot_bar(x = "new_label", fill="label")+
    geom_bar(stat="identity", position="stack")+
    scale_y_continuous("% Abundance")+
    scale_x_discrete("")+
    facet_grid(~facet_label, scales = "free", space = "free")+
    scale_fill_manual("ASV taxonomy", values = natparks.pals("Triglav", 9))+
    theme_bw()+ggtitle("")+
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1))
```

    ## Warning in psmelt(physeq, as = "data.frame"): The sample variables: 
    ## label
    ##  have been renamed to: 
    ## sample_label
    ## to avoid conflicts with taxonomic rank names.

``` r
oreida.gg
```

![](46_distributions_of_interesting_ASVs_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Interesting- the variation WITHIN ID amendment treatments is due to
whether compost was added or not. These taxa are heavily
compost-associated ones and are far less abundant in the green
manure-only amendment. This is probably why we didn’t see a
statistically significant OM increase with the general category
“Amended”.

### Notes on distributions

#### Where are treatment-yield ASVs much higher in abundance (or only detected in) sites where they were effective?

- MN2 amended eukaryote ASV493 pretty much only in MN2.  
- OR fumigated eukaryote ASV2495 pretty much only in OR.

#### Where were these ASVs abundant and prevalent to a similar degree in other sites?

Pretty much all other cases. Here are some examples:  
- MI amended eukaryote ASVs 80, 161 were also highly abundant in MN2 and
OR.  
- ME amended eukaryote ASVs 501 also highly abundant in MN1.  
- **ID amended bacterial ASVs were also higher in OR amended samples,
but were not yield-associated anywhere but ID.** These ASVs are not in
MN1, MN2, or WI samples, only ~ two are present in low abundance in MI
(and not responsive to amendment), and only three are present in ME but
they are higher in the amended samples. This example alone shows how
continental-scale distribution, treatment, and unexplained factors
influence target ASVs.  
(In contrast, the OR 2-yr summer amended bacterial ASVs were not
different across treatments, though were mostly all present across all
field sites. So, why did they work in OR for one rotation?)

#### session info

``` r
sessionInfo()
```

    ## R version 4.4.3 (2025-02-28)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS Sonoma 14.7.4
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
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
    ##  [1] NatParksPalettes_0.2.0 ggalluvial_0.12.5      patchwork_1.3.0       
    ##  [4] speedyseq_0.5.3.9021   phyloseq_1.50.0        lubridate_1.9.4       
    ##  [7] forcats_1.0.0          stringr_1.5.1          dplyr_1.1.4           
    ## [10] purrr_1.0.4            readr_2.1.5            tidyr_1.3.1           
    ## [13] tibble_3.2.1           ggplot2_3.5.1          tidyverse_2.0.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] ade4_1.7-23             tidyselect_1.2.1        farver_2.1.2           
    ##  [4] Biostrings_2.74.1       fastmap_1.2.0           digest_0.6.36          
    ##  [7] timechange_0.3.0        lifecycle_1.0.4         cluster_2.1.8          
    ## [10] survival_3.8-3          magrittr_2.0.3          compiler_4.4.3         
    ## [13] rlang_1.1.4             tools_4.4.3             igraph_2.1.4           
    ## [16] utf8_1.2.4              yaml_2.3.8              data.table_1.17.0      
    ## [19] knitr_1.47              labeling_0.4.3          bit_4.6.0              
    ## [22] plyr_1.8.9              withr_3.0.2             BiocGenerics_0.52.0    
    ## [25] grid_4.4.3              stats4_4.4.3            fansi_1.0.6            
    ## [28] multtest_2.62.0         biomformat_1.34.0       colorspace_2.1-1       
    ## [31] Rhdf5lib_1.28.0         scales_1.3.0            iterators_1.0.14       
    ## [34] MASS_7.3-64             cli_3.6.3               rmarkdown_2.27         
    ## [37] vegan_2.6-10            crayon_1.5.3            generics_0.1.3         
    ## [40] rstudioapi_0.17.1       httr_1.4.7              reshape2_1.4.4         
    ## [43] tzdb_0.4.0              ape_5.8-1               rhdf5_2.50.2           
    ## [46] zlibbioc_1.52.0         splines_4.4.3           parallel_4.4.3         
    ## [49] XVector_0.46.0          vctrs_0.6.5             Matrix_1.7-2           
    ## [52] jsonlite_1.8.8          IRanges_2.40.1          hms_1.1.3              
    ## [55] S4Vectors_0.44.0        bit64_4.6.0-1           foreach_1.5.2          
    ## [58] glue_1.7.0              codetools_0.2-20        stringi_1.8.4          
    ## [61] gtable_0.3.6            GenomeInfoDb_1.42.3     UCSC.utils_1.2.0       
    ## [64] munsell_0.5.1           pillar_1.9.0            htmltools_0.5.8.1      
    ## [67] rhdf5filters_1.18.1     GenomeInfoDbData_1.2.13 R6_2.5.1               
    ## [70] vroom_1.6.5             evaluate_1.0.3          lattice_0.22-6         
    ## [73] Biobase_2.66.0          highr_0.11              Rcpp_1.0.14            
    ## [76] nlme_3.1-167            permute_0.9-7           mgcv_1.9-1             
    ## [79] xfun_0.45               pkgconfig_2.0.3
