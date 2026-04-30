# Load Packages ----
library(ggplot2)
library(vegan)
library(dplyr)
library(phyloseq)
library(qiime2R)
library(tidyverse)
library(lme4)
library(car)
library(bbmle)
library(lmtest)
library(ape)
library(indicspecies)
library(MASS)
library(ecole)
library(ANCOMBC)
library(emmeans)
library(mirlyn)

# Install Packages ----

# #This is how I installed phyloseq. if it doesn't work, restart R studio and do this as the first thing in the new session
# source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
#        local = TRUE)
# 
# # and this is how I installed qiime2R
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TreeSummarizedExperiment")
# 
# install.packages("remotes")
# remotes::install_github("jbisanz/qiime2R")

#to install pairwiseAdonis
# install.packages('devtools')
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

#to install ecole
#install.packages("remotes")
#remotes::install_github("phytomosaic/ecole")

#to install mirlyn
#BiocManager::install("escamero/mirlyn")


# Load Metadata ----

###### 1. Load data. Alternative: once you run this once, you can then save it as an R file
#and load it directly. the code for this is in section 2d and 3e

### 1a. Metadata and elevation and slope/aspect
#the metadata

metadata<-readr::read_tsv("18sformatted_metadata.tsv")
#this is the elevation file
waypoints<- read.csv("waypoints.csv")
#slope and aspect file
slope_aspect<- read.csv("latrine_geog_info.csv")
#vicugna RAI
vicugnaRAI<- read.csv("vicugnaRAI_20260304.csv") #using data downloaded Mar 4 2026. 30 min IE used
#chronosequences
chrono<- read.csv('soil_chronosequence_points.csv') #using the dual Tang and Seimon method decided April 28 2026

### 1b. Other files, loaded into a phyloseq
# Load Other Data ----

#load it into a phyloseq object
#it wasn't loading the sample variables correctly so I had to open the metadata as an excel sheet and 
#just re-save it as a tsv so it doesn't cut off the column names. I also added NAs in the cells for the pos/neg
#controls that were initially blank. this is just a note and only needs to be addressed if a future
#metadata file is used and chops off the header names
phy <- qza_to_phyloseq("18S_0326_qiime/Peru_18S_0326_table.qza", 
                       "18S_0326_qiime/Peru_18S_0326_rooted-tree.qza", 
                       "18S_0326_qiime/Peru_18S_0326_taxonomy.qza",
                       "18sformatted_metadata.tsv")
#check that metadata didn't chop off names using sample_variables(phy). this should be the header names, not data
sample_variables(phy)


# Initial Processing ----
####### 2. Filtering

## some samples were redone but both samples are in the data (ex 12 and 12B)
#so we need to filter those out so we aren't dealing with duplicates

#look at the total number of ASVs each sample had
asv_all<- data.frame(otu_table(phy))
asv_tot<-data.frame(colSums(asv_all))

#filter only for the samples that had duplicates to compare
asv_tot<- asv_tot %>% 
  filter(row.names(asv_tot) %in% c('X47', 'X49','X48','X50','X19','X19B','X20','X20B','X4','X4B','X23','X23B','X24','X24B','X16','X16B','X11','X11B','X12','X12B','X8','X8B','X7','X7B', 'X15','X15B','X3','X3B'))

#see samples with duplicates with latrine names 
asv_tot2 <- asv_tot
asv_tot2$SampleID <- rownames(asv_tot2)  
asv_tot2$SampleID <- sub("^X", "", asv_tot2$SampleID)
asv_tot2_meta <- asv_tot2 %>%
  left_join(metadata, by = "SampleID")


#drop the samples that had the lower number of ASVs out of that duplicate pair
dup_phy<- subset_samples(phy, !row.names(phy@sam_data) %in% c('49', '48','19B','20','4B','23','24','16','11B','12B','8','7','15','3'))


### 2a. double check heading names 
taxa_data<- as.data.frame(tax_table(dup_phy)) #pull out taxa table
unique(taxa_data$Kingdom) #check column names are intact 
tax_table(dup_phy)<- as.matrix(taxa_data) #put it back into the phyloseq


### 2b. Filter out mitochondria, chloroplasts, and non eukaryotes from phyloseq object [UPDATED filtering schema below]
#filtered_phy <- dup_phy %>%
#  subset_taxa(   #the subset keeps rows only where the following operators are met
#   Kingdom == "Eukaryota" )

#filtered_phy2<- filtered_phy %>% 
#  subset_taxa(Family != 'Craniata' | is.na(Family))

#filtered_phy3<- filtered_phy2 %>% 
#  subset_taxa(Family != 'Embryophyceae' | is.na(Family))


### 2c. Filter out singletons 
pruned_filtered_phy<-prune_taxa(taxa_sums(dup_phy)>1, dup_phy) #this is similar to subset but it keeps only the taxa that had more than 1 occurence using the taxa sums function
dup_phy
pruned_filtered_phy
#print to compare, removed some taxa

### 2d. Filter out latrines that were switched AND the pos and neg controls
#final_filtered_phy<-subset_samples(pruned_filtered_phy, 
#    !latrine_trt_month %in% c("L83_latrine_wet", 'L83_control_wet', "L94_latrine_wet", 'L94_control_wet', "L62_control_wet", 'L62_latrine_wet') #the ! means not, so we are removing those latrines
#                                    !is.na(latrine)) #this removes NAs (the pos and neg controls)
#print the two to compare
#pruned_filtered_phy
#final_filtered_phy

####### CURRENT filtering

#filter out positive and negative controls 
pruned_filtered_phy <-subset_samples(pruned_filtered_phy, !is.na(latrine)) 

#filter out NA at the Kingdom level 
#pruned_filtered_phy <- subset_taxa(pruned_filtered_phy, !is.na(Kingdom)) [old, there actually aren't any NA!]
pruned_filtered_phy <- subset_taxa(pruned_filtered_phy, Kingdom != "Unassigned")

#filter out archea and bacteria, plus mitochondria and plastids 
pruned_filtered_phy <- subset_taxa(pruned_filtered_phy, !(Kingdom %in% c("Archaea", "Bacteria", "Eukaryota:mito", "Eukaryota:plas")))

#exclude vertebrates (family Craniata)
pruned_filtered_phy <- subset_taxa(pruned_filtered_phy, !(Family %in% c("Craniata")))

#remove cephalapods
pruned_filtered_phy <- subset_taxa(pruned_filtered_phy, !(Genus %in% c('Cephalopoda')))

#remove L62, 83, and 94
pruned_filtered_phy<-subset_samples(pruned_filtered_phy,  !latrine_trt_month %in% c("L83_latrine_wet", 'L83_control_wet', "L94_latrine_wet", 'L94_control_wet', "L62_control_wet", 'L62_latrine_wet'))

#print the two to compare
pruned_filtered_phy
dup_phy


#run some checks to make sure filtering worked correctly 
tax_df <- as.data.frame(tax_table(pruned_filtered_phy))
#sum(is.na(tax_df$Kingdom)) #should be 0 if NA removed at Kingdom level [OLD]
sum(tax_df$Kingdom == "Unassigned", na.rm = TRUE)
table(tax_df$Kingdom) #should just be Eukaryota 
table(tax_df$Family) #no Craniata 

#looks good, save to final object 
final_filtered_phy<- pruned_filtered_phy




#save it as R file so it can be easily loaded. at this point I recommend continuing
#through the rarefying step and save the rarefied file instead
saveRDS(final_filtered_phy, file="final_filtered_phy_silva_18s") #use whatever file path for where you want to save it
final_filtered_phy<-readRDS("final_filtered_phy_silva_18s")


######## 3. Rarefying. if you want to skip to the rarefying and not do all of these steps (3a-3d) 
#that visualize the number of reads and determine the rarefy value, skip to section 3e

###3a. Plot histogram of number of reads per sample
reads_per_sample<- data.frame(sum=sample_sums(pruned_filtered_phy))
ggplot(reads_per_sample, aes(x=sum))+
  geom_histogram(binwidth=2500)


### 3b. Determine the minimum number of reads
smin<- min(sample_sums(final_filtered_phy))
smin 
#top 5 minimums 
sums <- sample_sums(final_filtered_phy)
top8_min <- sort(sums)[1:8]
top8_min

### 3c. Plot the rarefaction curves to determine sampling depth

#For all the data
otu.matrix = otu_table(final_filtered_phy) #make data into data frame
otu.matrix = as.data.frame(t(otu.matrix))
sample_names = rownames(otu.matrix) #add sample names



#plot
otu.rarecurve = vegan:rarecurve(otu.matrix, step = 50, label = F, xlim=c(0,15000))
abline(v=1939)
abline(v=5856)
abline(v=7068)
abline(v=10578)
abline(v=11022)
abline(v=11396)
text(x = c(1939, 5856, 7068, 10578, 11022, 11396),
     y = rep(par("usr")[3] + 0.9 * diff(par("usr")[3:4]), 4),
     labels = c("1939", "5856", "7068", "10578", "11022", "11396"),
     srt = 90,
     adj = 1,
     col = "blue")

#we have looked at the curves and now decided which value to use to rarefy

### 3d. Rarefy using the chosen value. rngseed sets the seed for us within the function [OLD method - see new below]
#filt_rare_phy <- rarefy_even_depth(pruned_filtered_phy, rngseed = 200, sample.size=12987)


### 3d-2: Rarefying with mirlyn
#lib size should be 12987
#for 100 reps

#rarefy data
mirl_object_100<- mirl(final_filtered_phy, libsize=11396, set.seed=200, trimOTUs=T, replace=F, rep=100)

#identify which samples were dropped
reads <- sample_sums(final_filtered_phy)
reads[reads < 11396] #these are the dropped samples 

#make an empty object to put the ASV tables in
mirl_otu_100 <- vector("list", length(mirl_object_100))

#extract otu tables from each rarefied phyloseq and add to the empty object above
for (i in 1:length(mirl_object_100)){
  colnames(mirl_object_100[[i]]@otu_table) <- paste0(colnames(mirl_object_100[[i]]@otu_table))
  (mirl_otu_100[[i]] <- mirl_object_100[[i]]@otu_table)
}



#make metadata file with the correct samples (remove ones dropped during filtering)
sample_id<- data.frame(final_filtered_phy@sam_data) 
sample_id$Samples<- row.names(sample_id)
sample_id<- sample_id %>% 
  filter(!Samples %in% c(13, 144, 148, 2, 9))


sample_id <- sample_id$Samples

#make empty list for each sample
average_counts_100 <- vector("list", length(sample_id))

#give how many reps you will do
rep_100<-1:100
#make empty list to hold 100 dataframes
iter_list_100<- vector('list', length(rep_100))

#rewrite loop to select columns from each rep, then average them and put them in new otu table
for (i in 1:length(sample_id) ){
  for (j in rep_100){
    iter_list_100[[j]]<-dplyr::select(as.data.frame(mirl_otu_100[[j]]),i) #this selects each individual iteration's otu table  
    iter_list_100[[j]]$ASVname<- row.names(iter_list_100[[j]])
  }
  
  sample_df_100<- reduce(iter_list_100[rep_100], full_join, by='ASVname')
  sample_df_100[is.na(sample_df_100)]<-0
  row.names(sample_df_100)<- sample_df_100$ASVname
  sample_df_100<- sample_df_100[,c(1, 3:(1+length(rep_100)))]
  sample_average_100 <- data.frame(rowMeans(sample_df_100))
  colnames(sample_average_100) <- sample_id[[i]]
  average_counts_100[[i]] <- sample_average_100
}
average_count_df_100 <- do.call(cbind, average_counts_100)

write.csv(x=average_count_df_100, file="100rep_averaged_OTUtable_18S.csv")


######## do some checks
#check that they all have the rarefied number of ASVs
colSums(average_count_df_100)

#is this close to the number for the whole data frame?
sum(iter_list_100[[1]]$`99`!=0)
sum(average_count_df_100$`99` !=0)
#relatively close
#not exact because there could be asvs not present here that are present in other iterations

#add to phyloseq
mirl_phyloseq <- final_filtered_phy
mirl_phyloseq@otu_table@.Data <- as.matrix(average_count_df_100)

rowSums(mirl_phyloseq@otu_table)==rowSums(average_count_df_100)  #should print a bunch of "TRUE"

#compare the two phyloseqs just to see and confirm that the expected number of samples are present
final_filtered_phy
mirl_phyloseq

#save to final phyloseq name used for analyses 
filt_rare_phy<- mirl_phyloseq


# Final phyloseq (filtered and rarefied) ----
#save the data as an R file so it doesn't have to be loaded each time.
#now when you start R, you can load the metadata and waypoints in step 1a. and skip
#steps 1b-3e

saveRDS(filt_rare_phy, file="filt_rare_phy_18s.rds") #use whatever file path for where you want to save it
filt_rare_phy_18s<-readRDS("filt_rare_phy_18s.rds")






# Alpha Diversity ----
############
############

###### 4. Diversity Analysis 

## NECESSARY Calculate Diversity ----
### 4a. Calculate diversity for GENUS level test
# recode all NAs as incertae sedis so they are counted as the same
taxa_all<- data.frame(tax_table(filt_rare_phy_18s))
taxa_all$Phylum[is.na(taxa_all$Phylum)] <- "Incertae_Sedis"
taxa_all$Class[is.na(taxa_all$Class)] <- "Incertae_Sedis"
taxa_all$Order[is.na(taxa_all$Order)] <- "Incertae_Sedis"
taxa_all$Family[is.na(taxa_all$Family)] <- "Incertae_Sedis"
taxa_all$Genus[is.na(taxa_all$Genus)] <- "Incertae_Sedis"

filt_rare_phy_18s@tax_table<-tax_table(as.matrix(taxa_all))

# for now we are calculating diversity for just things IDed to genus
#this means richness is how many individual genera are in each sample and Shannon
#diversity is calculated where each individual genus is a "thing"

#filter out incertae sedis for genus level alpha diversity
taxa_noNA_genus<- taxa_all %>% 
  filter(Genus != 'Incertae_Sedis')
phy_noNA_genus<- filt_rare_phy_18s
phy_noNA_genus@tax_table<- tax_table(as.matrix(taxa_noNA_genus))

#glomerate by genus
phy_noNA_genus_glom<-tax_glom(phy_noNA_genus, taxrank='Genus', NArm=T)

#view taxa table
glom_taxa<- data.frame(tax_table(phy_noNA_genus_glom))

#take out asv table to calculate diversity
glom_asv<- data.frame(otu_table(phy_noNA_genus_glom))
glom_asv<- data.frame(t(glom_asv), check.names = F)

### All data Shannon
#calculate shannon's diversity
shan_div_glom<-data.frame(diversity(glom_asv, index='shannon'))

#rename column to Shannon
colnames(shan_div_glom)[1]<- 'Shannon'

#make column with sample id to merge with metadata
shan_div_glom$`SampleID`<- row.names(shan_div_glom)

#for some reason it added X to the beginning of the sample names so I removed it here:
shan_div_glom$`SampleID`<-sub('.', '', shan_div_glom$`SampleID`)

#merge with the metadata so we can run a model and filter out the stuff already filtered out
metadata_filt<-metadata %>%
  left_join(shan_div_glom, by='SampleID') %>%
  filter(!is.na(Shannon))

#All Richness
rich_glom<- data.frame(specnumber(glom_asv))

#rename column to Richness
colnames(rich_glom)[1]<- 'Observed'

#make column with sample id to merge with metadata
rich_glom$`SampleID`<- row.names(rich_glom)

# #for some reason it added X to the beginning of the sample names so I removed it here:
rich_glom$`SampleID`<-sub('.', '', rich_glom$`SampleID`)

# #merge with metadata
metadata_filt<- metadata_filt %>% 
  left_join(rich_glom, by='SampleID') %>% 
  filter(!is.na(Observed))

# all inverse simpson
#calculate shannon's diversity
invsimp_glom<-data.frame(diversity(glom_asv, index='invsimpson'))

#rename column to Shannon
colnames(invsimp_glom)[1]<- 'InvSimpson'

#make column with sample id to merge with metadata
invsimp_glom$`SampleID`<- row.names(invsimp_glom)

#for some reason it added X to the beginning of the sample names so I removed it here:
invsimp_glom$`SampleID`<-sub('.', '', invsimp_glom$`SampleID`)

metadata_filt<- metadata_filt %>% 
  left_join(invsimp_glom, by='SampleID')

# all Pielou evenness
#######CHANGE NUMBER 
metadata_filt$Pielou<- metadata_filt$Shannon/ log(146)

#add elevation
#format latrine names
waypoints$latrineF<-gsub("^.{0,4}", "", waypoints$latrine)
waypoints$latrine<- paste("L", waypoints$latrineF, sep='') 

#sometimes it makes weird column names so try just highlighting this code and rerunning it
metadata_filt<- metadata_filt %>% 
  left_join(waypoints, by='latrine') 


## NECESSARY Calculate Diversity ----
### 4a. Calculate diversity for ASV level test
# 
# ### All data shannon
all_shan_div<-estimate_richness(filt_rare_phy_18s, measures='Shannon')
all_shan_div$`SampleID`<- row.names(all_shan_div)

#merge with the metadata so we can run a model and filter out the stuff already filtered out
metadata_filt<-metadata %>%
  left_join(all_shan_div, by='SampleID') %>%
  filter(!is.na(Shannon))

#doesn't like non-integers from rarefaction averaging, so calc richness with vegan specnumber() instead
otu_mat <- as(otu_table(filt_rare_phy_18s), "matrix") 
if(taxa_are_rows(filt_rare_phy_18s)) {otu_mat <- t(otu_mat)} #transpose

richness <- specnumber(otu_mat) #calculate richness
sample_ids <- sample_names(filt_rare_phy_18s) #grab samp IDs
all_richness <- data.frame(SampleID = sample_ids, Observed = richness)

# #merge with metadata
metadata_filt<- metadata_filt %>% 
  left_join(all_richness, by='SampleID') %>% 
  filter(!is.na(Observed))

# all inverse simpson
all_simpson<- estimate_richness(filt_rare_phy_18S, measures='InvSimpson')
all_simpson$`SampleID`<- row.names(all_simpson)

metadata_filt<- metadata_filt %>% 
  left_join(all_simpson, by='SampleID')

# all Pielou evenness
overall_S <- sum(colSums(otu_mat) > 0) #calc total number of ASVs
overall_S
metadata_filt$Pielou<- metadata_filt$Shannon/ log(overall_S)


#add elevation
#format latrine names
waypoints$latrineF<-gsub("^.{0,4}", "", waypoints$latrine)
waypoints$latrine<- paste("L", waypoints$latrineF, sep='') 

#sometimes it makes weird column names so try just highlighting this code and rerunning it
metadata_filt<- metadata_filt %>% 
  left_join(waypoints, by='latrine') 

#also make an elevation reference column that is elevation from a certain point
min(metadata_filt$elevation)
metadata_filt$elevation_ref<- (metadata_filt$elevation-5111)
#instead of true elevation, this is elevation from ~5100m

#format slope and aspect data to be merged
slope_aspect$latrine<- slope_aspect$Latrine
#join them together
metadata_filt<- metadata_filt %>% 
  left_join(slope_aspect, by='latrine')


#add Vicuna RAI 
vicugnaRAI<- dplyr::select(vicugnaRAI, latrine, RAI_vicugna )

metadata_filt<- metadata_filt %>% 
  left_join(vicugnaRAI, by='latrine')

#add chronosequences
metadata_filt<- metadata_filt %>% 
  left_join(chrono, by='latrine')



############ make the things going into the models factors
## Wet Subset Models ----
### 4c. Run models wet data subset by soil age 

# Wet season metadata
metadata_wet<- metadata_filt %>% 
  filter(`month-collected`=='wet')

#factor the variables in the models and scale elevation
metadata_wet<-metadata_wet %>% 
  mutate(treatment=as.factor(treatment), soilAge=as.factor(soilAge), class=factor(class, levels=c('LIA', 'LIA-1931','1931-1962','1984-2024'))) %>% 
  mutate(elevation_sc= scale(elevation))
names(metadata_wet)
#when you view the data, the elevation scaled column has a weird name but R says its name is elevation_sc
# also it is centering it because that is the default
Wel_mean<-mean(metadata_wet$elevation)
Wel_sd<-sd(metadata_wet$elevation)
#to calculate an elevation's value on the scaled scale, subtract the mean and divide by the sd
# (5100-Wel_mean)/Wel_sd

#to filter into LIA and RGM 
metadata_wet_LIA <- metadata_wet %>%
  filter(soilAge == "lia")
metadata_wet_RGM <- metadata_wet %>%
  filter(soilAge == "rgm")


### Richness ----
#wet season richness using the reference elevation
m_wet_rich<- lmer(Observed~treatment*soilAge+elevation_sc*treatment+(1|latrine_trt_month)+(1|latrine), data=metadata_wet, na.action='na.fail')
summary(m_wet_rich)
Anova(m_wet_rich, type='III')
emmeans(m_wet_rich, pairwise~treatment*soilAge)
qqnorm(residuals(m_wet_rich))

#use exp to backtransform bc on log scale

#the coefficients work exactly how we'd expect. at low elevations, controls are move diverse but than
# at higher elevations, latrines are more diverse. I calculated this by doing the math first and then
# exp() of the final answer. I added treatment and soil age terms as necessary and then for elevation you multiply
# the coefficient times the scale(elevation) value (ex. -.34408*1.89) and if it's latrines you do that AND
# the latrine:elevation coefficient times the scale(elevation) value (ex -.34408*1.89 + .25140*1.89)

#do LIA to see if treatment is significant
wet_richLIA<- glmer.nb(Observed~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metadata_wet_LIA)
summary(wet_richLIA)
Anova(wet_richLIA)

#do RGM wet to see if treatment is significant
wet_richRGMw<- glmer.nb(Observed~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metadata_wet_RGM)
summary(wet_richRGMw)
Anova(wet_richRGMw)

#compare AIC to model without elevation
m_wet_rich<-lmer(Observed~treatment*soilAge+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_rich)
Anova(m_wet_rich)
ICtab(m_wet_richNB, m_wet_rich)

#compare to soil Age null model. this tests if having soil age at all in the model makes it better
m_wet_rich_nullS<- glmer.nb(Observed~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
lrtest(m_wet_richNB, m_wet_rich_nullS) #if p value is sig, then the regular model is better than the null model

#compare to interaction null model. this tests if the soil age interaction is significant
m_wet_rich_nullI<- glmer.nb(Observed~treatment*elevation_ref+soilAge+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
lrtest(m_wet_richNB, m_wet_rich_nullI)

#make vicuna rai model
m_wet_rich_rai<- lmer(Observed~treatment*RAI_vicugna+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_rich_rai)
Anova(m_wet_rich_rai, type='III')


#chronosequence model
m_wet_rich_chrono<- lmer(Observed~treatment*class+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_rich_chrono)
Anova(m_wet_rich_chrono, type='III')
emmeans(m_wet_rich_chrono, pairwise~treatment*class)
qqnorm(residuals(m_wet_rich_chrono))

#see plots for regressions / boxplots 



### Shannon's Diversity ----
#Wet season Shannon's diversity using reference elevation
m_wet_shan_div<-lmer(Shannon~treatment*soilAge+elevation_sc*treatment+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_shan_div)
Anova(m_wet_shan_div, type='III')
emmeans(m_wet_shan_div, pairwise~treatment*soilAge)

#check model assumptions
qqnorm(residuals(m_wet_shan_div)) #checking normality

#compare to model without elevation
m_wet_shan<-lmer(Shannon~treatment*soilAge+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_shan)
Anova(m_wet_shan)
ICtab(m_wet_shan, m_wet_shan_div)

#compare to soil Age null model
m_wet_shan_nullS<- lmer(Shannon~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
lrtest(m_wet_shan_div, m_wet_shan_nullS)

#compare to interaction null model
m_wet_shan_nullI<- lmer(Shannon~treatment*elevation_sc+soilAge+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
lrtest(m_wet_shan_div, m_wet_shan_nullI)


#make vicuna rai model
m_wet_shan_rai<- lmer(Shannon~treatment*RAI_vicugna+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_shan_rai)
Anova(m_wet_shan_rai, type='III')


#chronosequence model
m_wet_shan_chrono<- lmer(Shannon~treatment*class+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_shan_chrono)
Anova(m_wet_shan_chrono, type='III')
emmeans(m_wet_shan_chrono, pairwise~treatment*class)
qqnorm(residuals(m_wet_shan_chrono))



### Inv Simpson's Diversity ----
#wet season Simpson with reference elevation
m_wet_simp<-lmer(InvSimpson~treatment*soilAge+elevation_sc*treatment+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_simp)
Anova(m_wet_simp, type='III')
emmeans(m_wet_simp, pairwise~treatment*soilAge)

#check model assumptions
qqnorm(residuals(m_wet_simp)) #checking normality

#compare to soil Age null model
m_wet_simp_nullS<- lmer(InvSimpson~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
lrtest(m_wet_simp, m_wet_simp_nullS)

#compare to interaction null model
m_wet_simp_nullI<- lmer(InvSimpson~treatment*elevation_sc+soilAge+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
lrtest(m_wet_simp, m_wet_simp_nullI)

#make vicuna rai model
m_wet_simp_rai<- lmer(InvSimpson~treatment*RAI_vicugna+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_simp_rai)
Anova(m_wet_simp_rai, type='III')


#chronosequence model
m_wet_simp_chrono<- lmer(InvSimpson~treatment*class+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_simp_chrono)
Anova(m_wet_simp_chrono, type='III')
emmeans(m_wet_simp_chrono, pairwise~treatment*class)
qqnorm(residuals(m_wet_simp_chrono))



### Pielou evenness ----
#logit transform Pielou to use a linear model with it
m_wet_pie<- lmer(logit(Pielou)~treatment*soilAge+elevation_sc*treatment+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_pie)
Anova(m_wet_pie, type='III')
qqnorm(residuals(m_wet_pie))
emmeans(m_wet_pie, pairwise~treatment*soilAge)

#make vicuna rai model
m_wet_pie_rai<- lmer(Pielou~treatment*RAI_vicugna+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_pie_rai)
Anova(m_wet_pie_rai, type='III')


#chronosequence model
m_wet_pie_chrono<- lmer(Pielou~treatment*class+(1|latrine_trt_month)+(1|latrine), data=metadata_wet)
summary(m_wet_pie_chrono)
Anova(m_wet_pie_chrono, type='III')
emmeans(m_wet_pie_chrono, pairwise~treatment*class)
qqnorm(residuals(m_wet_pie_chrono))



## RGM Subset Models ----
### 4d. Run model for RGM data

#metadata file
metadata_RGM<-metadata_filt %>% 
  filter(soilAge=='rgm')

metadata_RGM<-metadata_RGM %>% 
  mutate(treatment=as.factor(treatment), `month-collected`=as.factor(`month-collected`)) %>% 
  mutate(elevation_sc=scale(elevation))
str(metadata_RGM)

### Richness----
#RGM model richness with scaled elevation
m_season_richNB<- glmer.nb(Observed~treatment*`month-collected`+elevation_sc*treatment+(1|latrine_trt_month)+(1|latrine), data=metadata_RGM)
summary(m_season_richNB)
Anova(m_season_richNB, type='III')
emmeans(m_season_richNB, pairwise~treatment*`month-collected`)

#compare to season null model. tests if having season at all is better than not having it
m_season_rich_nullS<- glmer.nb(Observed~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metadata_RGM)
lrtest(m_season_richNB, m_season_rich_nullS)

#compare to interaction null model. tests if interaction with season is better than non interaction
m_season_rich_nullI<- glmer.nb(Observed~treatment*elevation_sc+`month-collected`+(1|latrine_trt_month)+(1|latrine), data=metadata_RGM)
lrtest(m_season_richNB, m_season_rich_nullI) # if pvalue is sig, interaction is better than not

#plot the richness data to see interaction between treatment and elevation
ggplot(metadata_RGM, aes(x=elevation, y=Observed, color=trt_month))+
  geom_point()+
  geom_smooth(method='lm') #default 95% CI


### Shannon's Diversity----
# RGM model shannon with reference elevation
m_season_shan<- lmer(Shannon~treatment*`month-collected`+elevation_sc*treatment+(1|latrine_trt_month)+(1|latrine), data=metadata_RGM)
summary(m_season_shan)
Anova(m_season_shan)

#check model assumptions
qqnorm(residuals(m_season_shan)) #checking normality

#compare to season null model
m_season_shan_nullS<- lmer(Shannon~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metadata_RGM)
lrtest(m_season_shan, m_season_shan_nullS)

#compare to interaction null model
m_season_shan_nullI<- lmer(Shannon~treatment*elevation_sc+`month-collected`+(1|latrine_trt_month)+(1|latrine), data=metadata_RGM)
lrtest(m_season_shan, m_season_shan_nullI)

### Inv Simpson's Diversity----
#RGM simpson with reference elevation
m_RGM_simp<-lmer(InvSimpson~treatment*`month-collected`+elevation_sc*treatment+(1|latrine_trt_month)+(1|latrine), data=metadata_RGM)
summary(m_RGM_simp)
Anova(m_RGM_simp)


#check model assumptions
qqnorm(residuals(m_RGM_simp)) #checking normality


#compare to season null model
m_RGM_simp_nullS<- lmer(InvSimpson~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metadata_RGM)
lrtest(m_RGM_simp, m_RGM_simp_nullS)

#compare to interaction null model
m_RGM_simp_nullI<- lmer(InvSimpson~treatment*elevation_sc+`month-collected`+(1|latrine_trt_month)+(1|latrine), data=metadata_RGM)
lrtest(m_RGM_simp, m_RGM_simp_nullI)


### Pielou evenness ----
#logit transform Pielou to use a linear model with it
m_RGM_pie<- lmer(logit(Pielou)~treatment*`month-collected`+elevation_sc*treatment+(1|latrine_trt_month)+(1|latrine), data = metadata_RGM)
summary(m_RGM_pie)
Anova(m_RGM_pie, type='III')
qqnorm(residuals(m_RGM_pie))







### RGM Dry----
metaDryRGM_both<- metadata_filt %>% 
  filter(soilAge=='rgm' & `month-collected`=='dry')

metaDryRGM_both<-metaDryRGM_both %>% 
  mutate(treatment=as.factor(treatment), `month-collected`=as.factor(`month-collected`), class=as.factor(class)) %>% 
  mutate(elevation_sc=scale(elevation))

#richness
m_dry_rich<- lmer(Observed~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
summary(m_dry_rich)
Anova(m_dry_rich, type='III')
qqnorm(residuals(m_dry_rich))

#make vicuna rai model
m_dry_rich_rai<- lmer(Observed~treatment*RAI_vicugna+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
summary(m_dry_rich_rai)
Anova(m_dry_rich_rai, type='III')


ggplot(metaDryRGM_both, aes(x=RAI_vicugna, y=Observed))+
  geom_point(aes(color=treatment))

#chronosequence model
m_dry_rich_chrono<- lmer(Observed~treatment*class+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
summary(m_dry_rich_chrono)
Anova(m_dry_rich_chrono, type='III')
emmeans(m_dry_rich_chrono, pairwise~treatment*class)
qqnorm(residuals(m_dry_rich_chrono))



#Shannon
m_dry_shan<- lmer(Shannon~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
summary(m_dry_shan)
Anova(m_dry_shan)
qqnorm(residuals(m_dry_shan))

#make vicuna rai model
m_dry_shan_rai<- lmer(Shannon~treatment*RAI_vicugna+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
Anova(m_dry_shan_rai, type='III')
summary(m_dry_shan_rai)

#chronosequence model
m_dry_shan_chrono<- lmer(Shannon~treatment*class+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
summary(m_dry_shan_chrono)
Anova(m_dry_shan_chrono, type='III')
emmeans(m_dry_shan_chrono, pairwise~treatment*class)
qqnorm(residuals(m_dry_shan_chrono))

#Inv Simpson
m_dry_simp<- lmer(InvSimpson~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
summary(m_dry_simp)
Anova(m_dry_simp)
qqnorm(residuals(m_dry_simp))

#make vicuna rai model
m_dry_simp_rai<- lmer(InvSimpson~treatment*RAI_vicugna+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
Anova(m_dry_simp_rai, type='III')
summary(m_dry_simp_rai)

#chronosequence model
m_dry_simp_chrono<- lmer(InvSimpson~treatment*class+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
summary(m_dry_simp_chrono)
Anova(m_dry_simp_chrono, type='III')
emmeans(m_dry_simp_chrono, pairwise~treatment*class)
qqnorm(residuals(m_dry_simp_chrono))

# Pielou evenness 
#logit transform Pielou to use a linear model with it
m_dry_pie<- lmer(logit(Pielou)~treatment*elevation_sc+(1|latrine_trt_month)+(1|latrine), data = metaDryRGM_both)
summary(m_dry_pie)
Anova(m_dry_pie, type='III')
qqnorm(residuals(m_dry_pie))

#make vicuna rai model
m_dry_pie_rai<- lmer(Pielou~treatment*RAI_vicugna+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
Anova(m_dry_pie_rai, type='III')
summary(m_dry_pie_rai)

#chronosequence model
m_dry_pie_chrono<- lmer(Pielou~treatment*class+(1|latrine_trt_month)+(1|latrine), data=metaDryRGM_both)
summary(m_dry_pie_chrono)
Anova(m_dry_pie_chrono, type='III')
emmeans(m_dry_pie_chrono, pairwise~treatment*class)
qqnorm(residuals(m_dry_pie_chrono))




# Beta Diversity ASV LEVEL ----
##### 5. Beta Diversity analysis

############## This step is neccessary 
### NECESSARY Reroot the tree.----
# 5a. Reroot the tree
#It has to be binary but now it is not since we trimmed it
ps_tree<- phy_tree(filt_rare_phy_18s) #put tree into an object
is.binary(ps_tree) #asking if it is binary. if false, go to next step

phy_tree(filt_rare_phy_18s)<-multi2di(ps_tree) #fix the tree and put it back in the phyloseq
is.binary(phy_tree(filt_rare_phy_18s)) #check if it's binary, should be true



## NOT NECESSARY Compare subsamples (replicates) Permanova----

###permanova for replicates 1 and 2 to see if they differ

#factor the variables
metadata_factored<- metadata_filt
metadata_factored$replicate<- as.factor(metadata_factored$replicate)
metadata_factored$latrine_trt<- as.factor(metadata_factored$latrine_trt)

#reorder the metadata to match the order of the phyloseq
sampr<- sample_data(filt_rare_phy_18s) #pull out data from phyloseq

#order metadata to match that from phyloseq
metadata_factored_rep<-metadata_factored[ order(match(metadata_factored$`SampleID`, row.names(sampr))), ]

set.seed(200)
#run permanova
#permanova<- adonis2(distance(filt_rare_phy, method='wunifrac')~replicate, data=metadata_factored_rep, by='terms')
#permanova

#new way so that it is testing replicate but within each latrine. we think this is the proper way to test replicate variation
perm_rep<- adonis2(distance(filt_rare_phy_18s, method='wunifrac')~replicate*latrine_trt, data=metadata_factored_rep, by='terms')
perm_rep



## NECESSARY Run this step to get filtered phyloseq that is used for beta diversity----
# 5b. Change ASV names 
metadata_factored<- metadata_filt
metadata_factored$treatment<- as.factor(metadata_factored$treatment)
metadata_factored$soilAge<- as.factor(metadata_factored$soilAge)
metadata_factored$`month-collected`<- as.factor(metadata_factored$`month-collected`)
metadata_factored$trt_month<- as.factor(metadata_factored$trt_month)
metadata_factored$trt_soilAge<- as.factor(metadata_factored$trt_soilAge)

## filter out rep 2
filt_rare_rep2 <- subset_samples(filt_rare_phy_18s, replicate==1 |row.names(filt_rare_phy_18s@sam_data) %in% c('10', '14') )
#change to include samples dropped during rarefy
sample_names(filt_rare_rep2) #check

#### make a dataframe that has the original and new asv names for convenience
## the original names are a random long string of characters so this makes them easier
## to reference and the data frame saves the original and new name so we know what's what

#pull out taxa table
taxa<- as.data.frame(tax_table(filt_rare_rep2))
#pull out the tree
tree<- phy_tree(filt_rare_rep2)
#make sure tips are in same order as taxa
sum(tree$tip.label==row.names(taxa)) 
#use this number in the following steps where it has seq(1,#), should also match
#number of asvs in phyloseq

#put original names into df
#use both_names to look up the original qiime2 asv name
both_names<- data.frame(original=rownames(taxa))
#rename asvs in taxa table and add to df
rownames(taxa)<- paste('ASV', seq(1,21724,1), sep='_')
both_names$number<- rownames(taxa)
#take out asv table and rename that too
asvfull<- otu_table(filt_rare_rep2)
rownames(asvfull)<- paste('ASV', seq(1,21724,1), sep='_')
#convert them into matrix to put back into phyloseq
tax<- tax_table(as.matrix(taxa))
otu<- otu_table(as.matrix(asvfull), taxa_are_rows = T)
sample<- sample_data(filt_rare_rep2)
#rename the tree tips too
tree$tip.label<- paste('ASV', seq(1,21724,1), sep='_')

#put all this back into phyloseq so ASVs now have a normal number name
rep2_named_phy<- phyloseq(otu, tax, sample, tree)
#this is the phyloseq we will use


## NECESSARY Subset for Wet season ----
# 5c. Subset wet data
#make metadata. make sure that of the samples that were rarefied out, they don't belong to the replicate
# that is being chosen for the beta diversity stuff. so here, L70 control wet rep 1 was dropped when we 
# rarefied so if we choose replicate 1, there is no L70 control wet representation in our data so
#we have to make sure it is chosen

metadata_wetF <- metadata_factored %>% 
  filter((`month-collected` == 'wet' & replicate == 1))

#filter the phyloseq for only wet samples
filt_rare_wet2<- subset_samples(rep2_named_phy, `month.collected` %in% ('wet'))

#reorder the metadata to match the order of the phyloseq
samp<- sample_data(filt_rare_wet2) #pull out data from phyloseq

metadata_wetF<-metadata_wetF[ order(match(metadata_wetF$`SampleID`, row.names(samp))), ]

#check both match 
sample_names(filt_rare_wet2)
metadata_wetF$SampleID

## NECESSARY Subset for RGM data----
# 5d. Subset RGM data
#metadata factored
metadata_RGMF<- metadata_factored %>% 
  filter(soilAge=='rgm' & (replicate==1 | `SampleID` %in% c('10','14'))) 

#make phyloseq for RGM only
filt_rare_RGM2<- rep2_named_phy%>% 
  subset_samples(soilAge %in% ('rgm'))

#order samples
sampR<- sample_data(filt_rare_RGM2) #pull out data from phyloseq

metadata_RGMF<-metadata_RGMF[order(match(metadata_RGMF$`SampleID`, row.names(sampR))), ]

#check both match 
sample_names(filt_rare_RGM2)
metadata_RGMF$SampleID

## Wet Subset Permanova ----
### 5e. Permanova test with wet season data
set.seed(200) ###VERY IMPORTANT, always keep the same

#run permanova
permanova_wet<- adonis2(distance(filt_rare_wet2, method='wunifrac')~treatment*soilAge, data=metadata_wetF, by='terms')
permanova_wet

#pairwise permanova to see which groups are different from each other
permanova_pairwise(distance(filt_rare_wet2, method='wunifrac'), grp=metadata_wetF$trt_soilAge, padj='holm')

# see Plots_18S file for code to make plots

## RGM Subset Permanova----
# 5f. RGM Permanova
set.seed(200)
#permanova
permanova_rgm<- adonis2(distance(filt_rare_RGM2, method='wunifrac')~treatment*`month-collected`, data=metadata_RGMF, by='terms')
permanova_rgm

#pairwise permanova to see which groups are different
permanova_pairwise(distance(filt_rare_RGM2, method='wunifrac'), grp=metadata_RGMF$trt_month, padj='holm')

### see Plots_18S file for code on how to make the plots

## Dry RGM Permanova----
filt_rare_RGM_dry<- subset_samples(filt_rare_RGM2, `month.collected` %in% ('dry'))
metaDryRGM<- metadata_RGMF %>% 
  filter(`month-collected`=='dry') 

#order samples
sampRd<- sample_data(filt_rare_RGM_dry) #pull out data from phyloseq

metaDryRGM<-metaDryRGM[order(match(metaDryRGM$`SampleID`, row.names(sampRd))), ]

set.seed(200)
#permanova
permanova_rgmD<- adonis2(distance(filt_rare_RGM_dry, method='wunifrac')~treatment, data=metaDryRGM, by='terms')
permanova_rgmD

# Homogeneity of dispersions----

##Wet subset: 
wet_betadis <- betadisper(distance(filt_rare_wet2, method = 'wunifrac'), group = metadata_wetF$trt_soilAge, type = 'median') #create betadisper object with dispersion distances  
wet_permutest <- permutest(wet_betadis, permutations = 999) #test for differences in dispersions 
wet_permutest
boxplot(wet_betadis)

#adonis2(dist(wet_betadis$distances) ~ metadata_wetF$trt_soilAge)

#RGM Dry subset: 
dryRGM_betadis<-betadisper(distance(filt_rare_RGM_dry, method='wunifrac'), group=metaDryRGM$treatment, type='median')
permutest(dryRGM_betadis)
boxplot(dryRGM_betadis)

#adonis2(dist(dryRGM_betadis$distances)~metaDryRGM$treatment)



# Simper ----
##### 7. Try Simper for testing community difference

## LIA Simper ----
# filter phyloseq and metadata for lia
filt_lia2<- subset_samples(filt_rare_wet2, soilAge %in% ('lia'))
metalia<- metadata_wetF %>% 
  filter(soilAge=='lia')

#extract asv table and transpose
asvLIA<- as.data.frame(otu_table(filt_lia2))
tasvLIA <- data.frame(t(asvLIA), check.names = F)

# run simper
simper_lia<- simper(tasvLIA, metalia$treatment)
simper_lia

#see the top 10
s_lia<- summary(simper_lia)
top10_LIA<-head(s_lia$control_latrine, n = 10)
#if this is null after you run it, change the latrine_control to control_latrine. or
# view simper_lia and see what the name is at the top of the output

simpLIA_asv<- row.names(top10_LIA)

#get actual taxa info
taxa_lia <- as.data.frame(tax_table(filt_lia2)) #taxonomy
simperLIA_taxa<-taxa_lia[row.names(taxa_lia) %in% simpLIA_asv,]


## RGM Wet Simper----
#filter phyloseq and metadata for rgm
filt_wet_rgm2<- subset_samples(filt_rare_wet2, soilAge %in% ('rgm'))
metargmW<- metadata_wetF %>% 
  filter(soilAge=='rgm')

#extract asv table and transpose
asvWrgm<- as.data.frame(otu_table(filt_wet_rgm2))
tasvWrgm <- data.frame(t(asvWrgm), check.names = F)

# run simper
simper_Wrgm<- simper(tasvWrgm, metargmW$treatment)
simper_Wrgm

#see the top 10
s_Wrgm<- summary(simper_Wrgm)
top10_Wrgm<-head(s_Wrgm$latrine_control, n = 10)
#if this is null after you run it, change the latrine_control to control_latrine. or
# view simper_Wrgm and see what the name is at the top of the output

simpWrgm_asv<- row.names(top10_Wrgm)

#get actual taxa info
taxa_Wrgm <- as.data.frame(tax_table(filt_wet_rgm2)) #taxonomy
simperWrgm_taxa<-taxa_Wrgm[row.names(taxa_Wrgm) %in% simpWrgm_asv,]


# RGM dry Simper----
## filter phyloseq and metadata for rgm dry
filt_rare_RGM_dry<- subset_samples(filt_rare_RGM2, `month.collected` %in% ('dry'))
metaDryRGM<- metadata_RGMF %>% 
  filter(`month-collected`=='dry') 

#extract asvs and transpose
asvDrgm<- as.data.frame(otu_table(filt_rare_RGM_dry))
tasvDrgm <- data.frame(t(asvDrgm), check.names = F)

# run simper
simper_Drgm<- simper(tasvDrgm, metaDryRGM$treatment)

#see the top 10
s_Drgm<- summary(simper_Drgm)
top10_Drgm<-head(s_Drgm$control_latrine, n = 10)
#if this is null after you run it, change the latrine_control to control_latrine. or
# view simper_Wrgm and see what the name is at the top of the output

simpDrgm_asv<- row.names(top10_Drgm)

#get actual taxa info
taxa_dry <- as.data.frame(tax_table(filt_rare_RGM_dry)) #taxonomy
simperDrgm_taxa<-taxa_dry[row.names(taxa_dry) %in% simpDrgm_asv,]



### Simper for all wet----
asvs_wet <- as.data.frame(otu_table(filt_rare_wet2)) #ASVs
tASV_wet <- data.frame(t(asvs_wet), check.names = F)

simper_trt<- simper(tASV_wet, metadata_wetF$trt_soilAge)
simper_trt

#see top contributing taxa
summary(simper_trt)$control_rgm_control_lia %>%
  round(3) %>%
  head()
s<- summary(simper_trt)
top20<-head(s$control_rgm_control_lia, n = 20)
simpWASV<- row.names(top20)

#get actual taxa info
taxaW<- data.frame(tax_table(filt_rare_wet2))
simperW_taxa<-taxaW[row.names(taxaW) %in% simpWASV,]

#top20 for each group 
top20_list <- lapply(s, head, 20)
top20_combined <- do.call(rbind, lapply(names(s), function(name) {
  df <- head(s[[name]], 20)
  df$source <- name
  df
}))

simpWASVall<- row.names(top20_combined)
#get actual taxa info
taxaWall<- data.frame(tax_table(filt_rare_wet2))
simperWall_taxa<-taxaWall[row.names(taxaWall) %in% simpWASVall,]

#see only significant species
comparisons <- c("control_rgm_latrine_rgm" , "control_rgm_control_lia" , "control_rgm_latrine_lia" , "latrine_rgm_control_lia", "latrine_rgm_latrine_lia", "control_lia_latrine_lia")

simper.results <- c()

for(i in 1:length(comparisons)) {
  require(tidyverse)
  temp <- summary(simper_trt)[as.character(comparisons[i])] %>%
    as.data.frame()
  colnames(temp) <- gsub(
    paste(comparisons[i],".", sep = ""), "", colnames(temp))
  temp <- temp %>%
    mutate(Comparison = comparisons[i],
           Position = row_number()) %>%
    rownames_to_column(var = "Species")
  simper.results <- rbind(simper.results, temp)
}

simper.results %>%
  filter(p <= 0.05) %>%
  dplyr::select(Species, average, Comparison, Position)

#see sum for the groups
simper.results %>%
  group_by(Comparison) %>%
  summarize(sum.average = sum(average))

### See Plots_18S file for code that plots these ASVs as arrows









# Indicator Taxa----
## Wet Subset (LIA & wet RGM) ----
### 8. Indicator analysis for Wet data
#using rarefied data and just replicate 1

## 8a. separate lia and rgm samples

filt_lia2<- subset_samples(filt_rare_wet2, soilAge %in% ('lia'))
filt_wet_rgm2<- subset_samples(filt_rare_wet2, soilAge %in% ('rgm'))

## 8b. extract data from the phyloseq and format
#extract taxa table
taxa_lia <- as.data.frame(tax_table(filt_lia2)) #taxonomy
taxa_Wrgm <- as.data.frame(tax_table(filt_wet_rgm2)) #taxonomy

#extract the asvs 
asvLIA<- as.data.frame(otu_table(filt_lia2))
asvWrgm<- as.data.frame(otu_table(filt_wet_rgm2))

#transpose the asv matrix 
dim(asvLIA)
tasvLIA <- data.frame(t(asvLIA), check.names = F)
rownames(tasvLIA)
colnames(tasvLIA)

tasvWrgm <- data.frame(t(asvWrgm), check.names = F)

#make vector with treatment to use for the test
treatment_lia<- sample_data(filt_lia2)
treatment_lia<- treatment_lia$treatment
treatment_Wrgm<- sample_data(filt_wet_rgm2)
treatment_Wrgm<- treatment_Wrgm$treatment

### 8c. Run the test
set.seed(200) ### Very Important

ind_lia<- multipatt(tasvLIA, treatment_lia, func='IndVal.g')
summary(ind_lia, indvalcomp=T)

ind_Wrgm<- multipatt(tasvWrgm, treatment_Wrgm, func='IndVal.g')
summary(ind_Wrgm, indvalcomp=T)

#put output into data frame with sign and p value
output_lia<- data.frame(ind_lia$sign)
output_Wrgm<- data.frame(ind_Wrgm$sign)


### 8d. Extract the significant ASVs for each treatment
sigL_lia<-output_lia %>% 
  filter(p.value<=.05) %>% 
  filter(s.latrine==1)

sigC_lia<-output_lia %>% 
  filter(p.value<=.05) %>% 
  filter(s.control==1)

sigL_Wrgm<-output_Wrgm %>% 
  filter(p.value<=.05) %>% 
  filter(s.latrine==1)

sigC_Wrgm<-output_Wrgm %>% 
  filter(p.value<=.05) %>% 
  filter(s.control==1)

#make a dataframe with the taxonomic info of the significant asvs
ind_taxaL_lia <- taxa_lia[rownames(taxa_lia) %in% rownames(sigL_lia), ]  
ind_taxaC_lia <- taxa_lia[rownames(taxa_lia) %in% rownames(sigC_lia), ]  

ind_taxaL_Wrgm <- taxa_Wrgm[rownames(taxa_Wrgm) %in% rownames(sigL_Wrgm), ]  
ind_taxaC_Wrgm <- taxa_Wrgm[rownames(taxa_Wrgm) %in% rownames(sigC_Wrgm), ]  


#join taxonomic info to the output data 
sigL_lia$ASV<- row.names(sigL_lia)
ind_taxaL_lia$ASV<- row.names(ind_taxaL_lia)
ind_taxaL_lia<- ind_taxaL_lia %>% left_join(sigL_lia, by='ASV')

sigC_lia$ASV<- row.names(sigC_lia)
ind_taxaC_lia$ASV<- row.names(ind_taxaC_lia)
ind_taxaC_lia<- ind_taxaC_lia %>% left_join(sigC_lia, by='ASV')

sigL_Wrgm$ASV<- row.names(sigL_Wrgm)
ind_taxaL_Wrgm$ASV<- row.names(ind_taxaL_Wrgm)
ind_taxaL_Wrgm<- ind_taxaL_Wrgm %>% left_join(sigL_Wrgm, by='ASV')

sigC_Wrgm$ASV<- row.names(sigC_Wrgm)
ind_taxaC_Wrgm$ASV<- row.names(ind_taxaC_Wrgm)
ind_taxaC_Wrgm<- ind_taxaC_Wrgm %>% left_join(sigC_Wrgm, by='ASV')

# you can collapse it to just get one row for each unique ID
unique_ind_taxaL_lia<- unique(ind_taxaL_lia)
unique_ind_taxaC_lia<- unique(ind_taxaC_lia)
unique_ind_taxaL_rgm<- unique(ind_taxaL_Wrgm)
unique_ind_taxaC_rgm<- unique(ind_taxaC_Wrgm)


#find taxa that are shared and different
#this is at the family level but you can change it to be at any level by changing
#the column name
#you can also get the asvs instead of the taxonomic info by changing column name
#to the ASV column
#shared between latrine lia and rgm
unique(ind_taxaL_lia[ind_taxaL_lia$ASV %in% ind_taxaL_Wrgm$ASV,8])
#unique to lia
unique(ind_taxaL_lia[!ind_taxaL_lia$Family %in% ind_taxaL_Wrgm$Family,5])
#unique to rgm
unique(ind_taxaL_Wrgm[!ind_taxaL_Wrgm$Family %in% ind_taxaL_lia$Family,5])

#find taxa that are shared between wet rgm and lia controls
unique(ind_taxaC_lia[ind_taxaC_lia$Family %in% ind_taxaC_Wrgm$Family,5])
#unique to lia
unique(ind_taxaC_lia[!ind_taxaC_lia$Family %in% ind_taxaC_Wrgm$Family,5])
#unique to wet rgm
unique(ind_taxaC_Wrgm[!ind_taxaC_Wrgm$Family %in% ind_taxaC_lia$Family,5])



## Dry RGM----
### 8B. Indicator analysis for RGM data
#using rarefied data and just replicate 1

## 8a. separate dry samples
filt_rare_RGM_dry<- subset_samples(filt_rare_RGM2, `month.collected` %in% ('dry'))

## 8b. extract data from the phyloseq and format
#extract taxa table
taxa_dry <- as.data.frame(tax_table(filt_rare_RGM_dry)) #taxonomy

#extract the asvs 
asvDrgm<- as.data.frame(otu_table(filt_rare_RGM_dry))

#transpose the asv matrix 
dim(asvDrgm)
tasvDrgm <- data.frame(t(asvDrgm), check.names = F)
rownames(tasvDrgm)
colnames(tasvDrgm)

#make vector with treatment
treatment_dry<- sample_data(filt_rare_RGM_dry)
treatment_dry<- treatment_dry$treatment

### 8c. Run the test
set.seed(200) ### Very Important

ind_dry<- multipatt(tasvDrgm, treatment_dry, func='IndVal.g')
summary(ind_dry, indvalcomp=T)

output_dry<- data.frame(ind_dry$sign)


### 8d. Extract the significant ASVs for each treatment
#make data with just latrine significant species
sigL_dry<-output_dry %>% 
  filter(p.value<=.05) %>% 
  filter(s.latrine==1)

#control significant species
sigC_dry<-output_dry %>% 
  filter(p.value<=.05) %>% 
  filter(s.control==1)


#add taxanomic info
ind_taxaL_dry <- taxa_dry[rownames(taxa_dry) %in% rownames(sigL_dry), ]  
ind_taxaC_dry <- taxa_dry[rownames(taxa_dry) %in% rownames(sigC_dry), ]  


#join the taxonomy and outputs
sigL_dry$ASV<- row.names(sigL_dry)
ind_taxaL_dry$ASV<- row.names(ind_taxaL_dry)
ind_taxaL_dry<- ind_taxaL_dry %>% left_join(sigL_dry, by='ASV')

sigC_dry$ASV<- row.names(sigC_dry)
ind_taxaC_dry$ASV<- row.names(ind_taxaC_dry)
ind_taxaC_dry<- ind_taxaC_dry %>% left_join(sigC_dry, by='ASV')

#collapse it so I just get one row for each unique ID
unique_ind_taxaL_dry<- unique(ind_taxaL_dry)
unique_ind_taxaC_dry<- unique(ind_taxaC_dry)

#find taxa that are shared between wet and dry rgm latrines
unique(ind_taxaL_dry[ind_taxaL_dry$Family %in% ind_taxaL_Wrgm$Family,5])
#unique to dry
unique(ind_taxaL_dry[!ind_taxaL_dry$Family %in% ind_taxaL_Wrgm$Family,5])
#unique to wet (not in dry, doesn't take into account lia taxa)
unique(ind_taxaL_Wrgm[!ind_taxaL_Wrgm$Family %in% ind_taxaL_dry$Family,5])

#find taxa that are shared between wet and dry rgm controls 
unique(ind_taxaC_dry[ind_taxaC_dry$Family %in% ind_taxaC_Wrgm$Family,5])
#unique to dry
unique(ind_taxaC_dry[!ind_taxaC_dry$Family %in% ind_taxaC_Wrgm$Family,5])
#unique to wet
unique(ind_taxaC_Wrgm[!ind_taxaC_Wrgm$Family %in% ind_taxaC_dry$Family,5])





# Differential Abundance----
## Wet RGM Latrine vs Control DA----

RGM2_phy_ASV<- filt_rare_rep2 %>% 
  subset_samples(soilAge %in% ('rgm'))


#factor treatment
rgm2_sampdata<- sample_data(RGM2_phy_ASV)
rgm2_sampdata$treatment<- as.factor(rgm2_sampdata$treatment)
RGM2_phy_ASV@sam_data<- rgm2_sampdata
str(RGM2_phy_ASV@sam_data)

rgmW_rep2_phy<- subset_samples(RGM2_phy_ASV, month.collected=='wet')

# test with just one rep and no RE
rgmWetTreatmentDA<-ancombc2(data = rgmW_rep2_phy, tax_level = "Genus",
                            fix_formula = "treatment", rand_formula = NULL,
                            p_adj_method = "holm", pseudo_sens = TRUE,
                            prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                            group = "treatment", struc_zero = TRUE, neg_lb = TRUE,
                            alpha = 0.05, n_cl = 2, verbose = TRUE,
                            global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                            iter_control = list(tol = 1e-2, max_iter = 20, 
                                                verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),
                            lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

#put primary results in data frame
rgmWetT_prim<-rgmWetTreatmentDA$res

#save it as an rds file
saveRDS(rgmWetT_prim, file='rgmWetT_prim.rds')
rgmWetT_prim<-readRDS('rgmWetT_prim.rds') #this is the BH one

#filter for what's significant
rgmWetTSig<-rgmWetT_prim %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T)

#extract taxa from phyloseq
rgmWet_taxa<- data.frame(tax_table(rgmW_rep2_phy))



# Plot log fold change
rgmWetT_DAplot<- rgmWetTSig %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T) %>% 
  dplyr::arrange(desc(lfc_treatmentlatrine)) %>% 
  dplyr::mutate(direct = ifelse(lfc_treatmentlatrine> 0, "Positive LFC", "Negative LFC"))

#make taxon and direction factors
rgmWetT_DAplot$taxon<- factor(rgmWetT_DAplot$taxon, levels=rgmWetT_DAplot$taxon)
rgmWetT_DAplot$direct<- factor(rgmWetT_DAplot$direct, levels = c("Positive LFC", "Negative LFC"))


fig_rgmWetT = rgmWetT_DAplot %>%
  ggplot(aes(x = taxon, y = lfc_treatmentlatrine, fill=direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_treatmentlatrine - se_treatmentlatrine, ymax = lfc_treatmentlatrine + se_treatmentlatrine), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes") + 
  scale_fill_manual(values=c('purple3', 'cyan3'), name=NULL, labels=c('Positive LFC (more in latrine)','Negative LFC (more in control)'))+
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x= element_text(hjust=1, angle=45))
fig_rgmWetT


#do the test at the phylum level
rgmWetTDA_phylum<-ancombc2(data = rgmW_rep2_phy, tax_level = "Phylum",
                           fix_formula = "treatment", rand_formula = NULL,
                           p_adj_method = "holm", pseudo_sens = TRUE,
                           prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                           group = "treatment", struc_zero = TRUE, neg_lb = TRUE,
                           alpha = 0.05, n_cl = 2, verbose = TRUE,
                           global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                           iter_control = list(tol = 1e-2, max_iter = 20, 
                                               verbose = TRUE),
                           em_control = list(tol = 1e-5, max_iter = 100),
                           lme_control = lme4::lmerControl(),
                           mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

rgmWet_phylum_prim<- rgmWetTDA_phylum$res

rgmWet_phylumSig<- rgmWet_phylum_prim %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T)

# Plot log fold change
rgmWetphylum_DAplot<- rgmWet_phylumSig %>% 
  dplyr::arrange(desc(lfc_treatmentlatrine)) %>% 
  dplyr::mutate(direct = ifelse(lfc_treatmentlatrine> 0, "Positive LFC", "Negative LFC"))

#make taxon and direction factors
rgmWetphylum_DAplot$taxon<- factor(rgmWetphylum_DAplot$taxon, levels=rgmWetphylum_DAplot$taxon)
rgmWetphylum_DAplot$direct<- factor(rgmWetphylum_DAplot$direct, levels = c("Positive LFC", "Negative LFC"))

fig_rgmWetphylum = rgmWetphylum_DAplot %>%
  ggplot(aes(x = taxon, y = lfc_treatmentlatrine, fill=direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_treatmentlatrine - se_treatmentlatrine, ymax = lfc_treatmentlatrine + se_treatmentlatrine), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(hjust=1, angle=90))
fig_rgmWetphylum

## Dry RGM Latrine vs Control----

#got the phyloseq from the previous step

#make treatment a factor
rgm2_sampdata$treatment<- as.factor(rgm2_sampdata$treatment)
RGM2_phy_ASV@sam_data<- rgm2_sampdata
str(RGM2_phy_ASV@sam_data)

#make it just for the dry samples
rgmD_rep2_phy<- subset_samples(RGM2_phy_ASV, month.collected=='dry')

# test with just rep 2 and no RE
rgmDryTreatmentDA<-ancombc2(data = rgmD_rep2_phy, tax_level = "Genus",
                            fix_formula = "treatment", rand_formula = NULL,
                            p_adj_method = "holm", pseudo_sens = TRUE,
                            prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                            group = "treatment", struc_zero = TRUE, neg_lb = TRUE,
                            alpha = 0.05, n_cl = 2, verbose = TRUE,
                            global = TRUE, pairwise = F, dunnet = F, trend = F,
                            iter_control = list(tol = 1e-2, max_iter = 20, 
                                                verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),
                            lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

#put primary results in data frame
rgmDryT_prim<-rgmDryTreatmentDA$res

#save it as an rds file
saveRDS(rgmDryT_prim, file='rgmDryT_prim.rds')
rgmDryT_prim<-readRDS('rgmDryT_prim.rds')

#filter for what's significant
rgmDryTSig<-rgmDryT_prim %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T)

#extract taxa from phyloseq
rgmDry_taxa<- data.frame(tax_table(rgmD_rep2_phy))


# Plot log fold change
rgmDryT_DAplot<- rgmDryTSig %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T) %>% 
  dplyr::arrange(desc(lfc_treatmentlatrine)) %>% 
  dplyr::mutate(direct = ifelse(lfc_treatmentlatrine> 0, "Positive LFC", "Negative LFC"))

#make taxon and direction factors
rgmDryT_DAplot$taxon<- factor(rgmDryT_DAplot$taxon, levels=rgmDryT_DAplot$taxon)
rgmDryT_DAplot$direct<- factor(rgmDryT_DAplot$direct, levels = c("Positive LFC", "Negative LFC"))


fig_rgmDryT = rgmDryT_DAplot %>%
  ggplot(aes(x = taxon, y = lfc_treatmentlatrine, fill=direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_treatmentlatrine - se_treatmentlatrine, ymax = lfc_treatmentlatrine + se_treatmentlatrine), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Dry Latrine vs Control") + 
  scale_fill_manual(values=c( 'cyan3'), name=NULL, labels=c('Negative LFC (more in control)'))+
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x= element_text(hjust=1, angle=45))
fig_rgmDryT

#do the test at the phylum level
rgmDryTDA_phylum<-ancombc2(data = rgmD_rep2_phy, tax_level = "Phylum",
                           fix_formula = "treatment", rand_formula = NULL,
                           p_adj_method = "holm", pseudo_sens = TRUE,
                           prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                           group = "treatment", struc_zero = TRUE, neg_lb = TRUE,
                           alpha = 0.05, n_cl = 2, verbose = TRUE,
                           global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                           iter_control = list(tol = 1e-2, max_iter = 20, 
                                               verbose = TRUE),
                           em_control = list(tol = 1e-5, max_iter = 100),
                           lme_control = lme4::lmerControl(),
                           mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

rgmDry_phylum_prim<- rgmDryTDA_phylum$res


rgmDry_phylumSig<- rgmDry_phylum_prim %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T)



## LIA vs RGM Control but using only 4 rgm locations that we chose based on location and availability----
##51,60,56,58
#get the phyloseq with regular asv names
wet2_phy_ASV<- filt_rare_rep2 %>% 
  subset_samples(month.collected=='wet')

## make our test variables factors
wet2_sampdata<- sample_data(wet2_phy_ASV)
wet2_sampdata$soilAge<- as.factor(wet2_sampdata$soilAge)
wet2_phy_ASV@sam_data<- wet2_sampdata
str(wet2_phy_ASV@sam_data)

#select only the 4 rgm samples we want (as well as the LIA ones)
soilAgeCDAphy<- wet2_phy_ASV %>% 
  subset_samples(latrine %in% c('L51','L56','L58','L60','L100','L101','L102','L104')) %>% 
  subset_samples(treatment=='control')

#run the test
soilAgeCDA<-ancombc2(data = soilAgeCDAphy, tax_level = "Genus",
                     fix_formula = "soilAge", rand_formula =NULL,
                     p_adj_method = "holm", pseudo_sens = TRUE,
                     prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                     group = "soilAge", struc_zero = T, neg_lb = T,
                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                     global = TRUE, pairwise = F, dunnet = F, trend = F,
                     iter_control = list(tol = 1e-2, max_iter = 20, 
                                         verbose = TRUE),
                     em_control = list(tol = 1e-5, max_iter = 100),
                     lme_control = lme4::lmerControl(),
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

soilAgeC_prim<- soilAgeCDA$res

saveRDS(soilAgeC_prim, file='soilAgeC_prim.rds')
soilAgeC_prim<-readRDS('soilAgeC_prim.rds')

soilAgeC_sig<- soilAgeC_prim %>% 
  filter(q_soilAgergm<.05 & passed_ss_soilAgergm==T)

soilAgeC_zero<- soilAgeCDA$zero_ind
soilAgeC_zeroLIA<- soilAgeC_zero %>% 
  filter(`structural_zero (soilAge = lia)`==T & `structural_zero (soilAge = rgm)`==F)
soilAgeC_zeroRGM<- soilAgeC_zero %>% 
  filter(`structural_zero (soilAge = lia)`==F & `structural_zero (soilAge = rgm)`==T)

#### do a test at the phylum level
#get LIA samples
wet2_sampdata$treatment<- as.factor(wet2_sampdata$treatment)
wet2_phy_ASV@sam_data<- wet2_sampdata

LIA_phylum_DA<- wet2_phy_ASV %>% 
  subset_samples(soilAge== 'lia')

#run the test
soilAgeDAPhylum<-ancombc2(data = LIA_phylum_DA, tax_level = "Phylum",
                          fix_formula = "treatment", rand_formula =NULL,
                          p_adj_method = "holm", pseudo_sens = TRUE,
                          prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                          group = "treatment", struc_zero = T, neg_lb = T,
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          global = TRUE, pairwise = TRUE, dunnet = F, trend = F,
                          iter_control = list(tol = 1e-2, max_iter = 20, 
                                              verbose = TRUE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))
soilAgeDAPhylum_pair<- soilAgeDAPhylum$res
soilAgeDAPhylumSig<- soilAgeDAPhylum_pair %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T)

## Plot all 3 phylum tests in one figure----

#first run the 3 different phylum level tests under the Wet L vs C, Dry L vs C, and LIA L vs C

#######first process the wet season comparison
#select just the columns of interest
wet_phyl_fig <- rgmWet_phylum_prim %>%
  dplyr::select(taxon, contains('latrine')) 

#round LFC, choose only taxa that are significant and passed sensitivity test, make data tidy
wet_phyl_fig_lfc <- wet_phyl_fig %>%
  dplyr::filter(diff_treatmentlatrine == 1 & passed_ss_treatmentlatrine==T) %>%
  dplyr::mutate(lfc1 = ifelse(diff_treatmentlatrine == 1, 
                              round(lfc_treatmentlatrine, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon) %>% 
  dplyr::select(group, value, taxon)

# recode the group so instead of lfc1 it says what the comparison is
wet_phyl_fig_lfc$group <- dplyr::recode(wet_phyl_fig_lfc$group, 
                                        `lfc1` = "RGM Wet Season")

###### second process dry season comparison
dry_phyl_fig<- rgmDry_phylum_prim %>% 
  dplyr::select(taxon, contains('latrine')) 

#round LFC, choose only taxa that are significant and passed sensitivity test, make data tidy
dry_phyl_fig_lfc <- dry_phyl_fig %>%
  dplyr::filter(diff_treatmentlatrine == 1 & passed_ss_treatmentlatrine==T) %>%
  dplyr::mutate(lfc1 = ifelse(diff_treatmentlatrine == 1, 
                              round(lfc_treatmentlatrine, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon) %>% 
  dplyr::select(group, value, taxon)

# recode the group so instead of lfc1 it says what the comparison is
dry_phyl_fig_lfc$group <- dplyr::recode(dry_phyl_fig_lfc$group, 
                                        `lfc1` = "RGM Dry Season")

##### third process the LIA comparison
lia_phyl_fig<- soilAgeDAPhylum_pair %>% 
  dplyr::select(taxon, contains('latrine'))

#round LFC, choose only taxa that are significant and passed sensitivity test, make data tidy
lia_phyl_fig_lfc <- lia_phyl_fig %>%
  dplyr::filter(diff_treatmentlatrine == 1 & passed_ss_treatmentlatrine==T) %>%
  dplyr::mutate(lfc1 = ifelse(diff_treatmentlatrine == 1, 
                              round(lfc_treatmentlatrine, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon) %>% 
  dplyr::select(group, value, taxon)

# recode the group so instead of lfc1 it says what the comparison is
lia_phyl_fig_lfc$group <- dplyr::recode(lia_phyl_fig_lfc$group, 
                                        `lfc1` = "LIA Wet Season")

#join all the comparisons together
phyl_fig<-full_join(lia_phyl_fig_lfc, wet_phyl_fig_lfc)
phyl_fig<- full_join(dry_phyl_fig_lfc, phyl_fig)

# make the figure
lo = floor(min(phyl_fig$value))
up = ceiling(max(phyl_fig$value))

fig_phyl = phyl_fig %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "cyan3", high = "purple3", mid = "white", 
                       na.value = "white", midpoint = 0, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = 'Phylum LFC Latrine-Control') +
  theme_classic() +
  theme(axis.text.x=element_text())
fig_phyl

#### make that plot comparing simper and ancombc2 taxa----
## try it first with lia vs rgm control

LIA_RGM_taxa<- data.frame(tax_table(LIA_RGM_DAphy))

soilAgeCC_taxa1<-LIA_RGM_taxa %>% 
  filter(Genus=='Abditibacterium') %>% 
  row.names()

LIA_RGM_taxa %>% 
  filter(Genus=='Blastocatella') %>% 
  row.names()

soilAgeCC_taxa1names<-both_names %>% 
  filter(original %in% soilAgeCC_taxa1) 

soilAgeCC_simp<-s[["control_rgm_control_lia"]]

soilAgeCC_simp_taxa1<-soilAgeCC_simp %>% 
  filter(row.names(soilAgeCC_simp) %in% soilAgeCC_taxa1names$number)
#I looked at it to see if simper agrees that this taxa is more prevalent in control rgm
#and it does for the most part (except for 4 asvs, and about half of the asvs are 0 for all the simper stats)

sum(soilAgeCC_simp_taxa1$average)
#express this as a proportion out of the total dissimilarity (sum of average column in total simper results?)

## i need to do this for multiple taxa but im not sure how to simplify the process
# and not have to do it 1 by 1 for each taxa. and we could show multiple comparisons on 
# the plot so there is more going on but that makes the process even longer


# Beta Diversity GENUS LEVEL ----
##### 5. Beta Diversity analysis

############## This step is neccessary 
### NECESSARY Reroot the tree.----
# 5a. Reroot the tree
#It has to be binary but now it is not since we trimmed it
ps_tree_genus<- phy_tree(phy_noNA_genus_glom) #put tree into an object
is.binary(ps_tree_genus) #asking if it is binary. if false, go to next step


## NOT NECESSARY Compare subsamples (replicates) Permanova----

###permanova for replicates 1 and 2 to see if they differ

#factor the variables
metadata_factored<- metadata_filt
metadata_factored$replicate<- as.factor(metadata_factored$replicate)
metadata_factored$latrine_trt<- as.factor(metadata_factored$latrine_trt)

#reorder the metadata to match the order of the phyloseq
sampr_genus<- sample_data(phy_noNA_genus_glom) #pull out data from phyloseq

#order metadata to match that from phyloseq
metadata_factored_rep<-metadata_factored[ order(match(metadata_factored$`SampleID`, row.names(sampr_genus))), ]

set.seed(200)
#run permanova
#permanova<- adonis2(distance(filt_rare_phy, method='wunifrac')~replicate, data=metadata_factored_rep, by='terms')
#permanova

#new way so that it is testing replicate but within each latrine. we think this is the proper way to test replicate variation
perm_rep_genus<- adonis2(distance(phy_noNA_genus_glom, method='wunifrac')~replicate*latrine_trt, data=metadata_factored_rep, by='terms')
perm_rep_genus



## NECESSARY Run this step to get filtered phyloseq that is used for beta diversity----
# 5b. Change ASV names 
metadata_factored<- metadata_filt
metadata_factored$treatment<- as.factor(metadata_factored$treatment)
metadata_factored$soilAge<- as.factor(metadata_factored$soilAge)
metadata_factored$`month-collected`<- as.factor(metadata_factored$`month-collected`)
metadata_factored$trt_month<- as.factor(metadata_factored$trt_month)
metadata_factored$trt_soilAge<- as.factor(metadata_factored$trt_soilAge)

## filter out rep 2
filt_rare_rep2_genus <- subset_samples(phy_noNA_genus_glom, replicate==1 |row.names(phy_noNA_genus_glom@sam_data) %in% c('10', '14') )
#change to include samples dropped during rarefy
sample_names(filt_rare_rep2_genus) #check

#### make a dataframe that has the original and new asv names for convenience
## the original names are a random long string of characters so this makes them easier
## to reference and the data frame saves the original and new name so we know what's what

#pull out taxa table
taxa<- as.data.frame(tax_table(filt_rare_rep2_genus))
#pull out the tree
tree<- phy_tree(filt_rare_rep2_genus)
#make sure tips are in same order as taxa
sum(tree$tip.label==row.names(taxa)) 
#use this number in the following steps where it has seq(1,#), should also match
#number of asvs in phyloseq

#put original names into df
#use both_names to look up the original qiime2 asv name
both_names<- data.frame(original=rownames(taxa))
#rename asvs in taxa table and add to df
rownames(taxa)<- paste('ASV', seq(1,146,1), sep='_')
both_names$number<- rownames(taxa)
#take out asv table and rename that too
asvfull<- otu_table(filt_rare_rep2_genus)
rownames(asvfull)<- paste('ASV', seq(1,146,1), sep='_')
#convert them into matrix to put back into phyloseq
tax<- tax_table(as.matrix(taxa))
otu<- otu_table(as.matrix(asvfull), taxa_are_rows = T)
sample<- sample_data(filt_rare_rep2_genus)
#rename the tree tips too
tree$tip.label<- paste('ASV', seq(1,146,1), sep='_')

#put all this back into phyloseq so ASVs now have a normal number name
rep2_named_phy_genus<- phyloseq(otu, tax, sample, tree)
#this is the phyloseq we will use


## NECESSARY Subset for Wet season ----
# 5c. Subset wet data
#make metadata. make sure that of the samples that were rarefied out, they don't belong to the replicate
# that is being chosen for the beta diversity stuff. so here, L70 control wet rep 1 was dropped when we 
# rarefied so if we choose replicate 1, there is no L70 control wet representation in our data so
#we have to make sure it is chosen

metadata_wetF <- metadata_factored %>% 
  filter((`month-collected` == 'wet' & replicate == 1))

#filter the phyloseq for only wet samples
filt_rare_wet2_genus<- subset_samples(rep2_named_phy_genus, `month.collected` %in% ('wet'))

#reorder the metadata to match the order of the phyloseq
samp_genus<- sample_data(filt_rare_wet2_genus) #pull out data from phyloseq

metadata_wetF<-metadata_wetF[ order(match(metadata_wetF$`SampleID`, row.names(samp_genus))), ]

#check both match 
sample_names(filt_rare_wet2_genus)
metadata_wetF$SampleID

## NECESSARY Subset for RGM data----
# 5d. Subset RGM data
#metadata factored
metadata_RGMF<- metadata_factored %>% 
  filter(soilAge=='rgm' & (replicate==1 | `SampleID` %in% c('10','14'))) 

#make phyloseq for RGM only
filt_rare_RGM2_genus<- rep2_named_phy_genus%>% 
  subset_samples(soilAge %in% ('rgm'))

#order samples
sampR_genus<- sample_data(filt_rare_RGM2_genus) #pull out data from phyloseq

metadata_RGMF<-metadata_RGMF[order(match(metadata_RGMF$`SampleID`, row.names(sampR_genus))), ]

#check both match 
sample_names(filt_rare_RGM2_genus)
metadata_RGMF$SampleID

## Wet Subset Permanova ----
### 5e. Permanova test with wet season data
set.seed(200) ###VERY IMPORTANT, always keep the same

#run permanova
permanova_wet<- adonis2(distance(filt_rare_wet2_genus, method='wunifrac')~treatment*soilAge, data=metadata_wetF, by='terms')
permanova_wet

#pairwise permanova to see which groups are different from each other
permanova_pairwise(distance(filt_rare_wet2_genus, method='wunifrac'), grp=metadata_wetF$trt_soilAge)

# see Plots_18S file for code to make plots

## RGM Subset Permanova----
# 5f. RGM Permanova
set.seed(200)
#permanova
permanova_rgm<- adonis2(distance(filt_rare_RGM2_genus, method='wunifrac')~treatment*`month-collected`, data=metadata_RGMF, by='terms')
permanova_rgm

#pairwise permanova to see which groups are different
permanova_pairwise(distance(filt_rare_RGM2_genus, method='wunifrac'), grp=metadata_RGMF$trt_month)

### see Plots_18S file for code on how to make the plots

## Dry RGM Permanova----
filt_rare_RGM_dry_genus<- subset_samples(filt_rare_RGM2_genus, `month.collected` %in% ('dry'))
metaDryRGM<- metadata_RGMF %>% 
  filter(`month-collected`=='dry') 

#order samples
sampRd_genus<- sample_data(filt_rare_RGM_dry_genus) #pull out data from phyloseq

metaDryRGM<-metaDryRGM[order(match(metaDryRGM$`SampleID`, row.names(sampRd_genus))), ]

set.seed(200)
#permanova
permanova_rgmD<- adonis2(distance(filt_rare_RGM_dry_genus, method='wunifrac')~treatment, data=metaDryRGM, by='terms')
permanova_rgmD

# Homogeneity of dispersions----

##Wet subset: 
wet_betadis_genus <- betadisper(distance(filt_rare_wet2_genus, method = 'wunifrac'), group = metadata_wetF$trt_soilAge, type = 'median') #create betadisper object with dispersion distances  
wet_permutest <- permutest(wet_betadis_genus, permutations = 999) #test for differences in dispersions 
wet_permutest
boxplot(wet_betadis_genus)

#adonis2(dist(wet_betadis$distances) ~ metadata_wetF$trt_soilAge)

#RGM Dry subset: 
dryRGM_betadis_genus<-betadisper(distance(filt_rare_RGM_dry_genus, method='wunifrac'), group=metaDryRGM$treatment, type='median')
permutest(dryRGM_betadis_genus)
boxplot(dryRGM_betadis_genus)

#adonis2(dist(dryRGM_betadis$distances)~metaDryRGM$treatment)


# Differential Abundance----
## Wet RGM Latrine vs Control DA----

#make new phyloseq for rgm that is using original ASV names and not ASV named as a number
RGM2_phy_ASV<- filt_rare_rep2_genus %>% 
  subset_samples(soilAge %in% ('rgm'))

#order samples
sampR<- sample_data(RGM2_phy_ASV) #pull out data from phyloseq

metadata_RGMF<-metadata_RGMF[order(match(metadata_RGMF$`SampleID`, row.names(sampR))), ]

#make season a factor
rgm2_sampdata<- sample_data(RGM2_phy_ASV)
#factor treatment
rgm2_sampdata$treatment<- as.factor(rgm2_sampdata$treatment)
RGM2_phy_ASV@sam_data<- rgm2_sampdata
str(RGM2_phy_ASV@sam_data)

rgmW_rep2_phy<- subset_samples(RGM2_phy_ASV, month.collected=='wet')

# test with just one rep and no RE
rgmWetTreatmentDA<-ancombc2(data = rgmW_rep2_phy, tax_level = "Genus",
                            fix_formula = "treatment", rand_formula = NULL,
                            p_adj_method = "holm", pseudo_sens = TRUE,
                            prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                            group = "treatment", struc_zero = TRUE, neg_lb = TRUE,
                            alpha = 0.05, n_cl = 2, verbose = TRUE,
                            global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                            iter_control = list(tol = 1e-2, max_iter = 20, 
                                                verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),
                            lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

#put primary results in data frame
rgmWetT_prim<-rgmWetTreatmentDA$res

#save it as an rds file
saveRDS(rgmWetT_prim, file='rgmWetT_prim_genus.rds')
rgmWetT_prim_genus<-readRDS('rgmWetT_prim_genus.rds')

#filter for what's significant
rgmWetTSig_genus<-rgmWetT_prim_genus %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T)

#extract taxa from phyloseq
rgmWet_taxa_genus<- data.frame(tax_table(rgmW_rep2_phy))



# Plot log fold change
rgmWetT_DAplot_genus<- rgmWetTSig_genus %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T) %>% 
  dplyr::arrange(desc(lfc_treatmentlatrine)) %>% 
  dplyr::mutate(direct = ifelse(lfc_treatmentlatrine> 0, "Positive LFC", "Negative LFC"))

#make taxon and direction factors
rgmWetT_DAplot_genus$taxon<- factor(rgmWetT_DAplot_genus$taxon, levels=rgmWetT_DAplot_genus$taxon)
rgmWetT_DAplot_genus$direct<- factor(rgmWetT_DAplot_genus$direct, levels = c("Positive LFC", "Negative LFC"))


fig_rgmWetT_genus = rgmWetT_DAplot_genus %>%
  ggplot(aes(x = taxon, y = lfc_treatmentlatrine, fill=direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_treatmentlatrine - se_treatmentlatrine, ymax = lfc_treatmentlatrine + se_treatmentlatrine), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  scale_fill_manual(values=c('purple3', 'cyan3'), name=NULL, labels=c('Positive LFC (more in latrine)','Negative LFC (more in control)'))+
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes") + 
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x= element_text(hjust=1, angle=45))
fig_rgmWetT_genus

#do the test at the phylum level
rgmWetTDA_phylum<-ancombc2(data = rgmW_rep2_phy, tax_level = "Phylum",
                           fix_formula = "treatment", rand_formula = NULL,
                           p_adj_method = "BH", pseudo_sens = TRUE,
                           prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                           group = "treatment", struc_zero = TRUE, neg_lb = TRUE,
                           alpha = 0.05, n_cl = 2, verbose = TRUE,
                           global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                           iter_control = list(tol = 1e-2, max_iter = 20, 
                                               verbose = TRUE),
                           em_control = list(tol = 1e-5, max_iter = 100),
                           lme_control = lme4::lmerControl(),
                           mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

rgmWet_phylum_prim_genus<- rgmWetTDA_phylum$res

rgmWet_phylumSig<- rgmWet_phylum_prim_genus %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T)

# Plot log fold change
rgmWetphylum_DAplot<- rgmWet_phylumSig %>% 
  dplyr::arrange(desc(lfc_treatmentlatrine)) %>% 
  dplyr::mutate(direct = ifelse(lfc_treatmentlatrine> 0, "Positive LFC", "Negative LFC"))

#make taxon and direction factors
rgmWetphylum_DAplot$taxon<- factor(rgmWetphylum_DAplot$taxon, levels=rgmWetphylum_DAplot$taxon)
rgmWetphylum_DAplot$direct<- factor(rgmWetphylum_DAplot$direct, levels = c("Positive LFC", "Negative LFC"))

fig_rgmWetphylum = rgmWetphylum_DAplot %>%
  ggplot(aes(x = taxon, y = lfc_treatmentlatrine, fill=direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_treatmentlatrine - se_treatmentlatrine, ymax = lfc_treatmentlatrine + se_treatmentlatrine), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(hjust=1, angle=90))
fig_rgmWetphylum

## Dry RGM Latrine vs Control----

#got the phyloseq from the latrine by season step

#make treatment a factor
rgm2_sampdata$treatment<- as.factor(rgm2_sampdata$treatment)
RGM2_phy_ASV@sam_data<- rgm2_sampdata
str(RGM2_phy_ASV@sam_data)

#make it just for the dry samples
rgmD_rep2_phy<- subset_samples(RGM2_phy_ASV, month.collected=='dry')

# test with just rep 2 and no RE
rgmDryTreatmentDA<-ancombc2(data = rgmD_rep2_phy, tax_level = "Genus",
                            fix_formula = "treatment", rand_formula = NULL,
                            p_adj_method = "BH", pseudo_sens = TRUE,
                            prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                            group = "treatment", struc_zero = TRUE, neg_lb = TRUE,
                            alpha = 0.05, n_cl = 2, verbose = TRUE,
                            global = TRUE, pairwise = F, dunnet = F, trend = F,
                            iter_control = list(tol = 1e-2, max_iter = 20, 
                                                verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),
                            lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

#put primary results in data frame
rgmDryT_prim_genus<-rgmDryTreatmentDA$res

#save it as an rds file
saveRDS(rgmDryT_prim_genus, file='rgmDryT_prim_genus.rds')
rgmDryT_prim_genus<-readRDS('rgmDryT_prim_genus.rds')

#filter for what's significant
rgmDryTSig_genus<-rgmDryT_prim_genus %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T)

#extract taxa from phyloseq
rgmDry_taxa_genus<- data.frame(tax_table(rgmD_rep2_phy))


# Plot log fold change
rgmDryT_DAplot_genus<- rgmDryTSig_genus %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T) %>% 
  dplyr::arrange(desc(lfc_treatmentlatrine)) %>% 
  dplyr::mutate(direct = ifelse(lfc_treatmentlatrine> 0, "Positive LFC", "Negative LFC"))

#make taxon and direction factors
rgmDryT_DAplot_genus$taxon<- factor(rgmDryT_DAplot_genus$taxon, levels=rgmDryT_DAplot_genus$taxon)
rgmDryT_DAplot_genus$direct<- factor(rgmDryT_DAplot_genus$direct, levels = c("Positive LFC", "Negative LFC"))


fig_rgmDryT_genus = rgmDryT_DAplot_genus %>%
  ggplot(aes(x = taxon, y = lfc_treatmentlatrine, fill=direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_treatmentlatrine - se_treatmentlatrine, ymax = lfc_treatmentlatrine + se_treatmentlatrine), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Dry Latrine vs Control") + 
  scale_color_discrete(name = NULL) +
  scale_fill_manual(values=c( 'cyan3'), name=NULL, labels=c('Negative LFC (more in control)'))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x= element_text(hjust=1, angle=45))
fig_rgmDryT_genus

#do the test at the phylum level
rgmDryTDA_phylum<-ancombc2(data = rgmD_rep2_phy, tax_level = "Phylum",
                           fix_formula = "treatment", rand_formula = NULL,
                           p_adj_method = "BH", pseudo_sens = TRUE,
                           prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                           group = "treatment", struc_zero = TRUE, neg_lb = TRUE,
                           alpha = 0.05, n_cl = 2, verbose = TRUE,
                           global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                           iter_control = list(tol = 1e-2, max_iter = 20, 
                                               verbose = TRUE),
                           em_control = list(tol = 1e-5, max_iter = 100),
                           lme_control = lme4::lmerControl(),
                           mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

rgmDry_phylum_prim<- rgmDryTDA_phylum$res


rgmDry_phylumSig<- rgmDry_phylum_prim %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T)



## LIA vs RGM Control but using only 4 rgm locations that we chose based on location and availability----
##51,60,56,58
#get the phyloseq with regular asv names
wet2_phy_ASV<- filt_rare_rep2_genus%>% 
  subset_samples(month.collected %in% ('wet')) %>% 
  subset_samples(replicate==1) 

## make our test variables factors
wet2_sampdata<- sample_data(wet2_phy_ASV)
wet2_sampdata$soilAge<- as.factor(wet2_sampdata$soilAge)
wet2_phy_ASV@sam_data<- wet2_sampdata
str(wet2_phy_ASV@sam_data)

#select only the 4 rgm samples we want (as well as the LIA ones)
soilAgeCDAphy<- wet2_phy_ASV %>% 
  subset_samples(latrine %in% c('L51','L56','L58','L60','L100','L101','L102','L104')) %>% 
  subset_samples(treatment=='control')

#run the test
soilAgeCDA<-ancombc2(data = soilAgeCDAphy, tax_level = "Genus",
                     fix_formula = "soilAge", rand_formula =NULL,
                     p_adj_method = "holm", pseudo_sens = TRUE,
                     prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                     group = "soilAge", struc_zero = T, neg_lb = T,
                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                     global = TRUE, pairwise = F, dunnet = F, trend = F,
                     iter_control = list(tol = 1e-2, max_iter = 20, 
                                         verbose = TRUE),
                     em_control = list(tol = 1e-5, max_iter = 100),
                     lme_control = lme4::lmerControl(),
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

soilAgeC_prim_genus<- soilAgeCDA$res

saveRDS(soilAgeC_prim_genus, file='soilAgeC_prim_genus.rds')
soilAgeC_prim_genus<-readRDS('soilAgeC_prim_genus.rds')

soilAgeC_sig_genus<- soilAgeC_prim_genus %>% 
  filter(q_soilAgergm<.05 & passed_ss_soilAgergm==T)

soilAgeC_zero<- soilAgeCDA$zero_ind
soilAgeC_zeroLIA<- soilAgeC_zero %>% 
  filter(`structural_zero (soilAge = lia)`==T & `structural_zero (soilAge = rgm)`==F)
soilAgeC_zeroRGM<- soilAgeC_zero %>% 
  filter(`structural_zero (soilAge = lia)`==F & `structural_zero (soilAge = rgm)`==T)

#### do a test at the phylum level
#get LIA samples
wet2_sampdata$treatment<- as.factor(wet2_sampdata$treatment)
wet2_phy_ASV@sam_data<- wet2_sampdata

LIA_phylum_DA<- wet2_phy_ASV %>% 
  subset_samples(soilAge== 'lia')

#run the test
soilAgeDAPhylum<-ancombc2(data = LIA_phylum_DA, tax_level = "Phylum",
                          fix_formula = "treatment", rand_formula =NULL,
                          p_adj_method = "BH", pseudo_sens = TRUE,
                          prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                          group = "treatment", struc_zero = T, neg_lb = T,
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          global = TRUE, pairwise = TRUE, dunnet = F, trend = F,
                          iter_control = list(tol = 1e-2, max_iter = 20, 
                                              verbose = TRUE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))
soilAgeDAPhylum_pair<- soilAgeDAPhylum$res
soilAgeDAPhylumSig<- soilAgeDAPhylum_pair %>% 
  filter(q_treatmentlatrine<.05 & passed_ss_treatmentlatrine==T)

## Plot all 3 phylum tests in one figure----

#first run the 3 different phylum level tests under the Wet L vs C, Dry L vs C, and LIA L vs C

#######first process the wet season comparison
#select just the columns of interest
wet_phyl_fig <- rgmWet_phylum_prim %>%
  dplyr::select(taxon, contains('latrine')) 

#round LFC, choose only taxa that are significant and passed sensitivity test, make data tidy
wet_phyl_fig_lfc <- wet_phyl_fig %>%
  dplyr::filter(diff_treatmentlatrine == 1 & passed_ss_treatmentlatrine==T) %>%
  dplyr::mutate(lfc1 = ifelse(diff_treatmentlatrine == 1, 
                              round(lfc_treatmentlatrine, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon) %>% 
  dplyr::select(group, value, taxon)

# recode the group so instead of lfc1 it says what the comparison is
wet_phyl_fig_lfc$group <- dplyr::recode(wet_phyl_fig_lfc$group, 
                                        `lfc1` = "RGM Wet Season")

###### second process dry season comparison
dry_phyl_fig<- rgmDry_phylum_prim %>% 
  dplyr::select(taxon, contains('latrine')) 

#round LFC, choose only taxa that are significant and passed sensitivity test, make data tidy
dry_phyl_fig_lfc <- dry_phyl_fig %>%
  dplyr::filter(diff_treatmentlatrine == 1 & passed_ss_treatmentlatrine==T) %>%
  dplyr::mutate(lfc1 = ifelse(diff_treatmentlatrine == 1, 
                              round(lfc_treatmentlatrine, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon) %>% 
  dplyr::select(group, value, taxon)

# recode the group so instead of lfc1 it says what the comparison is
dry_phyl_fig_lfc$group <- dplyr::recode(dry_phyl_fig_lfc$group, 
                                        `lfc1` = "RGM Dry Season")

##### third process the LIA comparison
lia_phyl_fig<- soilAgeDAPhylum_pair %>% 
  dplyr::select(taxon, contains('latrine'))

#round LFC, choose only taxa that are significant and passed sensitivity test, make data tidy
lia_phyl_fig_lfc <- lia_phyl_fig %>%
  dplyr::filter(diff_treatmentlatrine == 1 & passed_ss_treatmentlatrine==T) %>%
  dplyr::mutate(lfc1 = ifelse(diff_treatmentlatrine == 1, 
                              round(lfc_treatmentlatrine, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon) %>% 
  dplyr::select(group, value, taxon)

# recode the group so instead of lfc1 it says what the comparison is
lia_phyl_fig_lfc$group <- dplyr::recode(lia_phyl_fig_lfc$group, 
                                        `lfc1` = "LIA Wet Season")

#join all the comparisons together
phyl_fig<-full_join(lia_phyl_fig_lfc, wet_phyl_fig_lfc)
phyl_fig<- full_join(dry_phyl_fig_lfc, phyl_fig)

# make the figure
lo = floor(min(phyl_fig$value))
up = ceiling(max(phyl_fig$value))

fig_phyl = phyl_fig %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "cyan3", high = "purple3", mid = "white", 
                       na.value = "white", midpoint = 0, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = 'Phylum LFC Latrine-Control') +
  theme_classic() +
  theme(axis.text.x=element_text())
fig_phyl
