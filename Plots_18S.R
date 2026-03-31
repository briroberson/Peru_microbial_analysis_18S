#save the theme so that it can be used for all the plots by just typing '+theme
# instead of having to include all of these lines
theme<-theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title.y = element_text(face="bold", size = 18), 
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 18, face = "bold",color = "black"),
      axis.text.x=element_text(size=16),
      plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

library(plyr)
library(heatmaply)
library(gplots)
library(dendextend)
library(colorspace)
library(scales)
library(reshape)
library(funrar)
library(mctoolsr)
library(rvg)
library(officer)


## Elevation plots----
#requires that you load in the metadata and phyloseq, then 
#run step 4a. Alpha diversity calculations and then the first line
# of 4c & 4d to get the metadata files 

# Note: this uses both replicates (three in some RGM cases)
#Wet subset Richness----
wetrichel<-ggplot(metadata_wet, aes(x=elevation, y=Observed, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_manual(values=c('cyan3','purple3'))+
  labs(title='18S Wet Richness')+
  theme+
  theme_bw()
wetrichel
#export the plot to a powerpoint to edit
fig_dml<- rvg::dml(ggobj = wetrichel)

officer::read_pptx() %>%
  # add slide 
officer::add_slide() %>%
  # specify object and location of object 
officer::ph_with(fig_dml, ph_location()) %>%
  # export slide 
base::print(
  target = "F:\\Research\\Plots\\wetrichel.pptx") #specify file path. powerpoint should not be made already, this makes it for you


#Wet subset Shannon----
ggplot(metadata_wet, aes(x=elevation, y=Shannon, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_manual(values=c('cyan3','purple3'))+
  theme_bw()+
  labs(title='18S Wet Shannon')

#Wet subset Inv Simpson----
ggplot(metadata_wet, aes(x=elevation, y=InvSimpson, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_manual(values=c('cyan3','purple3'))+
  labs(title='Wet Inv Simpson')

#Wet subset Pielou----
ggplot(metadata_wet, aes(x=elevation, y=Pielou, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_manual(values=c('cyan3','purple3'))+
  labs(title='Wet Pielou')

# RGM subset Richness----
ggplot(metadata_RGM, aes(x=elevation, y=Observed, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+ 
  scale_color_manual(values=c('cyan3','purple3'))+
  theme_bw()+
  labs(title='RGM Richness')

# RGM subset Shannon----
ggplot(metadata_RGM, aes(x=elevation, y=Shannon, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_manual(values=c('cyan3','purple3'))+
  theme_bw()+
  labs(title='RGM Shannon')

# RGM subset Inv Simpson----
ggplot(metadata_RGM, aes(x=elevation, y=InvSimpson, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_manual(values=c('cyan3','purple3'))+
  theme_bw()+
  labs(title='RGM Inv Simpson')


# RGM Dry Richness----
ggplot(metaDryRGM_both, aes(x=elevation, y=Observed, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+ 
  scale_color_manual(values=c('cyan3','purple3'))+
  theme_bw()+
  labs(title='RGM Dry Richness')

# RGM Dry Shannon----
ggplot(metaDryRGM_both, aes(x=elevation, y=Shannon, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_manual(values=c('cyan3','purple3'))+
  theme_bw()+
  labs(title='RGM Dry Shannon')

# RGM Dry Inv Simpson----
ggplot(metaDryRGM_both, aes(x=elevation, y=InvSimpson, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_manual(values=c('cyan3','purple3'))+
  theme_bw()+
  labs(title='RGM Dry Inv Simpson')

# RGM Dry Pielou----
ggplot(metaDryRGM_both, aes(x=elevation, y=InvSimpson, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_manual(values=c('cyan3','purple3'))+
  theme_bw()+
  labs(title='RGM Dry Pielou')




## Alpha diversity plot----
## Richness ----
wetrich<-ggplot(metadata_wet, aes(treatment, Observed)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) + #adds boxplot
  geom_point(size = 3, aes(color=elevation), alpha = .7) + #adds the individual points
  labs(x = NULL, y = "ASV Richness", title = "a) 18S  Wet Season Alpha Diversity") +
  scale_fill_manual(values=c('cyan3','purple3'), guide='none')+ #colors the two different treatments
  scale_color_gradient(low='lightgray', high='black')+ #colors elevation so low values are lighter
  theme_bw() +
  facet_wrap(~soilAge) #so that it shows two panels, one for each soil age
wetrich
#export the plot to a powerpoint to edit
fig_dml<- rvg::dml(ggobj = wetrich)

officer::read_pptx() %>%
  # add slide 
officer::add_slide() %>%
  # specify object and location of object 
officer::ph_with(fig_dml, ph_location()) %>%
  # export slide 
base::print(
  target = "F:\\Research\\github\\Peru_microbial_analysis_18s\\18s_Plots.pptx")





### Shannon----
ggplot(metadata_wet, aes(treatment, Shannon)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) +
  geom_point(size = 4, aes(color=elevation), alpha = .7) +
  labs(x = NULL, y = "Shannon's Diversity", title = "a) 18S Wet Season Alpha Diversity") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_gradient(low='lightgray', high='black')+
  theme_bw() +
  theme+
  facet_wrap(~soilAge)

### Inv Simpson ----
ggplot(metadata_wet, aes(treatment, InvSimpson)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) +
  geom_point(size = 4, aes(color=elevation), alpha = .7) +
  labs(x = NULL, y = "Inverse Simpson's", title = "a) 18S Wet Season Alpha Diversity") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_gradient(low='lightgray', high='black')+
  theme_bw() +
  facet_wrap(~soilAge)

### Pielou ----
ggplot(metadata_wet, aes(treatment, Pielou)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) +
  geom_point(size = 4, aes(color=elevation), alpha = .7) +
  labs(x = NULL, y = "Pielou", title = "a) 18S Wet Season Alpha Diversity") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_gradient(low='lightgray', high='black')+
  theme_bw() +
  facet_wrap(~soilAge)


### RGM subset----
## Richness ----
ggplot(metadata_RGM, aes(treatment, Observed)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) +
  geom_point(size = 4,aes(color=elevation), alpha = 0.5) +
  labs(x = NULL, y = "ASV Richness", title = "a) 18S RGM Alpha Diversity") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_gradient(low='lightgray', high='black')+
  #scale_x_discrete(limits = rev(levels(meta16$Description))) +
  theme_bw() +
  facet_wrap(~`month-collected`)

## Shannon ----
ggplot(metadata_RGM, aes(treatment, Shannon)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) +
  geom_point(size = 4,aes(color=elevation), alpha = 0.5) +
  labs(x = NULL, y = "Shannon's Diversity", title = "a) 18S RGM Alpha Diversity") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_gradient(low='lightgray', high='black')+
  #scale_x_discrete(limits = rev(levels(meta16$Description))) +
  theme_bw() +
  facet_wrap(~`month-collected`)

## Inv Simpson----
ggplot(metadata_RGM, aes(treatment, InvSimpson)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) +
  geom_point(size = 4,aes(color=elevation), alpha = 0.5) +
  labs(x = NULL, y = "Inverse Simpson's", title = "a) 18S RGM Alpha Diversity") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_gradient(low='lightgray', high='black')+
  #scale_x_discrete(limits = rev(levels(meta16$Description))) +
  theme_bw() +
  facet_wrap(~`month-collected`)



### RGM Dry----
## Richness ----
ggplot(metaDryRGM_both, aes(treatment, Observed)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) +
  geom_point(size = 4,aes(color=elevation), alpha = 0.5) +
  labs(x = NULL, y = "ASV Richness", title = "a) 18S RGM Dry") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_gradient(low='lightgray', high='black')+
  #scale_x_discrete(limits = rev(levels(meta16$Description))) +
  theme_bw() 

## Shannon ----
ggplot(metaDryRGM_both, aes(treatment, Shannon)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) +
  geom_point(size = 4,aes(color=elevation), alpha = 0.5) +
  labs(x = NULL, y = "Shannon", title = "a) 18S RGM Dry") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_gradient(low='lightgray', high='black')+
  #scale_x_discrete(limits = rev(levels(meta16$Description))) +
  theme_bw() 

## Inv Simpson----
ggplot(metaDryRGM_both, aes(treatment, InvSimpson)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) +
  geom_point(size = 4,aes(color=elevation), alpha = 0.5) +
  labs(x = NULL, y = "Inverse Simpson", title = "a) 18S RGM Dry") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_gradient(low='lightgray', high='black')+
  #scale_x_discrete(limits = rev(levels(meta16$Description))) +
  theme_bw() 

## Pielou ----
ggplot(metaDryRGM_both, aes(treatment, Pielou)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment)) +
  geom_point(size = 4,aes(color=elevation), alpha = 0.5) +
  labs(x = NULL, y = "Pielou", title = "a) 18S RGM Dry") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_gradient(low='lightgray', high='black')+
  #scale_x_discrete(limits = rev(levels(meta16$Description))) +
  theme_bw() 





## PCoA plot----
# requires that you run the NECESSARY steps under Beta Diversity to get the 
# subsetted wet and rgm data files

####### all wet----
#get asv table and transpose
asvW<- as.data.frame(otu_table(filt_rare_wet2))
tasvW <- data.frame(t(asvW), check.names = F)

#calculate the pcoa
pcoaW<-cmdscale(d=distance(filt_rare_wet2, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorW<-as.data.frame(wascores(x = pcoaW$points, w = tasvW))

#add the scores to the metadata
metadata_wetF$axis01<- vegan::scores(pcoaW)[,1]
metadata_wetF$axis02<- vegan::scores(pcoaW)[,2]

#use this function to calculate the hulls
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hullsW <- ddply(metadata_wetF, "trt_soilAge", find_hull)

#plot it
wetbeta<-ggplot(metadata_wetF, aes(axis01, axis02)) +
  geom_polygon(data = micro.hullsW, 
               aes(colour = trt_soilAge, fill = trt_soilAge), alpha = 0.1, show.legend = F) +
 # geom_segment(aes(x=0, xend=V1, y=0, yend=V2), data=spscorW, arrow=arrow())+
  geom_point(size = 3, aes(colour = trt_soilAge)) +
  scale_color_manual(labels=c('LIA Control','RGM Control','LIA Latrine','RGM Latrine'),
                       values=c('cyan4', 'cyan2', 'purple4', 'purple1'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  labs(colour = "Treatment & Soil Age", title='by Soil Age') +
  theme_bw() 
wetbeta

#export the plot to a powerpoint to edit
fig_dml<- rvg::dml(ggobj = wetbeta)

pres<-officer::read_pptx() %>%
  # add slide 
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(fig_dml, ph_location())
print(pres, target="F:\\Research\\github\\Peru_microbial_analysis_18s\\18s_Plots.pptx")

#the geom_segment code is used to add arrows to the plot. currently, it is plotting all
# asvs (the whole spscorW dataframe). to plot specific ones, you can subset the spscorW data frame
# by making a list of the ASVs that you want (from multipatt, simper, distance, etc.)
# and then do spscorW[names of ASVs you want,]

##### just LIA----
#get just lia phyloseq and metadata
filt_lia2<- subset_samples(filt_rare_wet2, soilAge %in% ('lia'))
metalia<- metadata_wetF %>% 
  filter(soilAge=='lia')
#get asv table and transpose it
asvLIA<- as.data.frame(otu_table(filt_lia2))
tasvLIA <- data.frame(t(asvLIA), check.names = F)

#calculate the pcoa
pcoaLIA<-cmdscale(d=distance(filt_lia2, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorLIA<-as.data.frame(wascores(x = pcoaLIA$points, w = tasvLIA))

#add scores to metadata
metalia$axis01<- vegan::scores(pcoaLIA)[,1]
metalia$axis02<- vegan::scores(pcoaLIA)[,2]

#retrieve just the top 10 Simper ASVs coordiantes
simpscor<- spscorLIA[simpLIA_asv,]

#order the taxa to match the coordiantes
simperLIA_taxa<-simperLIA_taxa[row.names(simpscor),]

# make a list that contains the lowest taxanomic assignment for each of the top 10
simpLIA_plottaxa<- c('Annelida_X
Annelida_XX','Turbellaria
Rhabdocoela','Annelida_X
Annelida_XX','Ochromonadales
Ochromonadaceae', 'Unassigned', 'Annelida_X
Annelida_XX','Olpidiales
Olpidiaceae','Annelida_X
Annelida_XX','Embryophyceae_X
Embryophyceae_XX','Pezizomycotina
Leotiomycetes')


#make hulls
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hullsLIA <- ddply(metalia, "treatment", find_hull)

#plot it
PCoALIA<-ggplot(metalia, aes(axis01, axis02)) +
  geom_polygon(data = micro.hullsLIA, 
               aes(colour = treatment, fill = treatment), alpha = 0.1, show.legend = F) +
  geom_point(size = 4, aes(color=treatment)) +
 geom_segment(data=simpscor, aes(x=0, xend=V1, y=0, yend=V2), arrow=arrow())+
 geom_point(data=simpscor, aes(x=V1, y=V2), shape=3)+
 annotate(geom="text", x=simpscor$V1, y=simpscor$V2, label=simpLIA_plottaxa,
           color="black") +
  scale_color_manual(labels=c('Control', 'Latrine'), values=c('cyan3','purple3'))+
  scale_fill_manual(values=c('cyan3','purple3'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  labs(colour = "Treatment", title='LIA', subtitle='') +
  theme_bw()
PCoALIA
#export the plot to a powerpoint to edit
fig_dml<- rvg::dml(ggobj = PCoALIA)

officer::read_pptx() %>%
  # add slide 
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(fig_dml, ph_location()) %>%
  # export slide 
  base::print(
    target = "F:\\Research\\16S_Soil\\Plots\\PCoALIA.pptx")



##### just RGM wet----
filt_wet_rgm2<- subset_samples(filt_rare_wet2, soilAge %in% ('rgm'))
metargmW<- metadata_wetF %>% 
  filter(soilAge=='rgm')


asvWrgm<- as.data.frame(otu_table(filt_wet_rgm2))
tasvWrgm <- data.frame(t(asvWrgm), check.names = F)

#calculate the pcoa
pcoaWrgm<-cmdscale(d=distance(filt_wet_rgm2, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorWrgm<-as.data.frame(wascores(x = pcoaWrgm$points, w = tasvWrgm))

simpscorWrgm<- spscorWrgm[simpWrgm_asv,]

#order the taxa to match the coordiantes
simperWrgm_taxa<-simperWrgm_taxa[row.names(simpscorWrgm),]

# make a list that contains the lowest taxanomic assignment for each of the top 10
simpWrgm_plottaxa<- c('Unassigned','Hypotrichia','Annelida_X
Annelida_XX','Olpidiales
Olpidiaceae','Heterotrichea_X
Blepharismidae','Chelicerata
Arachnida','Cryomonadida','Unassigned','Cryomonadida
Rhogostoma-lineage','Chrysophyceae')


#add the scores to the metadata
metargmW$axis01<- vegan::scores(pcoaWrgm)[,1]
metargmW$axis02<- vegan::scores(pcoaWrgm)[,2]
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hullsRGMW<- ddply(metargmW, "treatment", find_hull)


PCoARGMw<-ggplot(metargmW, aes(axis01, axis02)) +
  geom_polygon(data = micro.hullsRGMW, 
               aes(colour = treatment, fill = treatment), alpha = 0.1, show.legend = F) +
  geom_segment(data=simpscorWrgm, aes(x=0, xend=V1, y=0, yend=V2), arrow=arrow())+
  geom_point(size = 4, aes(color=treatment)) +
  geom_point(data=simpscorWrgm, aes(x=V1, y=V2), shape=3)+
  annotate(geom="text", x=simpscorWrgm$V1, y=simpscorWrgm$V2, label=simpWrgm_plottaxa,
           color="black") +
  scale_color_manual(labels=c('Control', 'Latrine'), values=c('cyan3','purple3'))+
  scale_fill_manual(values=c('cyan3','purple3'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  labs(colour = "Treatment", title='RGM Wet', subtitle='') +
  theme_bw()
PCoARGMw


#export the plot to a powerpoint to edit
fig_dml<- rvg::dml(ggobj = PCoARGMw)

officer::read_pptx() %>%
  # add slide 
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(fig_dml, ph_location()) %>%
  # export slide 
  base::print(
    target = "F:\\Research\\16S_Soil\\Plots\\PCoARGMw.pptx")



###### rgm dry----
filt_rare_RGM_dry<- subset_samples(filt_rare_RGM2, `month.collected` %in% ('dry'))
metaDryRGM<- metadata_RGMF %>% 
  filter(`month-collected`=='dry') 

asvDrgm<- as.data.frame(otu_table(filt_rare_RGM_dry))
tasvDrgm <- data.frame(t(asvDrgm), check.names = F)

#calculate the pcoa
pcoaDrgm<-cmdscale(d=distance(filt_rare_RGM_dry, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorDrgm<-as.data.frame(wascores(x = pcoaDrgm$points, w = tasvDrgm))
simpscorDrgm<- spscorDrgm[simpDrgm_asv,]

#order the taxa to match the coordinates
simperDrgm_taxa<-simperDrgm_taxa[row.names(simpscorDrgm),]

# make a list that contains the lowest taxanomic assignment for each of the top 10
simpDrgm_plottaxa<- c('Maxillopoda','Unassigned','Pezizomycotina Eurotiomycetes','Embryophyceae_XX','Unassigned','Unassigned','Unassigned','Arachnida','Dothideomycetes','Embryophyceae_XX')

metaDryRGM$axis01<- vegan::scores(pcoaDrgm)[,1]
metaDryRGM$axis02<- vegan::scores(pcoaDrgm)[,2]

find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hullsRGMD<- ddply(metaDryRGM, "treatment", find_hull)


PCoArgmD<-ggplot(metaDryRGM, aes(axis01, axis02)) +
  geom_polygon(data = micro.hullsRGMD, 
               aes(colour = treatment, fill = treatment), alpha = 0.1, show.legend = F) +
  geom_point(size = 4, aes(color=treatment)) +
  scale_color_manual(labels=c('Control', 'Latrine'), values=c('cyan3','purple3'))+
  scale_fill_manual(values=c('cyan3','purple3'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  labs(colour = "Treatment", title='RGM Dry', subtitle='') +
  theme_bw() 
PCoArgmD

#export the plot to a powerpoint to edit
fig_dml<- rvg::dml(ggobj = PCoArgmD)

officer::read_pptx() %>%
  # add slide 
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(fig_dml, ph_location()) %>%
  # export slide 
  base::print(
    target = "F:\\Research\\16S_Soil\\Plots\\PCoArgmD.pptx")



##### all RGM----
asvrgm<- as.data.frame(otu_table(filt_rare_RGM2))
tasvrgm <- data.frame(t(asvrgm), check.names = F)

#calculate the pcoa
pcoargm<-cmdscale(d=distance(filt_rare_RGM2, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorrgm<-as.data.frame(wascores(x = pcoargm$points, w = tasvrgm))

metadata_RGMF$axis01<- vegan::scores(pcoargm)[,1]
metadata_RGMF$axis02<- vegan::scores(pcoargm)[,2]

micro.hullsR <- ddply(metadata_RGMF, "trt_month", find_hull)

ggplot(metadata_RGMF, aes(axis01, axis02)) +
  geom_polygon(data = micro.hullsR, 
               aes(colour = trt_month, fill = trt_soilAge), alpha = 0.1, show.legend = F) +
  geom_point(size = 4, aes(color=trt_month)) +
  geom_segment(aes(x=0, y=0, xend=V1, yend=V2), data=spscorrgm, arrow=arrow())+
  scale_color_discrete(labels=c('Dry Control', 'Wet Control', 'Dry Latrine', 'Wet Latrine'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  labs(colour = "Treatment", title='by Season') +
  theme_bw() 




####### PCOA with both replicates labeled----

## LIA----
#get just lia phyloseq and metadata
filt_lia_both<- subset_samples(filt_rare_phy, soilAge %in% ('lia'))
metalia_both<- metadata_factored %>% 
  filter(soilAge=='lia')

samp_lia_both<- sample_data(filt_lia_both) #pull out data from phyloseq

metalia_both<-metalia_both[order(match(metalia_both$`SampleID`, row.names(samp_lia_both))), ]


#get asv table and transpose it
asvLIA<- as.data.frame(otu_table(filt_lia_both))
tasvLIA <- data.frame(t(asvLIA), check.names = F)

#calculate the pcoa
pcoaLIA<-cmdscale(d=distance(filt_lia_both, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorLIA<-as.data.frame(wascores(x = pcoaLIA$points, w = tasvLIA))

#add scores to metadata
metalia_both$axis01<- vegan::scores(pcoaLIA)[,1]
metalia_both$axis02<- vegan::scores(pcoaLIA)[,2]

metalia_both$latrine_rep<- paste(metalia_both$latrine_trt, metalia_both$replicate)
  
#make hulls
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hullsLIA <- ddply(metalia_both, "treatment", find_hull)

#plot it
PCoALIA<-ggplot(metalia_both, aes(axis01, axis02, label=latrine_rep)) +
  geom_polygon(data = micro.hullsLIA, 
               aes(colour = treatment, fill = treatment), alpha = 0.1, show.legend = F) +
  geom_point(size = 4, aes(color=treatment))  +
  geom_text()+
  scale_color_manual(labels=c('Control', 'Latrine'), values=c('cyan3','purple3'))+
  scale_fill_manual(values=c('cyan3','purple3'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  labs(colour = "Treatment", title='LIA', subtitle='with both reps') +
  theme_bw()
PCoALIA





##### just RGM wet----
filt_wet_rgm_both<- subset_samples(filt_rare_phy, soilAge %in% ('rgm') & month.collected=='wet')
metargmW_both<- metadata_factored %>% 
  filter(soilAge=='rgm' & `month-collected`=='wet')

samp_rgmW_both<- sample_data(filt_wet_rgm_both) #pull out data from phyloseq

metargmW_both<-metargmW_both[order(match(metargmW_both$`SampleID`, row.names(samp_rgmW_both))), ]


asvWrgm<- as.data.frame(otu_table(filt_wet_rgm_both))
tasvWrgm <- data.frame(t(asvWrgm), check.names = F)

#calculate the pcoa
pcoaWrgm<-cmdscale(d=distance(filt_wet_rgm_both, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorWrgm<-as.data.frame(wascores(x = pcoaWrgm$points, w = tasvWrgm))

#add the scores to the metadata
metargmW_both$axis01<- vegan::scores(pcoaWrgm)[,1]
metargmW_both$axis02<- vegan::scores(pcoaWrgm)[,2]
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hullsRGMW<- ddply(metargmW_both, "treatment", find_hull)

metargmW_both$treatmentabr<- ifelse(metargmW_both$treatment=='control', 'C', 'L')

metargmW_both$latrine_rep<- as.factor(paste(metargmW_both$latrine, metargmW_both$treatmentabr, metargmW_both$replicate))

library(ggnewscale)
library(ggrepel)

PCoARGMw<-ggplot(metargmW_both, aes(axis01, axis02)) +
  geom_polygon(data = micro.hullsRGMW, 
               aes(colour = treatment, fill = treatment), alpha = 0.1, show.legend = F) +
  geom_point(size = 4, aes(color=treatment)) +
  scale_color_manual(labels=c('Control', 'Latrine'), values=c('cyan3','purple3'))+
  ggnewscale::new_scale_color()+
  geom_text_repel(aes(label=latrine_rep, color=as.factor(replicate)))+
  scale_color_manual(values=c('blue','black'))+
  scale_fill_manual(values=c('cyan3','purple3'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  labs(colour = "Treatment", title='RGM Wet', subtitle='both reps') +
  theme_bw()
PCoARGMw


adonis2(distance(filt_wet_rgm_both, method='wunifrac')~treatment, data=metargmW_both, by='terms')


###### rgm dry----
filt_rare_RGM_dry_both<- subset_samples(filt_rare_phy, `month.collected` %in% ('dry') & soilAge=='rgm')
metaDryRGM_both<- metadata_factored %>% 
  filter(`month-collected`=='dry' & soilAge=='rgm') 

samp_rgmD_both<- sample_data(filt_rare_RGM_dry_both) #pull out data from phyloseq

metaDryRGM_both<-metaDryRGM_both[order(match(metaDryRGM_both$`SampleID`, row.names(samp_rgmD_both))), ]

asvDrgm<- as.data.frame(otu_table(filt_rare_RGM_dry_both))
tasvDrgm <- data.frame(t(asvDrgm), check.names = F)

#calculate the pcoa
pcoaDrgm<-cmdscale(d=distance(filt_rare_RGM_dry_both, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorDrgm<-as.data.frame(wascores(x = pcoaDrgm$points, w = tasvDrgm))

metaDryRGM_both$axis01<- vegan::scores(pcoaDrgm)[,1]
metaDryRGM_both$axis02<- vegan::scores(pcoaDrgm)[,2]

find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hullsRGMD<- ddply(metaDryRGM_both, "treatment", find_hull)

metaDryRGM_both$treatmentabr<- ifelse(metaDryRGM_both$treatment=='control', 'C', 'L')

metaDryRGM_both$latrine_rep<- as.factor(paste(metaDryRGM_both$latrine, metaDryRGM_both$treatmentabr, metaDryRGM_both$replicate))

library(ggrepel)
PCoArgmD<-ggplot(metaDryRGM_both, aes(axis01, axis02)) +
  geom_polygon(data = micro.hullsRGMD, 
               aes(colour = treatment, fill = treatment), alpha = 0.1, show.legend = F) +
  geom_point(size = 4, aes(color=treatment)) +
  scale_color_manual(labels=c('Control', 'Latrine'), values=c('cyan3','purple3'))+
  ggnewscale::new_scale_color()+
  geom_text_repel(aes(label=latrine_rep, color=as.factor(replicate)))+
  scale_color_manual(values=c('blue','black'))+
  scale_fill_manual(values=c('cyan3','purple3'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  labs(colour = "Treatment", title='RGM Dry', subtitle='both reps') +
  theme_bw() 
PCoArgmD


# Dendrogram and Heat maps----


############# dendrogram----

#lia dendrogram----
# get the distance matrix
distlia<- distance(filt_lia2, method='wunifrac')
# make a dendrogram using hclust on the distance
dendLia<- as.dendrogram(hclust(distlia, method='ward.D2'))
#color the branches so at the first split, one split is one color and the other split is the other
dendLia<- color_branches(dendLia, k=2, col=c("cyan3","purple3"))
col_labels<- get_leaves_branches_col(dendLia)
col_labels <- col_labels[order(order.dendrogram(dendLia))]
dendLia <- set(dendLia, "labels_cex", 0.8)
dendLia<-place_labels(dendLia, as.character(metalia$latrine_trt))
#this sets the dimensions of the plotting plane so that the dendrogram fits in it
par(mar = c(1,1,1,14))
plot_horiz.dendrogram(dendLia, side=F)

png("F:\\Research\\Plots\\dendlia.png", units='in', width=15, height=5, res=600)


#get the lables in order for use in the heat map
LIAlab<-labels(dendLia)

#export the plot to a powerpoint to edit
fig_dml<- rvg::dml(ggobj = dendLia)

officer::read_pptx() %>%
  # add slide 
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(fig_dml, ph_location()) %>%
  # export slide 
  base::print(
    target = "F:\\Research\\Plots\\dendlia.pptx") #specify file path. powerpoint should not be made already, this makes it for you


######### rgm wet dendrogram----
distrgmW<- distance(filt_wet_rgm2, method='wunifrac')
dendrgmW<- as.dendrogram(hclust(distrgmW, method='ward.D2'))
dendrgmW<- color_branches(dendrgmW, k=2, col=c("cyan3","purple3"))
col_labels<- get_leaves_branches_col(dendrgmW)
col_labels <- col_labels[order(order.dendrogram(dendrgmW))]
dendrgmW <- set(dendrgmW, "labels_cex", 0.8)
dendrgmW<-place_labels(dendrgmW, as.character(metargmW$latrine_trt))
par(mar = c(1,1,1,14))
plot_horiz.dendrogram(dendrgmW, side=F)

rgmWlab<-labels(dendrgmW)

#####rgm dry dendrogram----
distrgmD<- distance(filt_rare_RGM_dry, method='wunifrac')
dendrgmD<- as.dendrogram(hclust(distrgmD, method='ward.D2'))
dendrgmD<- color_branches(dendrgmD, k=3, col=c("cyan3",'cyan3', "purple3"))
col_labels<- get_leaves_branches_col(dendrgmD)
col_labels <- col_labels[order(order.dendrogram(dendrgmD))]
dendrgmD <- set(dendrgmD, "labels_cex", 0.8)
dendrgmD<-place_labels(dendrgmD, as.character(metaDryRGM$latrine_trt))
par(mar = c(1,1,1,14))
plot_horiz.dendrogram(dendrgmD, side=F)

rgmDlab<-labels(dendrgmD)




###############heatmap----

### plot ts heatmap code----
#I had to change the code because the mutate_() function doesn't work anymore and R will
#no longer ignore it and then it wasn't working so I had to remove a ~ that they had inthe code
#but basically I had to manually make this function instead of using the mctoolsr one to make the heatmaps
plot_ts_heatmap = function(tax_table, metadata_map, min_rel_abund, type_header,
                           scale_by = 'all', custom_sample_order,
                           rev_taxa = FALSE, custom_taxa_order, other_label,
                           remove_other = FALSE,
                           colors = c('blue', 'white', 'red')) {
  # group all taxa lower than threshold into other
  lt_thresh = tax_table[rowMeans(tax_table) < min_rel_abund,]
  gt_thresh = tax_table[rowMeans(tax_table) >= min_rel_abund,]
  Other = colSums(lt_thresh)
  sumtax_mod = rbind(gt_thresh, Other = Other)
  # remove other
  if (remove_other) {
    sumtax_mod = sumtax_mod[row.names(sumtax_mod) != 'Other',]
  }
  if (!missing(other_label)) {
    row.names(sumtax_mod)[row.names(sumtax_mod) == 'Other'] = other_label
  }
  # get means
  sumtax_smry = taxa_summary_by_sample_type(sumtax_mod, metadata_map,
                                            type_header, smry_fun = mean)
  sumtax_smry = round(sumtax_smry * 100, 1)
  melted = reshape2::melt(sumtax_smry)
  if (scale_by == 'sample_types') {
    to_plot = dplyr::mutate(dplyr::group_by(melted, "Var2"),
                            scaled =  scales::rescale(value))
  } else if (scale_by == 'taxa') {
    to_plot = dplyr::mutate(dplyr::group_by(melted, "Var1"),
                            scaled =  scales::rescale(value))
  } else if (scale_by == 'all') {
    to_plot = dplyr::mutate(melted, scaled =  scales::rescale(value))
  } else
    stop("scale_by must be one of: 'sample_types', 'taxa' or 'all'.")
  if (!missing(custom_sample_order)) {
    to_plot$Var2 = factor(to_plot$Var2, levels = rev(custom_sample_order))
  }
  if (!missing(custom_taxa_order)) {
    to_plot$Var1 = factor(to_plot$Var1, levels = custom_taxa_order)
  }
  # reverse order of taxa
  if (rev_taxa) {
    to_plot$Var1 = factor(to_plot$Var1, levels = rev(unique(to_plot$Var1)))
  }
  # plot
  # https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/
  p = ggplot2::ggplot(to_plot, ggplot2::aes_string("Var1", "Var2", 
                                                   fill = "scaled")) +
    ggplot2::geom_tile(color = 'black', size = 0.25) +
    ggplot2::scale_fill_gradientn(colours = colors,
                                  values = c(
                                    min(to_plot$scaled),
                                    mean(to_plot$scaled),
                                    max(to_plot$scaled)
                                  )) +
    ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::theme(
      legend.position = 'none',
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1, vjust = 0.25
      ),
      axis.text = ggplot2::element_text(color = 'gray20')
    ) +
    ggplot2::geom_text(data = to_plot, ggplot2::aes_string(label = "value"), 
                       size = 3) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0))
  p
}

#' @title Collapse taxonomy dataframe to character vector
#' @description A quick way to get taxonomy strings for each taxon
#' @param taxonomy_df The dataframe containing taxonomy as loaded by
#'   \code{\link{load_taxa_table}}.
#' @concept Taxonomy-based analyses
#' @examples 
#' collapse_taxonomy(fruits_veggies$taxonomy_loaded)
collapse_taxonomy = function(taxonomy_df) {
  apply(taxonomy_df, 1, function(x)
    paste0(x, collapse = '; '))
}

#########lia----
#glomerate by phylum
glomL <- tax_glom(filt_lia2, taxrank = 'Phylum')
#melt it
datL <- psmelt(glomL)
#make it a character
datL$Phylum <- as.character(datL$Phylum)

#aggregate each phylum for each sample
Phylum_abundanceL <- aggregate(Abundance~Sample+Phylum, datL, FUN=sum)
#cast it so it is long i think?
Phylum_abundanceL <- cast(Phylum_abundanceL, Sample ~ Phylum)

#change row names to match sample
row.names(Phylum_abundanceL)<- Phylum_abundanceL$Sample
#de-select first column
Phylum_abundanceL<- Phylum_abundanceL[, -1]
#make it a data frame bc right now it's a cast df bc we used the cast function
Phylum_abundanceL<- as.data.frame(Phylum_abundanceL)


#calculate relative abundance
Phy_relabL<-make_relative(as.matrix(Phylum_abundanceL))
Phy_relabL<-data.frame(Phy_relabL, check.names=F)
#transpose it
Phy_relabLt<- data.frame((t(Phy_relabL)), check.names=F)

#order labels to match dendrogram by making it a factor
metalia$latrine_trt<- factor(metalia$latrine_trt, levels=LIAlab)

#plot it
#this is the code to save it as an image, run this line then the plot then dev.off()
png("F:\\Research\\Plots\\liaheat2.pdf", units = "in", width = 10, height = 5, res = 600) 

dev.off()
heatlia<-plot_ts_heatmap(Phy_relabLt, metalia, 0.01, "latrine_trt", colors=c('#fcfdbf','#b73779', '#403c3c')) +
  theme(axis.text.y = element_text(size = 12, hjust = 0,
                                   margin = margin(c(0,-1,0,0))),
        axis.text.x = element_text(size = 10, angle = 30, hjust = 1, vjust = 1,
                                   margin = margin(c(-1.5,0,0,0))),
        plot.margin = unit(c(0.1,0.1,0.1,1), "cm"))+
  labs(title='LIA')

heatlia
#the colors part in the first line of the heatmap sets the gradient for the tile colors.
# you can set up to three

#export the plot to a powerpoint to edit
fig_dml<- rvg::dml(ggobj = heatlia)

officer::read_pptx() %>%
  # add slide 
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(fig_dml, ph_location()) %>%
  # export slide 
  base::print(
    target = "F:\\Research\\Plots\\heatlia2.pptx") #specify file path. powerpoint should not be made already, this makes it for you


########for wet rgm----
#aggregate the taxa at the phylum level
glomW <- tax_glom(filt_wet_rgm2, taxrank = 'Phylum')
datW <- psmelt(glomW)
datW$Phylum <- as.character(datW$Phylum)

Phylum_abundanceW <- aggregate(Abundance~Sample+Phylum, datW, FUN=sum)
Phylum_abundanceW <- cast(Phylum_abundanceW, Sample ~ Phylum)

row.names(Phylum_abundanceW)<- Phylum_abundanceW$Sample
Phylum_abundanceW<- Phylum_abundanceW[, -1]
Phylum_abundanceW<- as.data.frame(Phylum_abundanceW)


#calculate relative abundance
Phy_relabW<-make_relative(as.matrix(Phylum_abundanceW))
Phy_relabW<-data.frame(Phy_relabW, check.names=F)
#transpose it
Phy_relabWt<- data.frame((t(Phy_relabW)), check.names=F)

metargmW<- metadata_wetF %>% 
  filter(soilAge=='rgm')

#order latrine names
metargmW$latrine_trt<- factor(metargmW$latrine_trt, levels=rgmWlab)


#plot it
heatrgmw<-plot_ts_heatmap(Phy_relabWt, metargmW, 0.01, "latrine_trt", colors=c('#fcfdbf','#b73779', '#403c3c')) +
  theme(axis.text.y = element_text(size = 12, hjust = 0,
                                    margin = margin(c(0,-1,0,0))),
        axis.text.x = element_text(size = 10, angle = 30, hjust = 1, vjust = 1,
                                   margin = margin(c(-1.5,0,0,0))),
        plot.margin = unit(c(0.1,0.1,0.1,1), "cm"))+
  labs(title='Wet RGM')
heatrgmw

#export the plot to a powerpoint to edit
fig_dml<- rvg::dml(ggobj = heatrgmw)

officer::read_pptx() %>%
  # add slide 
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(fig_dml, ph_location()) %>%
  # export slide 
  base::print(
    target = "F:\\Research\\Plots\\heatrgmw3.pptx") #specify file path. powerpoint should not be made already, this makes it for you


### for dry rgm----
#glomerate it 
glomR <- tax_glom(filt_rare_RGM_dry, taxrank = 'Phylum')
datR <- psmelt(glomR)
datR$Phylum <- as.character(datR$Phylum)

Phylum_abundanceR <- aggregate(Abundance~Sample+Phylum, datR, FUN=sum)
Phylum_abundanceR <- cast(Phylum_abundanceR, Sample ~ Phylum)

row.names(Phylum_abundanceR)<- Phylum_abundanceR$Sample
Phylum_abundanceR<- Phylum_abundanceR[, -1]
Phylum_abundanceR<- as.data.frame(Phylum_abundanceR)


#calculate relative abundance
Phy_relabR<-make_relative(as.matrix(Phylum_abundanceR))
Phy_relabR<-data.frame(Phy_relabR, check.names=F)
#transpose it
Phy_relabRt<- data.frame((t(Phy_relabR)), check.names=F)

metaDryRGM<- metadata_RGMF %>% 
  filter(`month-collected`=='dry') # %>% 
  #filter(!`SampleID` %in% c('11B', '19B', '23', '49', '7')) #drop duplicate samples


metaDryRGM$latrine_trt<- factor(metaDryRGM$latrine_trt, levels=rgmDlab)


#plot it
plot_ts_heatmap(Phy_relabRt, metaDryRGM, 0.01, "latrine_trt",colors=c('#fcfdbf','#b73779', '#403c3c') ) +
  theme(axis.text.y = element_text(size = 12, hjust = 0,
                                   margin = margin(c(0,-1,0,0))),
        axis.text.x = element_text(size = 10, angle = 30, hjust = 1, vjust = 1,
                                   margin = margin(c(-1.5,0,0,0))),
        plot.margin = unit(c(0.1,0.1,0.1,1), "cm")) +
  labs(title='RGM Dry')








