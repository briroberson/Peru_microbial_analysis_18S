library(mirlyn)

#rarefy data
mirl_object<- mirl(final_filtered_phy, libsize=3705, set.seed=200, trimOTUs=T, replace=F, rep=5)

#make an empty object to put the ASV tables in
mirl_otu <- vector("list", length(mirl_object))

#extract otu tables from each rarefied phyloseq and add to the empty object above
for (i in 1:length(mirl_object)){
  colnames(mirl_object[[i]]@otu_table) <- paste0(colnames(mirl_object[[i]]@otu_table))
  (mirl_otu[[i]] <- mirl_object[[i]]@otu_table)
}



#make metadata file with the correct samples (remove ones dropped during filtering)
sample_id<- data.frame(final_filtered_phy@sam_data) 
sample_id$Samples<- row.names(sample_id)
sample_id<- sample_id %>% 
  filter(!Samples %in% c(5, 21, 127))

#get just the sample names
sample_id <- sample_id$Samples

#make empty list for each sample
average_counts <- vector("list", length(sample_id))

#give how many reps you will do
rep<-1:5
#make empty list to hold however many dataframes
iter_list<- vector('list', length(rep))

#rewrite loop to select columns from each rep, then average them and put them in new otu table
for (i in 1:length(sample_id) ){
  for (j in rep){
    iter_list[[j]]<-dplyr::select(as.data.frame(mirl_otu[[j]]),i) #this selects each individual iteration's otu table and 
    iter_list[[j]]$ASVname<- row.names(iter_list[[j]]) #this makes a column with asv names
  }

  sample_df<- reduce(iter_list[rep], full_join, by='ASVname') #this combines them into one dataframe
  sample_df[is.na(sample_df)]<-0 #make NAs into 0s
  row.names(sample_df)<- sample_df$ASVname #make row names the ASV names
  sample_df<- sample_df[,c(1, 3:(1+length(rep)))] #remove the ASV name column
  sample_average <- data.frame(rowMeans(sample_df)) #calculate the mean of each row (which is the avg abundance on each ASV across iterations)
  colnames(sample_average) <- sample_id[[i]] #make the column name the sample number
  average_counts[[i]] <- sample_average #put into list which has an element for each sample
}
  
#combine each list element (each sample) into one data frame
average_count_df <- do.call(cbind, average_counts)
write.csv(average_count_df, file='D:\\Soil\\16S\\5rep_avgd_data.csv')

######## do some checks
#check that they all have the rarefied number of ASVs
colSums(average_count_df)

#is this close to the number for the whole data frame?
sum(iter_list[[1]]$`99`!=0)
sum(average_count_df$`99` !=0)
#relatively close
#not exact because there could be asvs not present here that are present in other iterations

#choose a random ASV in sample 14 15390803410156ccce6a14bb37cff437

#rel abun of this ASV is 1.6 in the average_count_df and the total ASVs in sample 14 is 3705
1.6/3705

#the proportion of this ASV in sample 14 is about 0.0004318489
#find the rel abun of this ASV in the original data frame (this is pre filtering but that shouldn't affect this one)
#you can open asv_all and search for 1539 and manually scroll to find that the rel abun is 7
#find total number of ASVs in sample 14
sum(asv_all$X101)

7/15975 #find proportion of this ASV in original data frame
#about 0.0004381847 so very close to our new data


#add to phyloseq
mirl_phyloseq <- final_filtered_phy
mirl_phyloseq@otu_table@.Data <- as.matrix(average_count_df)

rowSums(mirl_phyloseq@otu_table)==rowSums(average_count_df)


#do for 100 reps
#rarefy data
mirl_object_100<- mirl(final_filtered_phy, libsize=3705, set.seed=200, trimOTUs=T, replace=F, rep=100)

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
  filter(!Samples %in% c(5, 21, 127))


sample_id <- sample_id$Samples

#make empty list for each sample
average_counts_100 <- vector("list", length(sample_id))

#give how many reps you will do
rep_100<-1:100
#make empty list to hold 5 dataframes
iter_list_100<- vector('list', length(rep_100))

#rewrite loop to select columns from each rep, then average them and put them in new otu table
for (i in 1:length(sample_id) ){
  for (j in rep_100){
    iter_list_100[[j]]<-dplyr::select(as.data.frame(mirl_otu_100[[j]]),i) #this selects each individual iteration's otu table and 
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

write.csv(x=average_count_df_100, file="D:\\Soil\\16S\\100rep_averaged_OTUtable.csv")
