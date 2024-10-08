
library(kableExtra)
library(gt)
library(gtsummary)
library(knitr)
library(khroma)
library(adegenet)
library(magrittr)
library(tidyverse)


marker_info<-read_csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Omykiss_Siletz_2023/GT-seq_genotyping_output/marker_info.csv")

# edit ind names and grab numeric a1 and a2 counts
marker_info %<>%
  mutate(a1_count =  as.numeric(substr(a1_count, 3, nchar(a1_count)))) %>%
  mutate(a2_count =  as.numeric(substr(a2_count, 3, nchar(a2_count)))) %>%
  mutate(ind = str_remove(ind, "^\\./")) %>%
  mutate(ind = str_remove(ind, "\\.genos"))

# read the raw genotypes file in to R
genos_0.1<-read_csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Omykiss_Siletz_2023/GT-seq_genotyping_output/OmykissSiletz_GTs_0.1.csv")


# add a field to mark controls
# here controls contained "positive," "negative" in their sample names so used simple pattern matching to create a new column, you can add you own here for known controls (e.g. known winter run steelhead)
genos_0.1 %<>%
  mutate(control = case_when(str_detect(Sample, "positive") ~ "positive",
                             str_detect(Sample, "negative") ~ "negative",
                             TRUE ~ "sample"))

# clean up sample name field
# split the sample name and the adapter sequence (note that replicates will have the same sample name, but we'll keep track with the adapter sequences)

genos_0.1 %<>%
  mutate(adapter = str_extract(Sample, "[ATCG]{6}-[ATCG]{6}")) %>%
  mutate(sample_simple = str_extract(Sample, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")) %>%
  relocate(Sample, sample_simple, adapter)


# great, prep is done, now lets make our first plot: distribution of reads between controls and samples
ggplot()+geom_histogram(data = genos_0.1, aes(x = `On-Target Reads`, fill= control)) + theme_classic()+scale_fill_viridis_d()

ggplot()+geom_histogram(data = genos_0.1[genos_0.1$control=="negative",], aes(x = `On-Target Reads`)) + theme_classic()

ggplot()+geom_histogram(data = genos_0.1, aes(x = `%On-Target`, fill= control)) + theme_classic()+scale_fill_viridis_d()


# here we filter out our known controls and create our next dataset genos_0.11
genos_0.11 <- genos_0.1 %>%
  filter(control == "sample") %>%
  select(-control) #get rid of this column


# sometimes we'll use more than a single replicate, this broke a previous version of this code chunk, for now we'll find the triplicate (or more) samples, and keep the two with greatest on target read depth (ignoring triplicates quadruplicates etc)

genos_0.11 %<>% 
  group_by(sample_simple) %>%
  slice_max(order_by = `On-Target Reads`, n = 2)

#now let's get duplicated samples
dups <- genos_0.11[genos_0.11$sample_simple %in% genos_0.11$sample_simple[duplicated(genos_0.11$sample_simple)],]
dups <- dups[order(dups$sample_simple),]


# next we'll calculate the percent concordance among replicates
# woof I don't see a good way around using a nested for loop here, maybe fix this in the future

dups_genos <- dups[,9:ncol(dups)] 
rep_info <- matrix(ncol=ncol(dups_genos), nrow=nrow(dups_genos)/2)
colnames(rep_info) <- colnames(dups_genos)
for (j in 1:(nrow(dups_genos)/2)) {
  for (i in 1:ncol(dups_genos)) {
    rep_info[j,i] <- sum(dups_genos[(j*2)-1,i]==dups_genos[(j*2),i])
  }
}

geno_concordance <- as.data.frame(as.matrix(rep_info)) %>%
  rowMeans()

rep_data <- as.data.frame(cbind(dups[c(1:length(geno_concordance))*2,1], geno_concordance))
ggplot(data=rep_data)+geom_histogram(aes(x=geno_concordance))+theme_classic()


#this writes a new dataset (0.2) by choosing the samples within duplicates and keeping the one with the highest genotyping success
genos_0.2 <- genos_0.11 %>%
  group_by(sample_simple) %>%
  filter(`On-Target Reads` == max(`On-Target Reads`))


# Plot sex marker counts

sex_marker_info <- marker_info %>%
  filter(str_detect(marker, "SEXY")) %>%
  mutate(called_geno = replace_na(called_geno, "00")) %>%
  mutate(called_geno = case_when(called_geno == "A1HOM" ~ "XX",
                                 called_geno == "HET" ~ "XY",
                                 called_geno == "00" ~ "00"))

ggplot(data = sex_marker_info)+geom_point(aes(a2_count, a1_count, color = called_geno))+scale_colour_viridis_d(name = "Called Sex Genotype")+theme_classic()+xlab("Y-specific probe count")+ylab("Theoretical X chromosome count")+geom_abline(aes(intercept = 0, slope = 0.1))+geom_abline(aes(intercept = 0, slope = 5))+geom_abline(aes(intercept = 0, slope = 10))+geom_abline(aes(intercept = 0, slope = 0.2))+xlim(0,max(c(sex_marker_info$a1_count, sex_marker_info$a2_count)))+ylim(0,max(c(sex_marker_info$a1_count, sex_marker_info$a2_count)))+geom_abline(aes(intercept = 0, slope = 1), color = "darkred")

kable(sex_marker_info %>% count(called_geno), caption = "Called Sex Ratio")


# First let's collect the number of on-target reads for each sample in the sex_marker_info dataframe
sex_marker_info$sample <- str_extract(sex_marker_info$ind, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")
sex_marker_info$adapter <- str_extract(sex_marker_info$ind, "[ATCG]{6}-[ATCG]{6}") 

sex_marker_info %<>%
  left_join(dplyr::select(genos_0.1, sample_simple, adapter, `On-Target Reads`), by = c("sample" = "sample_simple" , "adapter" = "adapter"))


# now here's the tricky part. we're going to do a first pass filter to eliminate any samples that coldn't possibly be males, then estimate the theoretical number of X chromosome counts for possible males
#first get the estimate
Xmale_proportion_error <- sex_marker_info %>%
  filter(a2_count > 1) %>% #this value is set to 1 if the actual count is zero in the omysexscript, so here we're are eliminating all individuals where it woulbe extremely unlikely that they are males 
  summarise(Xmale_proportion = mean(a2_count/`On-Target Reads`))

#then use the estimate to correct the data
sex_marker_info %<>%
  mutate(XY_count_estimate = Xmale_proportion_error[[1]]*`On-Target Reads`) %>%
  mutate(corrected_ratio = XY_count_estimate/a2_count) %>%
  mutate(corrected_sex_genotype = case_when(corrected_ratio > 10 ~ "XX",
                                            corrected_ratio <= 0.1 ~ "XY",
                                            corrected_ratio <= 0.2 ~ "00",
                                            corrected_ratio <= 5 ~ "XY",
                                            TRUE ~ "00"
  )) #note these ratios are taken directly from the omysex script but this seems like another bug, shouldn't XY < 0.2 be NA not XY???


ggplot(data = sex_marker_info)+geom_point(aes(a2_count, XY_count_estimate, color = corrected_sex_genotype), alpha = 0.7)+scale_colour_viridis_d(name = "Corrected Called\nSex Genotype")+theme_classic()+xlab("Y-specific probe count")+ylab("Corrected Theoretical X chromosome count")+geom_abline(aes(intercept = 0, slope = 0.1))+geom_abline(aes(intercept = 0, slope = 5))+geom_abline(aes(intercept = 0, slope = 10))+geom_abline(aes(intercept = 0, slope = 0.2))+xlim(0,max(c(sex_marker_info$XY_count_estimate, sex_marker_info$a2_count)))+ylim(0,max(c(sex_marker_info$XY_count_estimate, sex_marker_info$a2_count)))+geom_abline(aes(intercept = 0, slope = 1), color = "darkred")

kable(sex_marker_info %>% 
        group_by(corrected_sex_genotype) %>%
        summarise(count = n(), mean_depth = mean(`On-Target Reads`)))


# run this command for Omy
genos_0.2 %<>%
  left_join(select(sex_marker_info, sample, adapter, corrected_sex_genotype), by = c("sample_simple"="sample", "adapter"="adapter")) %>%
  mutate(OmyY1_2SEXY = corrected_sex_genotype) %>%
  select(-corrected_sex_genotype)

ggplot(genos_0.2)+geom_histogram(aes(x=IFI))+geom_vline(aes(xintercept= 2.5), color="red")+theme_classic()

ggplot(genos_0.2)+geom_histogram(aes(x=`%GT`))+geom_vline(aes(xintercept= 90), color="red")+theme_classic()

##column numbers differ from example script##
##remove column 9##
nine<-c(9)

genos_0.3<-genos_0.2[, -nine]

genos_0.2<-genos_0.3
rm(genos_0.3)

missingness <- (colSums(genos_0.2[,c(9:ncol(genos_0.2))] == "00" | genos_0.2[,c(8:(ncol(genos_0.2)-1))] == "0"))/nrow(genos_0.2) #warning hardcoding: "[,8:398]" is hardcoded to work on the example script using the Omy panel with 390 markers, these values will need to be changed to reflect the genotype columns of the genos r object that YOU are running. This excludes columns with metadata and genotyping results such as "sample name" "ifi" "on-target reads" etc

missing <- as.data.frame(missingness)
missing$marker <- row.names(missing)
ggplot(missing) + geom_histogram(aes(x=missingness))+geom_vline(aes(xintercept= 0.2), color="red")+geom_vline(aes(xintercept= 0.1), color="blue")+theme_classic()+xlab("missingness (loci)")

# collect bad markers
very_bad_markers <- missing[missing$missingness>0.5, 2]
print(paste(length(very_bad_markers), "markers with > 50% missing data"))


genos_0.3 <- genos_0.2 %>%
  dplyr::select(-one_of(very_bad_markers))


IFI <- marker_info %>%
  filter(marker %in% colnames(genos_0.3)) %>%
  group_by(ind) %>%
  summarize(back_count = sum(a1_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            hom_ct = sum(a1_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            ifi2 = (back_count/hom_ct)*100)


IFI$sample <- str_extract(IFI$ind, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")
IFI$adapter <- str_extract(IFI$ind, "[ATCG]{6}-[ATCG]{6}") 


genos_0.3 <- genos_0.3 %>%
  left_join(select(IFI, sample, adapter, ifi2), by = c("sample_simple" = "sample", "adapter" = "adapter")) %>%
  mutate(IFI = ifi2) %>%
  select(-one_of("ifi2"))

genos_0.3 %<>% #dont forget to change to the correct dataset
  rowwise %>%
  mutate(`%GT` = 100*(1-sum((cur_data() == "00"), na.rm = TRUE)/(ncol(.)-8))) %>% #change number
  ungroup() #always ungroup after rowwise

missingness3 <- (colSums(genos_0.3[,c(9:(ncol(genos_0.3)))] == "00" | genos_0.3[,c(9:(ncol(genos_0.3)))] == "0"))/nrow(genos_0.3) 
missing3 <- as.data.frame(missingness3)
missing3$marker <- row.names(missing3)


# collect bad markers
bad_markers <- missing3[missing3$missingness3>0.2, 2]
print(paste(length(bad_markers), "markers with > 20% missing data"))


genos_0.4 <- genos_0.3 %>%
  dplyr::select(-one_of(bad_markers))

IFI <- marker_info %>%
  filter(marker %in% colnames(genos_0.4)) %>%
  group_by(ind) %>%
  summarize(back_count = sum(a1_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            hom_ct = sum(a1_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            ifi2 = (back_count/hom_ct)*100)

IFI$sample <- str_extract(IFI$ind, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")
IFI$adapter <- str_extract(IFI$ind, "[ATCG]{6}-[ATCG]{6}") 


genos_0.4 <- genos_0.4 %>%
  left_join(select(IFI, sample, adapter, ifi2), by = c("sample_simple" = "sample", "adapter" = "adapter")) %>%
  mutate(IFI = ifi2) %>%
  select(-one_of("ifi2"))


#get marker names of markers with 0.1 > missingness > 0.2
miss0.1 <- missing3[missing3$missingness3 > 0.1,]
miss_mod <- miss0.1[miss0.1$missingness3 < 0.2, 2]


hets <- filter(marker_info, called_geno == "HET" | is.na(called_geno))


models <- hets %>%
  filter(marker %in% colnames(genos_0.4)) %>%
  filter(is.na(a1_count) == FALSE & is.na(a2_count) == FALSE) %>%
  group_by(marker) %>%
  group_map(~ lm(a1_count ~ a2_count, data= .))


# Apply coef to each model and return a list of allele count ratios
lms <- lapply(models, coef)
ggplot()+geom_histogram(aes(x = sapply(lms,`[`,2)))+theme_classic()+ggtitle("allele ratios for all NA and HET calls")+geom_vline(aes(xintercept = 1.5), color = "red", linetype = 2)+geom_vline(aes(xintercept = (2/3)), color = "red", linetype = 2)+xlab("allele ratio (a1/a2)")+geom_vline(aes(xintercept = 1), color = "black")


#list of p-values
lms_anova <- lapply(models, summary)


# collect info about each bad model
paralog_possible <- which(abs(sapply(lms,`[`,2)) > 1.5) #bad because a positively skewed allele ratio
paralog_possible2 <- which(abs(sapply(lms,`[`,2)) < (2/3)) # bad because a negative skewed allele ratio

paralog_possible3 <- which(sapply(lms_anova, function(x) x$coefficients[,4][2])> 0.01) # bad because too much variance in allele ratio, even if mean ratio is 1

paralog_possible <- c(paralog_possible, paralog_possible2, paralog_possible3)


plots <- marker_info %>%
  filter(marker %in% colnames(genos_0.4)) %>%
  filter(is.na(a1_count) == FALSE & is.na(a2_count) == FALSE) %>%
  group_by(marker) %>%
  do(plots=ggplot(data=.)+geom_point(aes(a1_count, a2_count, color = called_geno))+theme_classic()+geom_abline(aes(slope=1, intercept=0))+geom_abline(aes(slope = 10, intercept=0), color = "green")+geom_abline(aes(slope = 0.1, intercept=0), color = "red")+geom_abline(aes(slope = 0.2, intercept=0), color = "blue")+geom_abline(aes(slope = 5, intercept=0), color = "blue")+coord_equal(ratio=1)+geom_abline(slope = -1, intercept = 10)+ggtitle(unique(.$marker)))


mod_bad_plot_index <- which(plots$marker %in% miss_mod)
paralog_possible <- c(mod_bad_plot_index, paralog_possible)

plots$plots[[paralog_possible[1]]]


to_filt <- c("Omy_CRBF1-1", "Omy_RAD13073-16", "Omy_RAD69583-33", "Omy_RAD74691-49", "Omy_cd59b-112", "Omy_hus1-52", "OmyR40319Pearse", "Omy_RAD27740-55", "Omy_RAD2976-26", "Omy_RAD30243-74") # here list your bad marker names


genos_0.5 <- genos_0.4 %>%
  dplyr::select(-one_of(to_filt))


genos_1.0 <- genos_0.5 %>% 
  select_if(~ length(unique(.)) > 1)

##'related' r package not on CRAN anymore##
##Going off SOP to evaluate duplicates using genetic distance rather than the program Coancestry##

library(poppr)
library(spaa)


##remove "dots" from column names##

genos_2.0 <- genos_1.0
colnames(genos_2.0) <- gsub("\\.", "_", colnames(genos_2.0))


##convert only genotypes to matrix##
GT <- as.matrix(genos_2.0[,c(9:(ncol(genos_2.0)-1))])

##add sample names to rows##
row.names(GT) <- genos_2.0$Sample

##convert matrix to genind##
GT_genind<- df2genind(GT, sep ="", ploidy=2,NA.char = "0")

##write genalex csv## 
genind2genalex(GT_genind, filename = "Omykiss Siletz 2022 genalex.csv")


##read in genalex file##
GT_poppr<-read.genalex("C:/Users/olsenk2/Desktop/Omykiss_Siletz_2023/filtering_record_filtered_genotypes/Omykiss Siletz 2022 genalex.csv")


##genetic distances among fish including self comparisons##
GT_dist<-bitwise.dist(GT_poppr, percent = FALSE)
GT_list<-dist2list(GT_dist)

##remove "0" distances comparing same samples##
same<-which(GT_list$col == GT_list$row)
g_dist<-GT_list[-same, ]


range(g_dist$value)

##genetic distances range 59-162, so no duplicates##

genos_2.0<-genos_1.0

ggplot(genos_2.0)+geom_density(aes(x=`On-Target Reads`))+geom_vline(aes(xintercept=median(`On-Target Reads`)), color = "red") +theme_classic()

ggplot(genos_2.0)+geom_density(aes(x=`%On-Target`))+geom_vline(aes(xintercept=median(`%On-Target`)), color = "red") +theme_classic()

##depths##
marker_info %>%
  filter(marker %in% colnames(genos_2.0)) %>%
  filter(ind %in% genos_2.0$Sample) %>%
  mutate(sumdepth=a1_count+a2_count) %>%
  summarise(mean=mean(sumdepth, na.rm = TRUE), median=median(sumdepth, na.rm = TRUE), sd=sd(sumdepth, na.rm = TRUE))

##mean depth per locus per individual##
marker_info %>%
  filter(marker %in% colnames(genos_2.0)) %>%
  filter(ind %in% genos_2.0$Sample) %>%
  mutate(sumdepth=a1_count+a2_count) %>%
  ggplot + aes(x=sumdepth)+geom_histogram()+theme_classic()+xlab("Mean Depth Per Locus Per Individual")


marker_info %>%
  filter(marker %in% colnames(genos_2.0)) %>%
  filter(ind %in% genos_2.0$Sample) %>%
  mutate(sumdepth=a1_count+a2_count) %>%
  ggplot + aes(x=sumdepth)+geom_histogram()+theme_classic()+xlab("Mean Depth Per Locus Per Individual")+xlim(0,1000)+ggtitle("inset of plot above from 0 to 1000 reads")


save(genos_2.0, file ="Omykiss_Siletz_2022_filtered.R")
save(GT_genind, file= "Omykiss_Siletz_2022_filtered_genind.R")

write.csv(genos_2.0, file = "Omykiss_Siletz_2022_filtered_genos_2.0.csv")

