

library(kableExtra)
library(gt)
library(gtsummary)
library(knitr)
library(khroma)
library(adegenet)
library(magrittr)
library(tidyverse)



##Load genotype file, marker info, sampling info##

load("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Omykiss_Siletz_2023/filtering_record_filtered_genotypes/Omykiss_Siletz_2022_filtered_genotypes_2.2.R")
load("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Omykiss_Siletz_2023/filtering_record_filtered_genotypes/Omykiss_Siletz_2022_filtered_genind_2.0.R")




marker_mapping <- readxl::read_xlsx("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Omykiss_Siletz_2023/Omy_GTseq_panel/final_mapping_results.xlsx", sheet = 1)

marker_summary <- marker_mapping %>%
  mutate(marker = str_replace(marker, "Omy(\\d+)", "Chr\\1")) %>% # convert names to match
  filter(marker %in% colnames(genos_2.0)) %>% # exclude filtered loci
  mutate(neutral = if_else(str_detect(`Presumed Type`, 'neutral|Neutral'), "neutral", "adaptive")) %>%
  mutate(run_timing = if_else(str_detect(`Presumed Type`, 'run|Run'), "run_timing", "non-run_timing")) %>%
  dplyr::select(marker, chr, start, 'Presumed Type', neutral, run_timing, Source)


run_timing_loci_names <- marker_mapping %>%
  filter(chr == "28" | CRITFC_chromosome == "28") %>%
  filter(str_detect(`Presumed Type`, 'run|Run')) %>%
  select(marker)

run_timing_loci_names <- str_replace(run_timing_loci_names$marker, "Omy28", "Chr28")



####Grab Run time markers and sex marker##

genos_runtimesex<-genos_2.0[, run_timing_loci_names]


##two of the run time markers (Chr28_11632591 & Omy_RAD47080-54) are not in dataset### possibly dropped out during filtering##

present_markers<-c("sample_simple", "Chr28_11667578", "Omy_GREB1_09", "Chr28_11607954", "Omy_RAD52458-17", "Omy_GREB1_05", "Chr28_11625241", "Chr28_11658853", "Omy_RAD15709-53",
                   "Chr28_11671116", "Chr28_11676622", "Chr28_11683204", "Chr28_11773194", "Chr28_11702210", "OmyY1_2SEXY")

genos_runtimesex<-genos_2.0[, present_markers]


write.csv(genos_runtimesex, file = "Omykiss_Siletz_2022_run_sex_genos.csv")


##Grab Omy 5 resident vs anadromy markers##

##one (OmyR40319) of the 6 Omy-5 markers missing from filtered dataset## ##Was removed because it is monomorphic##



markers<-c("sample_simple", "OmyY1_2SEXY", "OmyR24370Pearse", "OmyR40252Pearse", "OmyR19198Pearse", "OmyR33562Pearse", "OmyR14589Pearse")

genos_omy5<-genos_2.0[, markers]


write.csv(genos_omy5, file = "Omykiss_Siletz_2022_Omy5.csv")






##Summary analyses##
##look at initial panel##

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


# here we filter out our known controls and create our next dataset genos_0.11
genos_0.11 <- genos_0.1 %>%
  filter(control == "sample") %>%
  select(-control) #get rid of this column


# sometimes we'll use more than a single replicate, this broke a previous version of this code chunk, for now we'll find the triplicate (or more) samples, and keep the two with greatest on target read depth (ignoring triplicates quadruplicates etc)

genos_0.11 %<>% 
  group_by(sample_simple) %>%
  slice_max(order_by = `On-Target Reads`, n = 2)


#this writes a new dataset (0.2) by choosing the samples within duplicates and keeping the one with the highest genotyping success
genos_0.2 <- genos_0.11 %>%
  group_by(sample_simple) %>%
  filter(`On-Target Reads` == max(`On-Target Reads`))



##load marker mapping to categorize number of neutral, adaptive markers##

marker_mapping <- readxl::read_xlsx("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Omykiss_Siletz_2023/Omy_GTseq_panel/final_mapping_results.xlsx", sheet = 1)

initial_marker_summary <- marker_mapping %>%
  mutate(marker = str_replace(marker, "Omy(\\d+)", "Chr\\1")) %>% # convert names to match
  filter(marker %in% colnames(genos_0.1)) %>% # exclude filtered loci
  mutate(neutral = if_else(str_detect(`Presumed Type`, 'neutral|Neutral'), "neutral", "adaptive")) %>%
  mutate(run_timing = if_else(str_detect(`Presumed Type`, 'run|Run'), "run_timing", "non-run_timing")) %>%
  dplyr::select(marker, chr, start, 'Presumed Type', neutral, run_timing, Source)


initial_neutral_markers<-subset(initial_marker_summary, initial_marker_summary$neutral == "neutral")

##251 neutral markers in initial marker set##

initial_adaptive_markers<-subset(initial_marker_summary, initial_marker_summary$neutral == "adaptive")

##132 adaptive markers## ##Two species ID markers and 1 sex marker in this list## 
##6 omy-5 pearse marker not included in this list but in marker panel##

initial_run_timing_markers<-subset(initial_marker_summary, initial_marker_summary$run_timing == "run_timing")

##14 Omy-28 run timing markers in initial panel##



##Load filtered dataset

load("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Omykiss_Siletz_2023/filtering_record_filtered_genotypes/Omykiss_Siletz_2022_filtered_genotypes_2.2.R")

##check that sample_simple orders match##

identical(genos_0.2[['sample_simple']], genos_2.0[['sample_simple']])


##Two markers of interest (OmyR40319) and (Chr28_11632591) were removed during filtering becasue they are monomorphic##
##add them back into filtered dataset##


genos_2.0<-as_tibble(cbind(genos_2.0, genos_0.2$OmyR40319Pearse, genos_0.2$Chr28_11632591))


##remove individual that is heterozygous at all 4 species ID markers## 


genos_2.0<-genos_2.0[-125, ]


##remove species ID markers##

species_ID<-c('Omy_Omyclmk438-96', 'Omy_myclarp404-111', 'Oki120255-113')

genos_2.0 <- genos_2.0 %>%
  select(-one_of(species_ID))

##for whatever reason need to remove the species ID marker (Ocl_gshpx-357) manually##

genos_2.0<-genos_2.0[, -355]


##load marker mapping to categorize number of neutral, adaptive markers##
marker_summary <- marker_mapping %>%
  mutate(marker = str_replace(marker, "Omy(\\d+)", "Chr\\1")) %>% # convert names to match
  filter(marker %in% colnames(genos_2.0)) %>% # exclude filtered loci
  mutate(neutral = if_else(str_detect(`Presumed Type`, 'neutral|Neutral'), "neutral", "adaptive")) %>%
  mutate(run_timing = if_else(str_detect(`Presumed Type`, 'run|Run'), "run_timing", "non-run_timing")) %>%
  dplyr::select(marker, chr, start, 'Presumed Type', neutral, run_timing, Source)

neutral_markers<-subset(marker_summary, marker_summary$neutral == "neutral")
adaptive_markers<-subset(marker_summary, marker_summary$neutral == "adaptive")
run_timing_markers<-subset(marker_summary, marker_summary$run_timing == "run_timing")

adaptive_nonruntiming_markers<-subset(marker_summary, marker_summary$neutral == "adaptive" & marker_summary$run_timing == "non-run_timing")



##subset omy5 markers, quantify genotypic frequencies, and plot##


markers<-c("sample_simple", "OmyY1_2SEXY", "OmyR24370Pearse", "OmyR40252Pearse", "OmyR19198Pearse", "OmyR33562Pearse", "OmyR14589Pearse", "genos_0.2$OmyR40319Pearse")

genos_omy5<-genos_2.0[, markers]


##OmyR24370Pearse##

locus_1<-genos_omy5 %>%
  summarize(n = n(), AA = sum(`OmyR24370Pearse` == "AA") / n, AR = sum(`OmyR24370Pearse` == "GA") / n, RR = sum(`OmyR24370Pearse` == "GG") / n, Missing = sum(`OmyR24370Pearse` == "00") / n)

locus_1 %<>%
  pivot_longer(cols = c(AA, AR, RR, Missing))

locus_1$locus<-rep("OmyR24370", times = nrow(locus_1))


##OmyR40252Pearse##

locus_2<-genos_omy5 %>% 
  summarize(n = n(), AA = sum(genos_omy5$OmyR40252Pearse == "AA") / n, AR = sum(genos_omy5$OmyR40252Pearse == "TA") / n, RR = sum(genos_omy5$OmyR40252Pearse == "TT") / n, Missing = sum(genos_omy5$OmyR40252Pearse == "00") / n)

locus_2 %<>%
  pivot_longer(cols = c(AA, AR, RR, Missing))

locus_2$locus<-rep("OmyR40252", times = nrow(locus_2))


##OmyR19198Pearse##


locus_3<-genos_omy5 %>%
  summarize(n = n(), AA = sum(`OmyR19198Pearse` == "AA") / n, AR = sum(`OmyR19198Pearse` == "GA") / n, RR = sum(`OmyR19198Pearse` == "GG") / n, Missing = sum(`OmyR19198Pearse` == "00") / n)

locus_3 %<>%
  pivot_longer(cols = c(AA, AR, RR, Missing))

locus_3$locus<-rep("OmyR19198", times = nrow(locus_3))


##OmyR33562Pearse##

locus_4<-genos_omy5 %>%
  summarize(n = n(), AA = sum(`OmyR33562Pearse` == "GG") / n, AR = sum(`OmyR33562Pearse` == "AG") / n, RR = sum(`OmyR33562Pearse` == "AA") / n, Missing = sum(`OmyR33562Pearse` == "00") / n)

locus_4 %<>%
  pivot_longer(cols = c(AA, AR, RR, Missing))

locus_4$locus<-rep("OmyR33562", times = nrow(locus_4))


##OmyR14589Pearse##

locus_5<-genos_omy5 %>%
  summarize(n = n(), AA = sum(`OmyR14589Pearse` == "AA") / n, AR = sum(`OmyR14589Pearse` == "TA") / n, RR = sum(`OmyR14589Pearse` == "TT") / n, Missing = sum(`OmyR14589Pearse` == "00") / n)

locus_5 %<>%
  pivot_longer(cols = c(AA, AR, RR, Missing))

locus_5$locus<-rep("OmyR14589", times = nrow(locus_5))


#OmyR40319Pearse##

locus_6<-genos_omy5 %>%
  summarize(n = n(), AA = sum(`genos_0.2$OmyR40319Pearse` == "TT") / n, AR = sum(`genos_0.2$OmyR40319Pearse` == "TC") / n, RR = sum(`genos_0.2$OmyR40319Pearse` == "CC") / n, Missing = sum(`genos_0.2$OmyR40319Pearse` == "00") / n)

locus_6 %<>%
  pivot_longer(cols = c(AA, AR, RR, Missing))

locus_6$locus<-rep("OmyR40319", times = nrow(locus_6))



all_loci<-rbind(locus_1, locus_2, locus_3, locus_4, locus_5, locus_6)

all_loci$locus<-as.character(all_loci$locus)
all_loci$locus<-factor(all_loci$locus, levels = c("OmyR24370", "OmyR40252", "OmyR19198", "OmyR33562", "OmyR14589", "OmyR40319"))

all_loci$name<-as.character(all_loci$name)
all_loci$name<-factor(all_loci$name, levels = c("AA", "AR", "RR", "Missing"))



ggplot(all_loci, aes(fill = name, y = value, x = locus)) +
  geom_bar(position = "stack", stat = "identity", width = 0.95, colour = 'black', size = 1.2) +
  scale_fill_manual(values = c("gold2", "coral3", "steelblue4", "gray70"), name = "Genotype") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") + 
  ylab("") +
  theme(panel.background = element_blank(),
        strip.background = element_rect(color = "black", linewidth = 2, fill = "white"),
        text = element_text(family = "sans"),
        strip.text.x = element_text(size = 16, face = "bold"),
        plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        axis.line = element_line(size = 2, color = 'black'),
        axis.ticks = element_line(size = 2, color = 'black'), 
        axis.ticks.length = unit(0.5, "cm"),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust= 1, size = 20, colour = 'black'), 
        axis.text.y = element_text(size = 20, colour = 'black')
  )


##subset omy28 markers, quantify genotypic frequencies, and plot##

omy_28_markers<-c("sample_simple", "OmyY1_2SEXY", "Chr28_11625241", "Chr28_11658853", "Chr28_11667578", "Omy_RAD15709-53", "Chr28_11676622", "Chr28_11683204", "Chr28_11702210")

genos_omy_28<-genos_2.0[, omy_28_markers]




##Chr28_11625241##

locus_1<-genos_omy_28 %>%
  summarize(n = n(), EE = sum(genos_omy_28$Chr28_11625241 == "AA") / n, EL = sum(genos_omy_28$Chr28_11625241 == "AG") / n, LL = sum(genos_omy_28$Chr28_11625241 == "GG") / n, Missing = sum(genos_omy_28$Chr28_11625241 == "00") / n)

locus_1 %<>%
  pivot_longer(cols = c(EE, EL, LL, Missing))

locus_1$locus<-rep("Omy28_11625241", times = nrow(locus_1))


##Chr28_11658853##

locus_2<-genos_omy_28 %>% 
  summarize(n = n(), EE = sum(genos_omy_28$Chr28_11658853 == "AA") / n, EL = sum(genos_omy_28$Chr28_11658853 == "AC") / n, LL = sum(genos_omy_28$Chr28_11658853 == "CC") / n, Missing = sum(genos_omy_28$Chr28_11658853 == "00") / n)

locus_2 %<>%
  pivot_longer(cols = c(EE, EL, LL, Missing))

locus_2$locus<-rep("Omy28_11658853", times = nrow(locus_2))


##Chr28_11667578##


locus_3<-genos_omy_28 %>%
  summarize(n = n(), EE = sum(genos_omy_28$Chr28_11667578 == "TT") / n, EL = sum(genos_omy_28$Chr28_11667578 == "TC") / n, LL = sum(genos_omy_28$Chr28_11667578 == "CC") / n, Missing = sum(genos_omy_28$Chr28_11667578 == "00") / n)

locus_3 %<>%
  pivot_longer(cols = c(EE, EL, LL, Missing))

locus_3$locus<-rep("Omy28_1167578", times = nrow(locus_3))


##Omy_RAD15709-53##

locus_4<-genos_omy_28 %>%
  summarize(n = n(), EE = sum(genos_omy_28$`Omy_RAD15709-53` == "AA") / n, EL = sum(genos_omy_28$`Omy_RAD15709-53` == "GA") / n, LL = sum(genos_omy_28$`Omy_RAD15709-53` == "GG") / n, Missing = sum(genos_omy_28$`Omy_RAD15709-53` == "00") / n)

locus_4 %<>%
  pivot_longer(cols = c(EE, EL, LL, Missing))

locus_4$locus<-rep("Omy_RAD15709-53", times = nrow(locus_4))


##Chr28_11676622##

locus_5<-genos_omy_28 %>%
  summarize(n = n(), EE = sum(genos_omy_28$Chr28_11676622 == "TT") / n, EL = sum(genos_omy_28$Chr28_11676622 == "TG") / n, LL = sum(genos_omy_28$Chr28_11676622 == "GG") / n, Missing = sum(genos_omy_28$Chr28_11676622 == "00") / n)

locus_5 %<>%
  pivot_longer(cols = c(EE, EL, LL, Missing))

locus_5$locus<-rep("Omy28_11676622", times = nrow(locus_5))


#Chr28_11683204##

locus_6<-genos_omy_28 %>%
  summarize(n = n(), EE = sum(genos_omy_28$Chr28_11683204 == "GG") / n, EL = sum(genos_omy_28$Chr28_11683204 == "GT") / n, LL = sum(genos_omy_28$Chr28_11683204 == "TT") / n, Missing = sum(genos_omy_28$Chr28_11683204 == "00") / n)

locus_6 %<>%
  pivot_longer(cols = c(EE, EL, LL, Missing))

locus_6$locus<-rep("Omy28_11683204", times = nrow(locus_6))


#Chr28_11702210##

locus_7<-genos_omy_28 %>%
  summarize(n = n(), EE = sum(genos_omy_28$Chr28_11702210 == "TT") / n, EL = sum(genos_omy_28$Chr28_11702210 == "TG") / n, LL = sum(genos_omy_28$Chr28_11702210 == "GG") / n, Missing = sum(genos_omy_28$Chr28_11702210 == "00") / n)

locus_7 %<>%
  pivot_longer(cols = c(EE, EL, LL, Missing))

locus_7$locus<-rep("Omy28_11702210", times = nrow(locus_7))


all_loci<-rbind(locus_1, locus_2, locus_3, locus_4, locus_5, locus_6, locus_7)

all_loci$locus<-as.character(all_loci$locus)
all_loci$locus<-factor(all_loci$locus, levels = c("Omy28_11625241", "Omy28_11658853", "Omy28_1167578", "Omy_RAD15709-53", "Omy28_11676622", "Omy28_11683204", "Omy28_11702210"))

all_loci$name<-as.character(all_loci$name)
all_loci$name<-factor(all_loci$name, levels = c("EE", "EL", "LL", "Missing"))


ggplot(all_loci, aes(fill = name, y = value, x = locus)) +
  geom_bar(position = "stack", stat = "identity", width = 0.95, colour = 'black', size = 1.2) +
  scale_fill_manual(values = c("gold2", "coral3", "steelblue4", "gray70"), name = "Genotype") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") + 
  ylab("") +
  theme(panel.background = element_blank(),
        strip.background = element_rect(color = "black", linewidth = 2, fill = "white"),
        text = element_text(family = "sans"),
        strip.text.x = element_text(size = 16, face = "bold"),
        plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        axis.line = element_line(size = 2, color = 'black'),
        axis.ticks = element_line(size = 2, color = 'black'), 
        axis.ticks.length = unit(0.5, "cm"),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust= 1, size = 20, colour = 'black'), 
        axis.text.y = element_text(size = 20, colour = 'black')
  )
