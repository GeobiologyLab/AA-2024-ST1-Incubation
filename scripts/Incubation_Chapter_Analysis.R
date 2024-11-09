# Note: ChatGPT (version 4) was used to help clean up and comment parts of the code.

library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(Polychrome)
library(dplyr)
library(ggforce)
library(ggrepel)
library(plotly)
library(shiny)
library(ggpubr)
library(ggpattern)
library(tidyr)
library(ggtext)
library(cowplot)
library(pals)
library(viridis)
library(forcats)
library(grid)
library(ggh4x)

########### Script for Incubation Thesis Chapter Figures ###########

# The dada0 object was created with pseudopool = FALSE.


# Initial phyloseq object creation

# dada0 and seqtab.nochim
da <- load("../data/dada2_May2.Rdata")
# annotations
da2 <- load("../data/dada2_May2_annotations.Rdata")

theme_set(theme_bw())

#create sample data for phyloseq object
samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(matrix(NA, nrow=75, ncol=2))
colnames(samdf)[1] <- "ID"
colnames(samdf)[2] <- "Sample_or_Control"

samplenames <- read.csv("../data/samplenames.csv", header=F)
samplenames2 <- read.csv("../data/samplenames2.csv", header=F)
expnames <- read.csv("../data/expnames.csv", header=F)
timepoints <- read.csv("../data/timepoints.csv", header=T)

data_string <- as.character(samplenames$V1)
rownames(seqtab.nochim) <- data_string
rownames(samdf) <- data_string
samdf$ID <- data_string
samdf$Sample_or_Control <- "sample"
samdf$Sample_or_Control[c(12,24,25,37,48,50,62,75)] <- "negative"
samdf <- cbind(samdf, timepoints)
sampnames <- as.character(samplenames2$V1)
samdf$Name <- sampnames
expnames <- as.character(expnames$V1)
samdf$Type <- expnames

#create phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa), sample_data(samdf))

#contamination 

taxonomydf = data.frame(taxa)
df = data.frame(seqtab.nochim)

df = df/rowSums(df)

cntrl_IDs = rownames(seqtab.nochim)[grepl("Blank", rownames(seqtab.nochim))]
controlsdf = seqtab.nochim[which(rownames(seqtab.nochim)%in%cntrl_IDs),]
controlsdf = controlsdf[,-which(colSums(controlsdf)==0)]

contaminants = colnames(controlsdf)


desulforudis = data.frame(rownames(taxonomydf[which(taxonomydf$Genus=="Candidatus Desulforudis"),]))
desulforudis$Location = "Sample"
colnames(desulforudis) = c("Sequence","Location")
desulforudis[which(desulforudis$Sequence%in%contaminants),]$Location = "Contaminant"

hyd = data.frame(rownames(taxonomydf[which(taxonomydf$Genus=="Hydrogenophaga"),]))
hyd$Location = "Sample"
colnames(hyd) = c("Sequence","Location")
hyd[which(hyd$Sequence%in%contaminants),]$Location = "Contaminant"


therm = data.frame(rownames(taxonomydf[which(taxonomydf$Class=="Thermodesulfovibrionia"),]))
therm$Location = "Sample"
colnames(therm) = c("Sequence","Location")

# removing ASVs that are contaminants (except Desulforudis)
cntrl_IDs = rownames(seqtab.nochim)[grepl("Blank", rownames(seqtab.nochim))]
controlsdf = seqtab.nochim[which(rownames(seqtab.nochim)%in%cntrl_IDs),]
controlsdf = controlsdf[,-which(colSums(controlsdf)==0)]
controlsdf = as.data.frame(controlsdf)
# don't include the two "contaminant" desulforudis ASVs
controlsdf = controlsdf[,-which(colnames(controlsdf)=="TACGTAGGGGGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGCTTTAGAAGTCGCAGGTGAAATCCCGCGGCTTAACCGTGGAACTGCCTGCGAAACCATTTAGCTTGAGGACAGGAGAGGGAAGCGGAATTCCTGGTGTAGCGGTGAAATGCGTAGATATCAGGAGGAACACCAGTGGCGAAGGCGGCTTTCTGGCCTGGTCCTGACGCTGAGGCGCGAAAGCTAGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCTAGCTGTAAACGATGGGCACTAGGTGTTGGGGGGTTCATACCCTTCAGTGCCGTAGGTAACCCAATAAGTGCCCCGCCTGGGGAGTACGGTCGCAAGACTG")]
controlsdf = controlsdf[,-which(colnames(controlsdf)=="TACGTAGGGGGCGAGCGTTGTCCGGAATCACTGGGCGTAGAGGGCGTGTAGGCGGCTAGAGTAGTCGCAGGTGAAATCCCACGGCTCAACCGTGGAACTGCCTGCGAAACCATCTAGCTTGAGGACAGGAGAGGAAAGCGGAATTCCTGGTGTAGCGGTGAAATGCGTAGATATCAGGAGGAACATCAGTGGCGAAGGCGGCTTTCTGGCCTGGTCCTGACGCTGAGGCGCGAAAGCTAGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCTAGCTGTAAACGATGGGTACTAGGTGTTGGTGGGCTCATACCCGTCAGTGCCGCAGTTAACACAATAAGTACCCCGCCTGGGGAGTACGGTCGCAAGACTG")]

contaminants = colnames(controlsdf)
df <- seqtab.nochim
common_cols <- intersect(colnames(df), colnames(controlsdf))
df <- df[ , !(colnames(df) %in% common_cols)]

# Filter samples based on metadata
sample_names_to_keep <- sample_data(ps)$ID[sample_data(ps)$Sample_or_Control == "sample"]
ps_filtered <- prune_samples(sample_names_to_keep, ps)

df <- otu_table(df, taxa_are_rows=F)
otu_table(ps_filtered) <- df


# rarefy samples
# Extracting the OTU table properly
otu_mat <- as.matrix(otu_table(ps_filtered))

# Check for any NA values
sum(is.na(otu_mat))

# Ensure all entries are numeric
all(is.numeric(otu_mat))


# Rarefaction function
simple_rarefy <- function(otu_mat, depth, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  apply(otu_mat, 1, function(row) {
    total = sum(row)
    if (total >= depth) {
      probs = row / total
      rmultinom(n = 1, size = depth, prob = probs)
    } else {
      warning("Sample size less than depth, returning original counts")
      row
    }
  })
}

set.seed(117) 

# Perform rarefaction with a specified seed
rarefied_otu <- simple_rarefy(otu_mat, 77000, seed = 117)

# Check the output after manual rarefaction
colSums(rarefied_otu)

# Update the phyloseq object
ps_rare <- ps_filtered
rownames(rarefied_otu) <- colnames(otu_table(ps_filtered))
ps_rare@otu_table <- otu_table(as(rarefied_otu, "matrix"), taxa_are_rows = T, errorIfNULL = TRUE)

# Identify the samples to keep (those that do not contain "1L" in their names)
samples_to_keep <- !grepl("1L", sample_names(ps_rare))

# Prune the phyloseq object to keep only the selected samples
ps_rare_no1L <- prune_samples(samples_to_keep, ps_rare)

# Check the sample names in the pruned phyloseq object to ensure "1L" samples are removed
sample_names(ps_rare_no1L)




############ FIGURE 3.1 - Microbial Dynamics during Incubation ###########

############A###########
dataTrip <- read.csv("../data/IncubationResultsTrips.csv", sep = ";")

# Replace negative values with the small positive value
dataTrip_mod <- dataTrip %>% mutate(H2O2 = ifelse(H2O2 < 0, 0.00625, H2O2))
dataTrip_mod2 <- dataTrip_mod[-c(64:67),]
dataTrip_mod2$Time[2] <- 0
dataTrip_mod2$Time[3] <- 0

cbPalette2 <- c("#FFCC00","#CC0000", "#990099", "#009E73", "#56B4E9")

color_palette <- c(  
  rep("#FFCC00", 3),
  rep("#56B4E9", 15),
  rep("#009E73", 15),
  rep("#CC0000", 15),
  rep("#990099", 15)
)

desired_order <- c(
  "t0.1","t0.2","t0.3", 
  "W1.1", "W1.2", "W1.3", "W2.1", "W2.2",
  "W2.3", "W3.1", "W3.2", "W3.3", "W4.1",
  "W4.2", "W4.3", "W5.1", "W5.2", "W5.3",
  "G1.1", "G1.2", "G1.3", "G2.1", "G2.2",
  "G2.3", "G3.1", "G3.2", "G3.3", "G4.1",
  "G4.2", "G4.3", "G5.1", "G5.2", "G5.3",
  "TAn1.1", "TAn1.2", "TAn1.3", "TAn2.1", "TAn2.2",
  "TAn2.3", "TAn3.1", "TAn3.2", "TAn3.3", "TAn4.1",
  "TAn4.2", "TAn4.3", "TAn5.1", "TAn5.2", "TAn5.3",
  "TO1.1", "TO1.2", "TO1.3", "TO2.1", "TO2.2",
  "TO2.3", "TO3.1", "TO3.2", "TO3.3", "TO4.1",
  "TO4.2", "TO4.3", "TO5.1", "TO5.2", "TO5.3")


dataTrip_mod2$Sample <- factor(dataTrip_mod2$Sample, levels = desired_order)
dataTrip_mod2 <- dataTrip_mod2 %>%
  mutate(width = ifelse(Time %in% c(0), 0.3, 0.8))

# I removed the legend because I cant get it to work. I can edit this to make fill=Experiment to get the legend I want and append it to my other plot
ggplot(dataTrip_mod2, aes(x = as.factor(Time), fill = Sample)) +
  geom_bar(aes(y = H2O2, width=width), stat="identity", position=position_dodge(width = 0.8), width = 0.8, color="black") +
  labs(x = "Time (Hours)", y = "Hydrogen Peroxide (μM)", title = "Hydrogen Peroxide - Incubation V2") +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    legend.margin = margin(t = 10, r = 10, b = 10, l = 30),
    axis.title.x = element_text(size = 16, margin = margin(t = 20)),  
    axis.title.y = element_text(size = 16, margin = margin(r = 20)),  
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  scale_fill_manual(values = color_palette) +
  #scale_color_manual(name='Experiment',
  #breaks=c('Glass', 'Rock Anoxic', 'Rock Oxic', 'Water'),
  #values=c('Glass'='#009E73', 'Rock Anoxic'='#CC0000', 'Rock Oxic'='#990099', 'Water'='#56B4E9'))+
  scale_y_continuous(labels = function(x) format(x, scientific = F)) 


###########B###########
# Define the desired order of the types
desired_order <- c("Initial", "Water", "Glass", "RockAnoxic", "RockOxic")

# Convert the 'Type' column to a factor with the desired order
dataTrip_mod2$Experiment <- factor(dataTrip_mod2$Experiment, levels = desired_order)

color_palette3 <- c(
  rep("#FFCC00", 3),
  rep("#56B4E9", 15),
  rep("#009E73", 15),
  rep("#CC0000", 15),
  rep("#990099", 15)
)

cbPalette4 <- c("#FFCC00","#009E73", "#CC0000", "#990099", "#56B4E9")

dataTrip_mod2 <- dataTrip_mod2 %>%
  mutate(width = ifelse(Time %in% c(0), 0.3, 1))

g <- ggplot(dataTrip_mod2, aes(x = as.factor(Time), fill = Sample)) +
  geom_bar(aes(y = log10(qPCR_Count), width=width), stat="identity", position=position_dodge(width = 0.9),color="black") +
  #geom_errorbar(aes(ymin = qPCR_Count - qPCR_StdDev, ymax = qPCR_Count + qPCR_StdDev), position = position_dodge(width = 0.8), width = 0.5) +
  labs(x = "Time (Hours)", y = expression(log("Copy Number mL"^"-1"*"")), title = "qPCR - Incubation V2") +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),legend.margin = margin(t = 10, r = 10, b = 10, l = 30),
        axis.title.x = element_text(size = 16, margin = margin(t = 20)),  
        axis.title.y = element_text(size = 16, margin = margin(r = 20)),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_fill_manual(values = color_palette3) +
  scale_y_continuous(labels = function(x) format(x, scientific = F)) 

###########C###########

# Convert phyloseq object to matrix format
otu_matrix <- as(otu_table(ps_rare), "matrix")
otu_matrix <- t(otu_matrix)


desired_order <- c(
  "t0_1","t0_2","t0_3", "G1_1", "G1_2", "G1_3", "G2_1", "G2_2",
  "G2_3", "G3_1", "G3_2", "G3_3", "G4_1",
  "G4_2", "G4_3", "G5_1", "G5_2", "G5_3",
  "G1L", "W1_1", "W1_2", "W1_3", "W2_1", "W2_2",
  "W2_3", "W3_1", "W3_2", "W3_3", "W4_1",
  "W4_2", "W4_3", "W5_1", "W5_2", "W5_3",
  "W1L", "TO1_1", "TO1_2", "TO1_3", "TO2_1", "TO2_2",
  "TO2_3", "TO3_1", "TO3_2", "TO3_3", "TO4_1",
  "TO4_2", "TO4_3", "TO5_1", "TO5_2", "TO5_3",
  "TO1L", "TAn1_1", "TAn1_2", "TAn1_3", "TAn2_1", "TAn2_2",
  "TAn2_3", "TAn3_1", "TAn3_2", "TAn3_3", "TAn4_1",
  "TAn4_2", "TAn4_3", "TAn5_1", "TAn5_2", "TAn5_3",
  "TAn1L")

# Apply the desired order to the 'Sample' column in both datasets
distance_long$Sample1 <- factor(distance_long$Sample1, levels = desired_order)

# Function to determine Type based on Sample name prefixes
getType <- function(sample_name) {
  if (grepl("^G", sample_name)) {
    "Glass"
  } else if (grepl("^W", sample_name)) {
    "Water"
  } else if (grepl("^TO", sample_name)) {
    "Rock Oxic"
  } else if (grepl("^TAn", sample_name)) {
    "Rock Anoxic"
  } else {
    "Other"  
  }
}

# Extracting OTU table and Taxonomy table from the phyloseq object
otu_data <- as.data.frame(otu_table(ps_rare))
tax_data <- as.data.frame(tax_table(ps_rare))

# Merge OTU data with taxonomic information
merged_data <- cbind(otu_data, tax_data)
merged_data <- merged_data[,-68]
merged_data <- pivot_longer(merged_data, cols = -c(Phylum, Class, Order, Family, Genus), names_to = "Sample", values_to = "Abundance")

# Convert Abundance to numeric
merged_data$Abundance <- as.numeric(as.character(merged_data$Abundance))

# For Genera

genera_data <- merged_data %>%
  group_by(Genus, Sample) %>%
  summarize(Abundance = sum(Abundance), .groups = 'drop') %>%
  mutate(Type = sapply(Sample, getType))  

# Identifying taxa of interest and top taxa
taxa_of_interest <- c("Candidatus Desulforudis", "Hydrogenophaga", "Thermodesulfovibrionia")
top_genera <- genera_data %>%
  filter(!Genus %in% taxa_of_interest) %>%
  group_by(Genus) %>%
  summarize(TotalAbundance = sum(Abundance), .groups = 'drop') %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:20) %>%
  pull(Genus)


# Update the Taxonomy assignment step to handle potential NA or unexpected cases robustly
genera_data <- genera_data %>%
  mutate(Genus = ifelse(is.na(Genus), "Unknown", Genus),  
         Taxonomy = case_when(
           Genus %in% taxa_of_interest ~ Genus,
           Genus %in% top_genera ~ Genus,
           TRUE ~ "Other"  
         ))

# calculate relative abundances 
genera_data <- genera_data %>%
  group_by(Sample) %>%
  mutate(TotalAbundance = sum(Abundance, na.rm = TRUE),
         RelativeAbundance = Abundance / TotalAbundance * 100)

# Include top taxa and taxa of interest for plotting
genera_data$Taxonomy <- ifelse(genera_data$Genus %in% taxa_of_interest, genera_data$Genus, ifelse(genera_data$Genus %in% top_genera, genera_data$Genus, "Other"))

# Calculate relative abundances
calculate_relative_abundance <- function(data) {
  data %>%
    group_by(Sample) %>%
    mutate(TotalAbundance = sum(Abundance)) %>%
    ungroup() %>%
    mutate(RelativeAbundance = Abundance / TotalAbundance * 100)
}

genera_data <- calculate_relative_abundance(genera_data)

# Plot 

c25genera <- c(
  "dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00","red4", "gold1","skyblue2",
  "#FB9A99", "palegreen2","#CAB2D6", "blue1", "gray70", "khaki2","maroon", "orchid1",
  "deeppink1", "#FDBF6F", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")


# Apply the desired order to the 'Sample' column in both datasets
genera_data$Sample <- factor(genera_data$Sample, levels = desired_order)

# Remove rows where the "Sample" column contains "1L"
genera_data_no1L <- genera_data %>% filter(!grepl("1L", Sample))

genera_data_no1L <- genera_data_no1L %>% mutate(Type = ifelse(Type == "Other", "Initial", Type))

type_order <- c("Initial", "Water", "Glass", "Rock Anoxic", "Rock Oxic")

genera_data_no1L$Type<- factor(genera_data_no1L$Type, levels = type_order)

sulfuritalea <- genera_data_no1L %>%
  filter(Genus == "Sulfuritalea")

thiobacillus <- genera_data_no1L %>%
  filter(Genus == "Thiobacillus")


# Genera
ggplot(genera_data_no1L_1, aes(x = Sample, y = RelativeAbundance, fill = Taxonomy)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c25genera) +
  facet_wrap(~Type, scales="free_x", nrow=1) +
  labs(title = "Relative Abundance of Genera",
       x = "Sample",
       y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


########### FIGURE 3.2 - ASVs shared with initial community ###########

# Function to get ASVs for a given sample or subset of samples
get_asvs <- function(ps, sample_names) {
  # Ensure sample_names is a character vector
  sample_names <- as.character(sample_names)
  # Prune samples
  pruned_ps <- prune_samples(sample_names(ps_rare) %in% sample_names, ps_rare)
  
  # Extract OTU table
  otu_table_pruned <- otu_table(pruned_ps)
  
  # Ensure OTU table is in matrix format
  if (!is.matrix(otu_table_pruned)) {
    otu_mat <- as(otu_table_pruned, "matrix")
  } else {
    otu_mat <- otu_table_pruned
  }
  
  # Filter rows where the sum of counts is greater than 0
  non_zero_otus <- otu_mat[rowSums(otu_mat) > 0, ]
  
  # Get the ASVs (row names) where the sum of the row is greater than zero
  asvs <- rownames(non_zero_otus)
  
  # Print the number of ASVs found
  cat("Number of ASVs found:", length(asvs), "\n")
  
  return(asvs)
}

# Define timepoints and treatments
initial_timepoints <- c("t0_1", "t0_2", "t0_3")
later_timepoints <- list(
  rockoxic = c("TO1_1", "TO1_2", "TO1_3"),
  rockanoxic = c("TAn1_1", "TAn1_2", "TAn1_3"),
  water = c("W1_1", "W1_2", "W1_3"),
  glass = c("G1_1", "G1_2", "G1_3")
)

# Extract ASVs for the initial timepoints
initial_asvs_list <- lapply(initial_timepoints, function(sample) get_asvs(ps_rare, sample))
initial_asvs <- Reduce(intersect, initial_asvs_list)

# Extract ASVs for each later timepoint
later_asvs_list <- lapply(later_timepoints, function(tp) lapply(tp, function(sample) get_asvs(ps_rare, sample)))

# Function to count shared ASVs between initial and later timepoints
count_shared_asvs <- function(initial_asvs, later_asvs) {
  length(intersect(initial_asvs, later_asvs))
}

# Calculate the number of shared ASVs for each treatment
shared_asvs_counts <- lapply(later_asvs_list, function(later_asvs) {
  sapply(later_asvs, function(rep_asvs) {
    count_shared_asvs(initial_asvs, rep_asvs)
  })
})

# Print the shared ASVs counts to check
print(shared_asvs_counts)

# Combine the results into a data frame
shared_asvs_df <- do.call(rbind, lapply(names(shared_asvs_counts), function(treatment) {
  data.frame(
    Treatment = treatment,
    Shared_ASVs = shared_asvs_counts[[treatment]]
  )
}))

# Print the data frame to check
print(shared_asvs_df)


# Rename the levels of the Treatment column
shared_asvs_df$Treatment <- factor(shared_asvs_df$Treatment,
                                   levels = c("rockoxic", "rockanoxic", "glass", "water"),
                                   labels = c("Rock Oxic", "Rock Anoxic", "Glass", "Water"))

# Rename the levels of the Treatment column
shared_asvs_df$Treatment <- factor(shared_asvs_df$Treatment,
                                   levels = c("water", "glass", "rockanoxic", "rockoxic"),
                                   labels = c("Water", "Glass", "Rock Anoxic", "Rock Oxic"))

# Print the data frame to check the changes
print(shared_asvs_df)


# Define custom colors for the treatments
custom_colors <- c("Rock Oxic" = "#990099", "Rock Anoxic" = "#CC0000", "Water" = "#56B4E9", "Glass" = "#009E73")

# Create the box plot with statistical comparisons
ggplot(shared_asvs_df, aes(x = Treatment, y = Shared_ASVs, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Shared ASVs Between 0 and 24 hours",
       x = "Treatment",
       y = "Number of Shared ASVs") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 16, margin = margin(r = 15)),
        axis.title.x = element_text(size = 16, margin = margin(t = 15)),  
        axis.text.y = element_text(size = 14)) +  
  geom_hline(yintercept = 272, linetype = "dashed", color = "#FFCC00", size=1) +
  stat_compare_means(comparisons = list(
    c("Rock Oxic", "Water"),
    c("Rock Oxic", "Glass"),
    c("Rock Anoxic", "Water"),
    c("Rock Anoxic", "Glass")
  ), label = "p.signif", method = "t.test")

########### FIGURE 3.3 - Relative abundance heatmap ###########

# Calculate average abundance of each ASV in the initial samples
initial_samples <- subset_samples(ps_rare, sample_data(ps_rare)$Type == "Initial")
initial_abundance <- apply(otu_table(initial_samples), 1, mean)

# Function to calculate percent change
percent_change <- function(abundance, initial_abundance) {
  return((abundance - initial_abundance) / initial_abundance * 100)
}

# Initialize a dataframe to store the percent change for each time point and treatment
percent_change_df <- data.frame(ASV = rownames(otu_table(ps_rare)))

# Treatments and time points
treatments <- c("RockOxic", "RockAnoxic", "Glass", "Water")
timepoints <- 2:6

# Calculate percent change for each treatment and time point
for (treatment in treatments) {
  for (timepoint in timepoints) {
    sample_ids <- sample_data(ps_rare)$ID[
      sample_data(ps_rare)$Type == treatment &
        sample_data(ps_rare)$Timepoint == timepoint
    ]
    treatment_samples <- prune_samples(sample_ids, ps_rare)
    treatment_abundance <- apply(otu_table(treatment_samples), 1, mean)
    change <- percent_change(treatment_abundance, initial_abundance)
    col_name <- paste0(substr(treatment, 1, 5), timepoint)
    percent_change_df[col_name] <- change
  }
}

# Add ASV IDs and taxonomy to the dataframe
percent_change_df$ASV <- rownames(percent_change_df)
tax_table_df <- as.data.frame(tax_table(ps_rare))
tax_table_df$ASV <- rownames(tax_table_df)
rownames(tax_table_df) <- percent_change_df$ASV
percent_change_df <- left_join(percent_change_df, tax_table_df, by = "ASV")

# Columns for treatments and timepoints 
percent_change_columns <- c("RockO2", "RockO3", "RockO4", "RockO5", "RockO6", 
                            "RockA2", "RockA3", "RockA4", "RockA5", "RockA6", 
                            "Glass2", "Glass3", "Glass4", "Glass5", "Glass6", 
                            "Water2", "Water3", "Water4", "Water5", "Water6")


# Replace Inf and NaN with NA
percent_change_df[percent_change_columns] <- lapply(
  percent_change_df[percent_change_columns], 
  function(x) { x[is.infinite(x) | is.nan(x)] <- NA; return(x) }
)
# Group by taxonomy and calculate the mean percent change for each group
collapsed_df <- percent_change_df %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>%
  summarise(across(all_of(percent_change_columns), mean, na.rm = TRUE))


# Keeping only the genera shown in the stacked bar plot

# List of genera to keep
genera_to_keep <- c("Aquabacterium", "Azonexus", "Bellilinea", "Dechloromonas", 
                    "Dehalobacterium", "Dethiobacter", "Limnohabitans", "Malikia", 
                    "Methylotenera", "Ramlibacter", "Rhodobacter", "Roseateles", 
                    "Roseococcus", "Rubrivivax", "Smithella", "Sulfuritalea", 
                    "Syntrophomonas", "Thiobacillus")

# Subset the dataframe to include only the specified genera
subset_df <- percent_change_df %>% filter(Genus %in% genera_to_keep)


# Columns for treatments and timepoints 
percent_change_columns <- c("RockO2", "RockO3", "RockO4", "RockO5", "RockO6", 
                            "RockA2", "RockA3", "RockA4", "RockA5", "RockA6", 
                            "Glass2", "Glass3", "Glass4", "Glass5", "Glass6", 
                            "Water2", "Water3", "Water4", "Water5", "Water6")



# Group by taxonomy and calculate the mean percent change for each group
collapsed_subset_df <- subset_df %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>%
  summarise(across(all_of(percent_change_columns), mean, na.rm = TRUE))



genera_to_keep2 <- c(genera_to_keep, "Candidatus Desulforudis", "Hydrogenophaga")

# Now the same but with relative abundances

subset_genera_data <- genera_data_no1L %>% filter(Genus %in% genera_to_keep2)
subset_genera_data <- subset_genera_data[,c(-5,-6)]


# Filter initial samples to calculate the average relative abundance
initial_avg <- subset_genera_data %>%
  filter(Type == "Initial") %>%
  group_by(Genus) %>%
  summarise(InitialRelAbundance = mean(RelativeAbundance, na.rm = TRUE))

# Function to calculate percent change
percent_change <- function(current, initial) {
  return((current - initial) / initial * 100)
}

# Calculate percent change for each sample compared to the initial average
percent_change_df <- subset_genera_data %>%
  left_join(initial_avg, by = "Genus") %>%
  mutate(PercentChange = percent_change(RelativeAbundance, InitialRelAbundance)) %>%
  select(Genus, Sample, PercentChange)

# Separate Sample into Treatment and Timepoint
percent_change_df <- percent_change_df %>%
  mutate(Treatment = substr(Sample, 1, 1),
         Timepoint = substr(Sample, 2, nchar(as.character(Sample))))

# Combine Treatment and Timepoint into a single column for wide format
percent_change_df <- percent_change_df %>%
  unite("SampleID", Treatment, Timepoint, sep = "")

percent_change_df <- percent_change_df[,-4]

# Reshape the data to wide format
percent_change_wide <- percent_change_df %>%
  pivot_wider(names_from = Sample, values_from = PercentChange)

# Melt the wide dataframe into a long format suitable for ggplot2
percent_change_long <- percent_change_wide %>%
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "PercentChange")

# Filter out NaN and Inf values
percent_change_long <- percent_change_long %>%
  filter(!is.na(PercentChange) & !is.infinite(PercentChange))

# Log-transform the percent change values (adding a small constant to avoid log(0))
percent_change_long <- percent_change_long %>%
  mutate(LogPercentChange = log10(abs(PercentChange) + 1) * sign(PercentChange))


# Add the Type column based on Sample1 using the getType function
percent_change_long$Type <- sapply(percent_change_long$Sample, getType)

# Filter out rows with Type "Other"
percent_change_long <- percent_change_long %>% filter(Type != "Other")

# Convert the 'Type' column to a factor with the desired order
type_order2 <- c("Rock Anoxic", "Rock Oxic", "Glass", "Water")
percent_change_long$Type <- factor(percent_change_long$Type, levels = type_order2)

# Convert the 'Type' column to a factor with the desired order
heatmap_order <- c("Aquabacterium", "Bellilinea", "Candidatus Desulforudis", "Dechloromonas", 
                   "Dehalobacterium", "Dethiobacter", "Hydrogenophaga", "Limnohabitans", "Malikia", 
                   "Methylotenera", "Ramlibacter", "Rhodobacter", 
                   "Roseococcus", "Rubrivivax", "Smithella", "Sulfuritalea", 
                   "Syntrophomonas", "Thiobacillus")

percent_change_long$Genus <- factor(percent_change_long$Genus, levels = heatmap_order)

# Reverse the levels of the Genus factor
percent_change_long <- percent_change_long %>%
  mutate(Genus = fct_rev(factor(Genus)))


# Plot heatmap with log-transformed values
ggplot(percent_change_long, aes(x = Sample, y = Genus, fill = LogPercentChange)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "royalblue4", mid = "white", high = "firebrick", midpoint = 0, 
                       limits = c(min(percent_change_long$LogPercentChange, na.rm = TRUE), 
                                  max(percent_change_long$LogPercentChange, na.rm = TRUE)),
                       name = expression(log[10]~"(Percent Change)")) +
  theme_minimal() +
  facet_wrap(~Type, scales="free_x", nrow=1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Log-Transformed Percent Change in Relative Abundances compared to initial average", 
       x = "Sample", y = "Genus")

########### FIGURE 3.4 - Growth of selected organisms ###########

max <- read.csv("../data/copyno3.csv", sep = ",")
max$Time <- factor(max$Time, levels = c("0","24","Max"))
max$Treatment <- factor(max$Treatment, levels = c("Water","Glass","Rock Anoxic", "Rock Oxic"))

# first Total

max_tot <- max[,c(1:6)]

# Add columns for the log-transformed values and error bounds
max_tot$log_Total_CopyNo <- log10(max_tot$Total_CopyNo)
max_tot$log_Total_CopyNo_lower <- log10(max_tot$Total_CopyNo - max_tot$Total_StdDev)
max_tot$log_Total_CopyNo_upper <- log10(max_tot$Total_CopyNo + max_tot$Total_StdDev)

# Plotting the log-transformed data
ggplot(max_tot, aes(x = Time, y = log_Total_CopyNo, group = Treatment, color = Treatment)) +
  geom_line(size = 1.2) +  # Add lines for each Type
  geom_point(size = 3) +  # Add points for each Type
  geom_errorbar(aes(ymin = log_Total_CopyNo_lower, ymax = log_Total_CopyNo_upper), 
                width = 0.1, size = 1) +  # Add error bars
  scale_x_discrete(labels = c("Initial", "24 hours", bquote(t[~"max total abundance"]))) +  # Rename x-axis labels
  labs(x = "Time", y = expression(Log[10]("Total Copy Number ")), title = "Total") +
  theme_bw() +
  scale_color_manual(values = c("Water" = "#56B4E9", "Glass" = "#009E73", 
                                "Rock Anoxic" = "#CC0000", "Rock Oxic" = "#990099")) +  # Customize colors
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.title.y = element_text(size = 14, margin = margin(r = 15)),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "none"
  ) +
  ylim(0,5)


# CDO

max_cdo <- max[,c(1:4,7,8)]

# Add columns for the log-transformed values and error bounds
max_cdo$log_Desulfo_CopyNo <- log10(max_cdo$Desulfo_CopyNo)
max_cdo$log_Desulfo_CopyNo_lower <- log10(max_cdo$Desulfo_CopyNo - max_cdo$Desulfo_StdDev)
max_cdo$log_Desulfo_CopyNo_upper <- log10(max_cdo$Desulfo_CopyNo + max_cdo$Desulfo_StdDev)

# Plotting the log-transformed data
ggplot(max_cdo, aes(x = Time, y = log_Desulfo_CopyNo, group = Treatment, color = Treatment)) +
  geom_line(size = 1.2) +  # Add lines for each Type
  geom_point(size = 3) +  # Add points for each Type
  geom_errorbar(aes(ymin = log_Desulfo_CopyNo_lower, ymax = log_Desulfo_CopyNo_upper), 
                width = 0.1, size = 1) +  # Add error bars
  scale_x_discrete(labels = c("Initial", "24 hours", bquote(t[~"max total abundance"]))) +  # Rename x-axis labels
  labs(x = "Time", y = expression(Log[10]("Total Copy Number ")), title = "Candidatus Desulforudis onstottii") +
  theme_bw() +
  scale_color_manual(values = c("Water" = "#56B4E9", "Glass" = "#009E73", 
                                "Rock Anoxic" = "#CC0000", "Rock Oxic" = "#990099")) +  # Customize colors
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.title.y = element_text(size = 14, margin = margin(r = 15)),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "none"
  ) +
  ylim(0,5)

# Hydrogenophaga

max_hyd <- max[,c(1:4,9:10)]

# Add columns for the log-transformed values and error bounds
max_hyd$log_Hydrogeno_CopyNo <- log10(max_hyd$Hydrogeno_CopyNo)
max_hyd$log_Hydrogeno_CopyNo_lower <- log10(max_hyd$Hydrogeno_CopyNo - max_hyd$Hydrogeno_StdDev)
max_hyd$log_Hydrogeno_CopyNo_lower[max_hyd$log_Hydrogeno_CopyNo_lower < 0] <- 0
max_hyd$log_Hydrogeno_CopyNo_upper <- log10(max_hyd$Hydrogeno_CopyNo + max_hyd$Hydrogeno_StdDev)

# Plotting the log-transformed data
ggplot(max_hyd, aes(x = Time, y = log_Hydrogeno_CopyNo, group = Treatment, color = Treatment)) +
  geom_line(size = 1.2) +  # Add lines for each Type
  geom_point(size = 3) +  # Add points for each Type
  geom_errorbar(aes(ymin = log_Hydrogeno_CopyNo_lower, ymax = log_Hydrogeno_CopyNo_upper), 
                width = 0.1, size = 1) +  # Add error bars
  scale_x_discrete(labels = c("Initial", "24 hours", "Max Growth Rate")) +  # Rename x-axis labels
  labs(x = "Time", y = "Log10(Total Copy Number)", title = "Hydrogenophaga") +
  theme_bw() +
  scale_color_manual(values = c("Water" = "#56B4E9", "Glass" = "#009E73", 
                                "Rock Anoxic" = "#CC0000", "Rock Oxic" = "#990099")) +  # Customize colors
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.title.y = element_text(size = 14, margin = margin(r = 15)),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "none"
  )+
  ylim(0,5)

# Unclassified Thermodesulfovibrionia

max_thermo <- max[,c(1:4,11:12)]

# Add columns for the log-transformed values and error bounds
max_thermo$log_Thermo_CopyNo <- log10(max_thermo$Thermo_CopyNo)
max_thermo$log_Thermo_CopyNo_lower <- log10(max_thermo$Thermo_CopyNo - max_thermo$Thermo_StdDev)
max_thermo$log_Thermo_CopyNo_upper <- log10(max_thermo$Thermo_CopyNo + max_thermo$Thermo_StdDev)

# Plotting the log-transformed data
ggplot(max_thermo, aes(x = Time, y = log_Thermo_CopyNo, group = Treatment, color = Treatment)) +
  geom_line(size = 1.2) +  # Add lines for each Type
  geom_point(size = 3) +  # Add points for each Type
  geom_errorbar(aes(ymin = log_Thermo_CopyNo_lower, ymax = log_Thermo_CopyNo_upper), 
                width = 0.1, size = 1) +  # Add error bars
  scale_x_discrete(labels = c("Initial", "24 hours", bquote(t[~"max total abundance"]))) +  # Rename x-axis labels
  labs(x = "Time", y = expression(Log[10]("Total Copy Number ")), title = "Unclassified Thermodesulfovibrionia") +
  theme_bw() +
  scale_color_manual(values = c("Water" = "#56B4E9", "Glass" = "#009E73", 
                                "Rock Anoxic" = "#CC0000", "Rock Oxic" = "#990099")) +  # Customize colors
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.title.y = element_text(size = 14, margin = margin(r = 15)),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "none"
  )+
  ylim(0,5)

# Sulfuritalea

max_sulfu <- max[,c(1:4,13:14)]

# Add columns for the log-transformed values and error bounds
max_sulfu$log_sulfu_CopyNo <- log10(max_sulfu$Sulfuritalea_CopyNo)
max_sulfu$log_sulfu_CopyNo_lower <- log10(max_sulfu$Sulfuritalea_CopyNo - max_sulfu$Sulfuritalea_StdDev)
max_sulfu$log_sulfu_CopyNo_upper <- log10(max_sulfu$Sulfuritalea_CopyNo + max_sulfu$Sulfuritalea_StdDev)

max_sulfu <- max_sulfu %>%
  mutate(log_sulfu_CopyNo_lower = ifelse(log_sulfu_CopyNo_lower < 1, 0, log_sulfu_CopyNo_lower))

# Plotting the log-transformed data
ggplot(max_sulfu, aes(x = Time, y = log_sulfu_CopyNo, group = Treatment, color = Treatment)) +
  geom_line(size = 1.2) +  # Add lines for each Type
  geom_point(size = 3) +  # Add points for each Type
  geom_errorbar(aes(ymin = log_sulfu_CopyNo_lower, ymax = log_sulfu_CopyNo_upper), 
                width = 0.1, size = 1) +  # Add error bars
  scale_x_discrete(labels = c("Initial", "24 hours", bquote(t[~"max total abundance"]))) +  # Rename x-axis labels
  labs(x = "Time", y = expression(Log[10]("Total Copy Number ")), title = "Sulfuritalea") +
  theme_bw() +
  scale_color_manual(values = c("Water" = "#56B4E9", "Glass" = "#009E73", 
                                "Rock Anoxic" = "#CC0000", "Rock Oxic" = "#990099")) +  # Customize colors
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.title.y = element_text(size = 14, margin = margin(r = 15)),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "none"
  )+
  ylim(0,5)

# Thiobacillus

max_thio <- max[,c(1:4,15:16)]

# Add columns for the log-transformed values and error bounds
max_thio$log_thio_CopyNo <- log10(max_thio$Thiobacillus_CopyNo)
max_thio$log_thio_CopyNo_lower <- log10(max_thio$Thiobacillus_CopyNo - max_thio$Thiobacillus_StdDev)
max_thio$log_thio_CopyNo_upper <- log10(max_thio$Thiobacillus_CopyNo + max_thio$Thiobacillus_StdDev)

# Plotting the log-transformed data
ggplot(max_thio, aes(x = Time, y = log_thio_CopyNo, group = Treatment, color = Treatment)) +
  geom_line(size = 1.2) +  # Add lines for each Type
  geom_point(size = 3) +  # Add points for each Type
  geom_errorbar(aes(ymin = log_thio_CopyNo_lower, ymax = log_thio_CopyNo_upper), 
                width = 0.1, size = 1) +  # Add error bars
  scale_x_discrete(labels = c("Initial", "24 hours", bquote(t[~"max total abundance"]))) +  # Rename x-axis labels
  labs(x = "Time", y = expression(Log[10]("Total Copy Number ")), title = "Thiobacillus") +
  theme_bw() +
  scale_color_manual(values = c("Water" = "#56B4E9", "Glass" = "#009E73", 
                                "Rock Anoxic" = "#CC0000", "Rock Oxic" = "#990099")) +  # Customize colors
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.title.y = element_text(size = 14, margin = margin(r = 15)),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "none"
  )+
  ylim(0,5)



########### FIGURE 3.5 - Functionality of MAGs ###########
functions_df <- data.frame(
  organism = c('Ca. D. onstottii', 'Ca. D. sp002294005', 'Ca. D. audaxviator', 'Hydrogenophaga bedretto',
               'Hydrogenophaga sp018333605', 'Hydrogenophaga ASM1878007', 'Hydrogenophaga ASM1878038', "Unclassified Thermodesulfovibrionia"),
  `Wood-Ljungdahl` = c(2, 2, 2, 0, 0, 0, 0, 2),   
  `Glycolysis` = c(2, 2, 2, 2, 2, 2, 1, 2), 
  `TCA Cycle` = c(1, 1, 1, 2, 2, 2, 2, 1),
  `Beta oxidation of fatty acids` = c(0, 0, 0, 2, 2, 2, 2, 0),
  `CBB Cycle` = c(0, 0, 1, 2, 2, 2, 2, 0),
  `Carbohydrate degradation` = c(1, 2, 2, 0, 0, 0, 1, 0),
  `Sulfate reduction` = c(2, 2, 2, 0, 0, 0, 0, 0),
  `Sulfur oxidation` = c(0, 0, 0, 2, 2, 2, 2, 0),
  `Sulfur reduction` = c(0, 0, 0, 0, 0, 0, 0, 2),
  `Denitrification` = c(0, 0, 0, 2, 2, 2, 2, 0),
  `Hydrogen oxidation` = c(2, 2, 2, 2, 2, 2, 0, 1),
  `Oxidative stress response` = c(2, 1, 1, 2, 2, 2, 2, 2),
  `CRISPR` = c(2, 2, 2, 2, 2, 0, 0, 2),
  `Restriction/Modification` = c(2, 2, 2, 2, 2, 2, 0, 2),
  `Sporulation` = c(2, 2, 2, 0, 0, 0, 0, 0)
)


# Convert to long format and replace dots with spaces in pathway names
functions_long <- functions_df %>%
  pivot_longer(cols = -organism, names_to = "pathway", values_to = "presence") %>%
  mutate(pathway = gsub("\\.", " ", pathway))

# Define colors for families and pathways
family_colors <- c("Family1" = "#FF7F00", "Family2" = "#FB9A99", "Family3" = "yellowgreen") 
group_colors <- c("darkorchid4","goldenrod1","forestgreen","firebrick","deepskyblue3")
#pathway_colors <- pals::cols25(14)
pathway_colors <- rep(group_colors, times = c(4,1,1,3,6))


# Add family information and ensure the organism order is preserved
functions_long <- functions_long %>%
  mutate(family = ifelse(grepl("^Ca. D.", organism), "Family1",
                         ifelse(grepl("^Hydrogenophaga", organism), "Family2",
                                ifelse(organism == "Unclassified Thermodesulfovibrionia", "Family3", "OtherFamily"))),
         presence = factor(presence, levels = c(0, 1, 2)),
         organism = factor(organism, levels = c('Ca. D. onstottii', 'Ca. D. sp002294005', 'Ca. D. audaxviator', 
                                                'Hydrogenophaga bedretto', 'Hydrogenophaga sp018333605', 
                                                'Hydrogenophaga ASM1878007', 'Hydrogenophaga ASM1878038', 'Unclassified Thermodesulfovibrionia')),
         pathway = factor(pathway, levels = c('Wood Ljungdahl','Glycolysis','TCA Cycle','Beta oxidation of fatty acids',
                                              'CBB Cycle','Carbohydrate degradation','Sulfate reduction',
                                              'Sulfur oxidation', 'Sulfur reduction', 'Denitrification','Hydrogen oxidation',
                                              'Oxidative stress response','CRISPR','Restriction Modification','Sporulation')))

functions_long <- functions_long %>%
  mutate(pathway = fct_rev(factor(pathway)))

# Plot 
functions_long2 <- functions_long[-c(106:120),]
functions_long2 <- functions_long2 %>%
  filter(pathway != "Sulfur reduction")
pathway_colors2 <- rep(group_colors, times = c(4,1,1,2,6))

ggplot(functions_long2, aes(x = organism, y = pathway, fill = pathway)) +
  geom_tile(data = functions_long2 %>% filter(presence == 2), color = "black", size = 0.5) +
  geom_tile_pattern(data = functions_long2 %>% filter(presence == 1), aes(pattern = presence), color = "black", size = 0.5) +
  geom_tile(data = functions_long2 %>% filter(presence == 0), fill = "white", color = "black", size = 0.5) +
  scale_pattern_manual(values = c(`1` = "stripe")) +
  scale_fill_manual(values = pathway_colors2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 1, face = "bold", color = rep(family_colors, times = c(3, 4, 1))),
        axis.text.y = element_text(size = 10, face = "bold", color = rep(group_colors, times = c(4,1,1,2,6))),
        axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"),
        legend.position = "none") +
  scale_x_discrete(position = "top", labels = levels(functions_long2$organism))

########### FIGURE S.1.2 - Hydrochemistry ###########

###########A###########
dataTrip <- read.csv("../data/IncubationResultsTrips.csv", sep = ";")
dataTrip_no1L <- dataTrip[c(1:63),]

cbPalette <- c("#FFCC00", "#56B4E9", "#009E73","#CC0000", "#990099")

# Convert 'Group' and 'Experiment' and 'Sample'to factor variables
dataTrip_no1L$Group <- factor(dataTrip_no1L$Group, levels = c("t0",  
                                                              "W1", "W2", "W3", "W4", "W5", 
                                                              "G1", "G2", "G3", "G4", "G5",
                                                              "TAn1", "TAn2", "TAn3", "TAn4", "TAn5", 
                                                              "TO1", "TO2", "TO3", "TO4", "TO5"))  
#dataTrip_no1L$Group <- factor(dataTrip_no1L$Group, levels = c("t0", "TAn1", "TO1", "G1", "W1", "TAn2", "TO2", "G2", "W2", "TAn3", "TO3", "G3", "W3", "TAn4", "TO4", "G4", "W4", "TAn5", "TO5", "G5", "W5")) dataTrip_no1L$Group <- factor(dataTrip_no1L$Group, levels = c("t0", "TAn1", "TO1", "G1", "W1", "TAn2", "TO2", "G2", "W2", "TAn3", "TO3", "G3", "W3", "TAn4", "TO4", "G4", "W4", "TAn5", "TO5", "G5", "W5")) # Group by time point
dataTrip_no1L$Name <- factor(dataTrip_no1L$Name, levels = c("t0", "Water_24hr", "Water_48hr", "Water_96hr", "Water_168hr", "Water_240hr", 
                                                            "Glass_24hr", "Glass_48hr", "Glass_96hr", "Glass_168hr", "Glass_240hr", 
                                                            "RockAnoxic_24hr", "RockAnoxic_48hr", "RockAnoxic_96hr", "RockAnoxic_168hr", "RockAnoxic_240hr", 
                                                            "RockOxic_24hr", "RockOxic_48hr", "RockOxic_96hr", "RockOxic_168hr", "RockOxic_240hr"))  
dataTrip_no1L$Experiment <- factor(dataTrip_no1L$Experiment)
dataTrip_no1L$Sample <- factor(dataTrip_no1L$Sample)

# Define the desired order of legend items
experiment_order <- c("Initial", "Water", "Glass", "RockAnoxic", "RockOxic")
# Set the levels of the 'Experiment' variable to the desired order
dataTrip_no1L$Experiment <- factor(dataTrip_no1L$Experiment, levels = experiment_order)


# Create a grouped barplot with triplicates split into separate bars
ggplot(dataTrip_no1L, aes(x = Name, y = EC, fill = Experiment, group = Sample)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Sample", y = "Electrical Conductivity (μS/cm)", title = "Electrical Conductivity") +
  theme_bw() +
  scale_fill_manual(values = cbPalette) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(size=14),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(size=16, margin = margin(t = 15)),
        axis.title.y = element_text(size=16, margin = margin(r = 15))) +
  guides(fill = guide_legend(title = "Treatment")) +
  coord_cartesian(ylim = c(400, 750))

###########B###########
ggplot(dataTrip_no1L, aes(x = Name, y = pH, fill = Experiment, group = Sample)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Sample", y = "pH", title = "pH") +
  theme_bw() +
  scale_fill_manual(values = cbPalette) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(size=14),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(size=16, margin = margin(t = 15)),
        axis.title.y = element_text(size=16, margin = margin(r = 15))) +
  guides(fill = guide_legend(title = "Treatment")) +  
  coord_cartesian(ylim = c(9, 10.25))

###########C###########
ggplot(dataTrip_no1L, aes(x = Name, y = ORP, fill = Experiment, group = Sample)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Sample", y = "ORP (mV)", title = "Redox Potential") +
  theme_bw() +
  scale_fill_manual(values = cbPalette) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(size=14),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(size=16, margin = margin(t = 15)),
        axis.title.y = element_text(size=16, margin = margin(r = 15))) +
  guides(fill = guide_legend(title = "Treatment"))  +  
  coord_cartesian(ylim = c(-180, -100))

###########D###########
sulf <- read.csv("../data/sulfate.csv", sep = ",")
sulf_no1L <- sulf[-c(7,13,19,25),]
cbPalette <- c("#FFCC00", "#56B4E9", "#009E73","#CC0000", "#990099" )

sulf_no1L$Type<- factor(sulf_no1L$Type, levels = type_order)

ggplot(sulf_no1L, aes(x = as.factor(Time), y = SO4_ppm, fill = Type)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.8, color="black") +
  geom_errorbar(aes(ymin = SO4_ppm - StdDev, ymax = SO4_ppm + StdDev), position = position_dodge(width = 1), width = 0.5) +
  labs(x = "Time (Hours)", y = "Sulfate (ppm)", title = "Sulfate - Incubation V2") +
  facet_wrap(~Type, nrow = 1) +
  theme_bw() +
  scale_fill_manual(values=cbPalette) +
  theme(legend.position = "none", 
        axis.title.y = element_text(size=16,margin = margin(r=15)),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=11),
        axis.title.x = element_text(size=16, margin = margin(t = 15)),
        strip.text.x = element_text(size = 16))+
  scale_y_continuous(labels = function(x) format(x, scientific = F)) +  
  coord_cartesian(ylim = c(150, 250))