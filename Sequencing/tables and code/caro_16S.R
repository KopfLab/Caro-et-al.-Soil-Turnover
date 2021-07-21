library(mctoolsr)

tax_table_fp = 'seqtab_wTax_mctoolsr_16s_caro.txt'
map_fp = 'caro_16s_mapping.txt' 
input = load_taxa_table(tax_table_fp, map_fp)

# Remove mitochondrial and chloroplast sequences, Remove reads assigned as eukaryotes
input_filt <- filter_taxa_from_input(input, taxa_to_remove = c("Chloroplast","Mitochondria", "Eukaryota"))
# Remove reads that are unassigned at domain level
input_filt <- filter_taxa_from_input(input_filt, at_spec_level = 2, taxa_to_remove = "NA")

# Normalize or rarefy your ASV table
#input_filt <- single_rarefy(input = input_filt, depth = 5000) # CHANGE ME to desired depth.

#### Now following Cliff's R Tutorial for Alpha & Beta Diversity and other Visualization Techniques
library(plyr) # For hulls
library(tidyverse) # For data manipulation and plotting
library(RColorBrewer) # For color palettes
library(vegan) # For alpha and beta diversity analyses
library(indicspecies) # For indicator species
library(car) # For Levene Test and anova
library(PMCMR) # For Nemenyi posthoc test


# We now have an object called input. It contains 3 files - a sequence table, a mapping file, and a taxonomic file.
head(input$data_loaded) # Rows are OTUs, columns are samples
head(input$taxonomy_loaded) # Rows are OTUs (now ASVs), columns are different taxonomic levels
head(input$map_loaded) # Rows are samples, columns are variables

# Initial data examination, rarefying, filtering
# One of the first things we want to do is see how many sequences per sample there are
# This is done by getting the column sums of the sequence table, and we'll sort it too.
sort(colSums(input$data_loaded))
# We could also save this as an object and plot it
seqcounts <- as.data.frame(sort(colSums(input$data_loaded)))
# This puts the sample names as row names and the sequence counts column is titled "sort(colSums(input$data_loaded)))
# We'll use the pipeline %>% which is very useful for dataframe management with the dplyr package
# We'll rename the column and make a new column all at once
seqcounts <- as.data.frame(sort(colSums(input$data_loaded))) %>%
  rename("seqs" = "sort(colSums(input$data_loaded))") %>%
  rownames_to_column(var = "sampleID")
# Now we have a better dataframe with two columns, seqs and sampleID which we can plot
ggplot(seqcounts, aes(reorder(sampleID, seqs, mean), seqs)) + # Dataframe and variables
  geom_bar(stat = "identity") + # Type of graph
  labs(y = "# Reads", x = "Sample") + # Axes labels
  coord_flip() + # Flip axes
  theme_classic() + # Plot style. I also use theme_bw() a lot.
  theme(axis.text.y = element_text(size = 2)) # Tweak things about the text and legend in theme()

seqcounts
# This shows that we have anywhere from 27,000-60,000 reads per sample

# Moving on, let's rarefy the data at the lowest count per sample. Sometimes you may want to drop samples and choose a higher number.
# Rarefy at 27,000 seqs/sample
input_rar = single_rarefy(input, 27000) # This makes a new mctoolsr dataset called input_rar. The blank sample dropped out!

# Check seq counts in new dataset
colSums(input_rar$data_loaded) # Good - all 27000!


# Filtering - we can filter out samples or taxa
# We often filter out anything not assigned to a phylum and any control samples. You may also want to do analyses on only a subset of your data.
# For 16S, filter out mitochondrial, chloroplast, eukaryote DNA, and things unassigned at phylum level 

# Done above 

# Note: normally after you've rarefied and filtered it's recommended to save the dataset as a .rds file. This insures you use the same dataset and analysis is reproducible (because rarefaction involves random sampling). Then, if you quit your session but want to reanalyze something, you would ignore the rarefaction/filtering and just start right here!
saveRDS(input_rar, file = "input_rar.rds")
input_rar <- readRDS("input_rar.rds")


#### Alpha Diversity ####
# not extremely useful, but people report it
# Let's analyze and graph the number of OTUs (now ASVs) in our sample types
# Two categories: t-test (parametric) or Wilcoxon Test (non-parametric)
# Three + categories: ANOVA (parametric) or Kruskal-Wallis Test (non-parametric)

# First let's get the number of ASVs (richness) per sample and add it to the mapping file
input_rar$map_loaded$rich <- specnumber(input_rar$data_loaded, MARGIN = 2)
# Note: since the data_loaded file has ASVs as rows, MARGIN must be set to 2.


# Three or more categories
leveneTest(input_rar$map_loaded$rich ~ input_rar$map_loaded$combo_label)
# Variance homogeneous (p > 0.05)

m <- aov(input_rar$map_loaded$rich ~ input_rar$map_loaded$combo_label)
shapiro.test(m$residuals)
# Residuals normally distributed (p > 0.05)

# Other diagnostics - learn more here (https://data.library.virginia.edu/diagnostic-plots/)
plot(m) # Click in the console and hit Return to see each diagnostic plot
summary(m)
# Significant effect. Run posthoc test to see which groups are different from each other.
# Significant effect of ASV richness on soil type, but which ones (six groups possible)
TukeyHSD(m)
# the degraded and pitted soils had different richness in their soils vs. their slurries. Healthy ones did not. 

# Plot. Just copy the previous plot, change the variable name and the colour argument!
ggplot(input_rar$map_loaded, aes(combo_label, rich, colour = combo_label)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.5) +
  labs(x = "Soil Type", y = "ASV Richness", colour = "combo_label") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14),
        legend.position = "none")

# Cliff's Tutorial shows how to add labels to the chart to indicate significant differences. Add these later

#### Beta Diversity ############
# Let's calculate Bray-Curtis dissimilarity, plot some PCoA ordinations, and run some stats.
# Bray-Curtis Matrix
bc <- calc_dm(input_rar$data_loaded)
# can also do NMDS plots

# Principle Coordinates Analysis (PCoA)
pcoa <- cmdscale(bc, k = nrow(input_rar$map_loaded) - 1, eig = T)
# Variation Explained
eigenvals(pcoa)/sum(eigenvals(pcoa)) # 25.7, 18.9 % variation explained
# can plot in mctoolsr but Cliff prefers ggplots, so save the axes info

# Save Axis 1 and 2 to the mapping file
input_rar$map_loaded$Axis01 <- scores(pcoa)[,1]
input_rar$map_loaded$Axis02 <- scores(pcoa)[,2]

# Function for making a convex hull
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]

# Calculate hulls and save to dataframe
micro.hulls <- ddply(input_rar$map_loaded, c("combo_label"), find_hull)
ggplot(input_rar$map_loaded, aes(Axis01, Axis02, colour = combo_label)) +
  geom_polygon(data = micro.hulls, aes(colour = combo_label),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "PC1: 25.7% Variation Explained", 
       y = "PC2: 18.9% Variation Explained",
       colour = "Soil Type") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))

# PERMANOVA
# test for differences in centroid of different pCoA hull groups
set.seed(1223) # To make reproducible 
adonis(bc ~ input_rar$map_loaded$combo_label)
# Significant effects of type.

# PERMDISP
# multivariate version of Levenne Test. Difference in variation in each factor level
m1 <- betadisper(bc, input_rar$map_loaded$combo_label)
anova(m1) # Dispersion homogeneous for soil type


#### Taxonomic Analyses ####
# Indicator taxa - SIMPER (Similarity Percentages) or MULTIPATT
# both take into account the number of samples that the taxa is present in within each group as well as abundance
# SIMPER (list how much each ASV contributes to dissimilarity among groups)
sim <- simper(t(input_rar$data_loaded), 
              input_rar$map_loaded$combo_label)
s <- summary(sim)
# Let's look at the top 5 contributing to dissimilarity between the two types with significantly different soils
head(s$GGm0_GGm3, n = 5)

# average is the proportion contribution, cumsum is cumulative, ava and avb are mean sequence abundances per group.

# MULTIPATT (list ASVs associated with each group)
# difference between this and SIMPER is in what it is reporting
# MULTIPATT doesn't give a value for every OTU. It only outputs important OTUs
set.seed(1223) # For reproducibility
mp <- multipatt(t(input_rar$data_loaded), 
                input_rar$map_loaded$combo_label, 
                func = "IndVal.g", 
                control = how(nperm=999))
summary(mp)

# Main taxa at different levels
# Let's quickly just use the mctoolsr base functions
# Summarize at family level. Choose different levels with the level argument
tax_sum_families <- summarize_taxonomy(input_rar, level = 5, report_higher_tax = FALSE) #level 5 = families, level 6 = genus

# mctoolsr Default heatmap, which uses ggplot2 in the background so we can use normal ggplot2 syntax
# numbers are percent relative abundances
plot_ts_heatmap(tax_sum_families, 
                input_rar$map_loaded, 
                0.01, # minimum abundance threshold
                'combo_label')

# Add some customization
plot_ts_heatmap(tax_sum_families, 
                input_rar$map_loaded, 
                0.01, 
                'combo_label',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1))

# Stacked bar plot for mean abundance levels
# Default (choose how many taxa you want with num_taxa)
plot_taxa_bars(tax_sum_families,
               input_rar$map_loaded,
               "combo_label",
               num_taxa = 10) # this limits the number of taxa shown in the barplot and the rest go into "other"

# Add some customization
plot_taxa_bars(tax_sum_families,
               input_rar$map_loaded,
               "combo_label",
               num_taxa = 10) +
  labs(x = "Soil Type", y = "Relative Abundance", fill = "Family") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))

# Add even more customization. plot_taxa_bars has a data_only argument. So you can save the dataset created and then write your own plot.
bars <- plot_taxa_bars(tax_sum_families,
                       input_rar$map_loaded,
                       "combo_label",
                       num_taxa = 10,
                       data_only = TRUE) # allows us to save the data and not make the plot
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = "black", size = 0.25) +
  labs(x = "Soil Type", y = "Relative Abundance", fill = "Family") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))

# Test (run a Kruskal-Wallis test on all families with mean rel abund. > filter_level in at least one of the factor levels)
# which families are different between the sample types (testing all at once)
taxa_summary_by_sample_type(tax_sum_families, 
                            input_rar$map_loaded, 
                            type_header = 'combo_label', 
                            filter_level = 0.05, # cut off for relative abundance
                            test_type = 'KW') # Bonferroni corection imporant since you did many tests at once

# Two of these are significantly different...and they are because of the desert soils

# Next steps...filtering the ASV table to get rid of the NA and other categories.

# Trying at the Phylum level = 2
# Summarize at family level. Choose different levels with the level argument
tax_sum_phyla <- summarize_taxonomy(input_rar, level = 2, report_higher_tax = FALSE) #level 5 = families, level 6 = genus

# mctoolsr Default heatmap, which uses ggplot2 in the background so we can use normal ggplot2 syntax
# numbers are percent relative abundances
plot_ts_heatmap(tax_sum_phyla, 
                input_rar$map_loaded, 
                0.01, # minimum abundance threshold
                'combo_label')

# Add some customization
plot_ts_heatmap(tax_sum_phyla, 
                input_rar$map_loaded, 
                0.01, 
                'combo_label',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1))

# Stacked bar plot for mean abundance levels
# Default (choose how many taxa you want with num_taxa)
plot_taxa_bars(tax_sum_phyla,
               input_rar$map_loaded,
               "combo_label",
               num_taxa = 10) # this limits the number of taxa shown in the barplot and the rest go into "other"

# Add some customization
plot_taxa_bars(tax_sum_phyla,
               input_rar$map_loaded,
               "combo_label",
               num_taxa = 10) +
  labs(x = "Soil Type", y = "Relative Abundance", fill = "Phyla") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))

# Add even more customization. plot_taxa_bars has a data_only argument. So you can save the dataset created and then write your own plot.
bars <- plot_taxa_bars(tax_sum_phyla,
                       input_rar$map_loaded,
                       "combo_label",
                       num_taxa = 10,
                       data_only = TRUE) # allows us to save the data and not make the plot

ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = "black", size = 0.25) +
  labs(x = "Soil Type", y = "Relative Abundance", fill = "Phyla") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))

# Test (run a Kruskal-Wallis test on all families with mean rel abund. > filter_level in at least one of the factor levels)
# which phyla are different between the sample types (testing all at once)
taxa_summary_by_sample_type(tax_sum_phyla, 
                            input_rar$map_loaded, 
                            type_header = 'combo_label', 
                            filter_level = 0.05, # cut off for relative abundance
                            test_type = 'KW') # Bonferroni correction important since you did many tests at once

# Verrucomicrobia is significantly different

########### Graphing individual samples, not averages ############################

# Main taxa at different levels
# Let's quickly just use the mctoolsr base functions
# Summarize at family level. Choose different levels with the level argument
tax_sum_phyla <- summarize_taxonomy(input_rar, level = 2, report_higher_tax = FALSE) #level 5 = families, level 6 = genus

# mctoolsr Default heatmap, which uses ggplot2 in the background so we can use normal ggplot2 syntax
# Add some customization
plot_ts_heatmap(tax_sum_phyla, 
                input_rar$map_loaded, 
                0.01, 
                'names',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1))

# Stacked bar plot for abundance levels
# Add even more customization. plot_taxa_bars has a data_only argument. So you can save the dataset created and then write your own plot.
bars <- plot_taxa_bars(tax_sum_phyla,
                       input_rar$map_loaded,
                       "names",
                       num_taxa = 11,
                       data_only = TRUE) # allows us to save the data and not make the plot
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = "black", size = 0.25) +
  labs(x = "Soil Type", y = "Relative Abundance", fill = "Phyla") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14, angle = 90))