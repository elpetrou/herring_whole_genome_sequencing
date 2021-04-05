# The purpose of this script is to plot the distribution of raw sequencing reads for each herring sample. The input data are data frames specifying the raw number of reads per sample (saved as .csv files). These data were provided to me by the NW Genomics Sequencing Center.

# Load libraries
library(tidyverse)

# Specify path to input files
MYPATH <- "./sample_metadata/"
setwd(MYPATH)
list.files()

# Specify file names
file1 <- "Eleni_Lane1_WA_herring.stats.csv"
file2 <- "Eleni_Lane2_AK_herring.stats.csv"

# Read in the data
wa_df <- read.csv(file1)
ak_df <- read.csv(file2)

plot_df <- rbind(wa_df, ak_df)

# plot the data

distro_plot <- ggplot(data = plot_df) +
  geom_histogram(aes(x = reads_millions), bins = 20) +
  #facet_wrap(~lane) +
  ylab("Number of samples") +
  xlab("Number of raw sequences (millions)") +
  theme_bw()
  
distro_plot 

# Calculate some summary statistics
(mymean <- mean(plot_df$reads_millions)) # mean = 12.88 million reads
(myrange <- range(plot_df$reads_millions)) #range = 0 to 80.06 million reads
(mysd <-sd(plot_df$reads_millions)) # sd = 4.41 million reads

# Flad some outlier samples
(lower_threshold <- mymean - 2*mysd) 
(upper_threshold <- mymean + 2*mysd)
  
outlier_df <- plot_df %>%
  filter(reads_millions < lower_threshold | reads_millions > upper_threshold )

# save output of analyses
ggsave("plot_raw_seq_distro.pdf", distro_plot)
