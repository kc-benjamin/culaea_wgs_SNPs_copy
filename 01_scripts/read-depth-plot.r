install.packages(c("zoo", "ggplot2", "tidyr", "dplyr", "rollmean"), repos="https://cran.r-project.org")

library('ggplot2')
library('tidyr')
library('zoo')
library('rollmean')
library('dplyr')
setwd("/scratch/kcb95328/ShundaLakeBrooks/11_read_depth")

#read in data
df <- read.delim("Shunda-11to19-chr20-depth.txt")

#reformat data
colnames(df) <- gsub("X.scratch.kcb95328.ShundaLakeBrooks.09_bam_copies."," ",colnames(df)) 
names(df) <- gsub(".trimmed.fastq.gz.sorted.bam"," ",names(df))
colnames(df) <- trimws(colnames(df))

#reformat for combination with sex
depth_long<-pivot_longer(df,cols=starts_with("SRR"))
depth_long$POS<-as.numeric(depth_long$POS)
depth_long$value<-as.numeric(depth_long$value)

#read in sex metadata
sex_shunda<-read.csv("/home/kcb95328/Info-Shunda/SraRunTable-metadata-SL.csv", header=TRUE)
sex_shunda<-sex_shunda[, c("Run","sex")]

#normalize depth values
depth_long$value <- depth_long$value / mean(depth_long$value)

#merge
names(depth_long)[3]<-"Run"
df_with_sex<-merge(depth_long,sex_shunda, by="Run")

#smooth out the data:
window_size <- 10000

df_with_sex$smooth <- rollmean(df_with_sex$value, k = window_size, fill = NA, align = "center")

#plot
largeplot<-ggplot(df_with_sex, aes(x = POS)) +
  geom_line(aes(y = value), color = "gray80", alpha = 0.5) +  # Raw data in the background
  geom_line(aes(y = smooth), color = "blue", size = 1) + # Smoothed line
  labs(title = "Sliding Window Smoothing (Rolling Mean)", x = "Position (bp)", y = "Depth")
#save it
ggsave(filename="Norm10000-Shunda-chr20-11to19-all_8-2026.pdf", plot=largeplot)

#let's get rid of that crazy region:
#define threshold
threshold <- 3
#keep only those below the threshold
small_df_with_sex<- df_with_sex %>%
    filter(value <= threshold)

#plot
thresholdplot<-ggplot(small_df_with_sex, aes(x = POS)) +
  geom_line(aes(y = smooth), color = "blue", size = 1) + # Smoothed line
  labs(title = "Sliding Window Smoothing (Rolling Mean)", x = "Position (bp)", y = "Depth") +
  theme(text = element_text(family = "Arial"))
#save
ggsave(filename="Threshold3-Norm10000-Shunda-chr20-11to19-all_8-2026.pdf", plot=thresholdplot)
