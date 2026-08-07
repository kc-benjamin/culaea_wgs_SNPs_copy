# Define required packages
required_packages <- c("ggplot2", "dplyr", "zoo", "tidyr")

# Check for missing packages
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

# Install missing ones only
if(length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://r-project.org")
}

# Load all packages safely
invisible(lapply(required_packages, library, character.only = TRUE))

# Set working directory
setwd("/scratch/kcb95328/ShundaLakeBrooks/11_read_depth")

#read in data (make for loop)
sex_shunda <- read.csv("/home/kcb95328/Info-Shunda/SraRunTable-metadata-SL.csv", header=TRUE)
sex_shunda <- sex_shunda[, c("Run", "sex")]

for(i in sex_shunda$Run) {
  df <- read.delim(paste0(i, "_chr20_11.7to19.2_read-depth.txt"), header = TRUE)

  infile <- paste0(i, "_chr20_11.7to19.2_read-depth.txt")
  message("Reading: ", infile)

  #reformat data
  colnames(df) <- gsub("X.scratch.kcb95328.ShundaLakeBrooks.06_bam_files.", " ", colnames(df))
  names(df) <- gsub(".trimmed.fastq.gz.sorted.bam.bai", " ", names(df))
  colnames(df) <- trimws(colnames(df))

  #reformat for combination with sex
  depth_long <- pivot_longer(df, cols = starts_with("SRR"))
  depth_long$POS <- as.numeric(depth_long$POS)
  depth_long$value <- as.numeric(depth_long$value)
  rm(df)

  #normalize depth values
  depth_long$value <- depth_long$value / mean(depth_long$value)

  #merge
  names(depth_long)[3] <- "Run"
  df_with_sex <- merge(depth_long,sex_shunda, by = "Run")

  #smooth out the data:
  window_size <- 10000

  df_with_sex$smooth <- rollmean(df_with_sex$value, k = window_size, fill = NA, align = "center")

  #plot (males only)
  plot <- ggplot(df_with_sex, aes(x = POS)) +
    geom_line(aes(y = smooth), color = "blue", linewidth = 1) + # Smoothed line
    labs(
      title = paste("Read depth from 11M to 19M", i, "Normalized 10,000 bp windows"),
      x = "Position (bp)",
      y = "Depth")
  #save it
  pdf(paste0("Norm10000-Shunda-chr20-11to19", i, ".pdf"))
  print(plot)
  dev.off()
}
#plot (females only)
# femaleplot<-ggplot(df_female, aes(x = POS)) +
#   geom_line(aes(y = smooth), color = "blue", linewidth = 1) + # Smoothed line
#   labs(title = "Read depth from 11M to 19M, Shunda females only - Normalized 10,000 bp windows", x = "Position (bp)", y = "Depth")
# #save it
# pdf("Norm10000-Shunda-chr20-11to19-all_8-2026.pdf")
# print(femaleplot)
# dev.off()

# #let's get rid of that crazy region:
# #define threshold
# threshold <- 3
# #keep only those below the threshold
# small_df_with_sex<- df_with_sex %>%
#     filter(value <= threshold)

# #plot
# thresholdplot<-ggplot(small_df_with_sex, aes(x = POS)) +
#   geom_line(aes(y = smooth), color = "blue", linewidth = 1) + # Smoothed line
#   labs(title = "Sliding Window Smoothing (Rolling Mean)", x = "Position (bp)", y = "Depth") +
#   theme(text = element_text(family = "Arial"))
# #save
# pdf("Threshold3-Norm10000-Shunda-chr20-11to19-all_8-2026.pdf")
# print(thresholdplot)
# dev.off()