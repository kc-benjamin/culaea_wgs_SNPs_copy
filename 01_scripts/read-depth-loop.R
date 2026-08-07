#testing script
setwd("/scratch/kcb95328/ShundaLakeBrooks/11_read_depth")

#load packages
required_packages <- c("ggplot2", "dplyr", "zoo", "tidyr")

missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if(length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://r-project.org")
}

invisible(lapply(required_packages, library, character.only = TRUE))

#read in data (make for loop)
sex_shunda <- read.csv("/home/kcb95328/Info-Shunda/SraRunTable-metadata-SL.csv", header=TRUE)
sex_shunda <- sex_shunda[, c("Run", "sex")]

for(i in sex_shunda$Run) {
  df <- read.delim(paste0(i, "_chr20_11.7to19.2_read-depth.txt"), header = TRUE)

#no need to reformat for combo with sex
#make sure its in the right orientation and the column names are good
  colnames(df) <- gsub("X.scratch.kcb95328.NoFiltering-Shunda.01_aligned_bams.", " ", colnames(df))
  names(df) <- gsub(".trimmed.fastq.gz.sorted.bam.bai", " ", names(df))
  colnames(df) <- trimws(colnames(df))

  df$POS <- as.numeric(df$POS)
  df[, 3] <- as.numeric(df[, 3])

#normalize by average read depth across the chromosome
  df[, 3] <- df[, 3] / mean(df[, 3])

#smooth out by averaging across 15,000bp windows
  window_size <- 15000

  df$smooth <- rollmean(df[, 3], k = window_size, fill = NA, align = "center")

#cut down the dataset to 15 to 16,000,000
  df_15to16 <- df[df$POS >= 15000000 & df$POS <= 16000000, ]

#plot
  plot <- ggplot(df_15to16, aes(x = POS)) +
    geom_line(aes(y = smooth), color = "blue", linewidth = 1) + # Smoothed line
    labs(
	    title = paste("Read depth from 15M to 16M", i, "Normalized 15,000 bp windows"),
	    x = "Position (bp)",
	    y = "Depth")

  pdf(paste0("Norm15000-NFShunda-chr20-15to16-", i, ".pdf"))
  print(plot)
  dev.off()
}
