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
  colnames(df) <- gsub("X.scratch.kcb95328.ShundaLakeBrooks.06_bam_files.", " ", colnames(df))
  names(df) <- gsub(".trimmed.fastq.gz.sorted.bam.bai", " ", names(df))
  colnames(df) <- trimws(colnames(df))

  df$POS <- as.numeric(df$POS)
  df[, 3] <- as.numeric(df[, 3])

#normalize by average read depth across the chromosome
  df[, 3] <- df[, 3] / mean(df[, 3])

#smooth out by averaging across 20,000bp windows
  window_size <- 20000

  df$smooth <- rollmean(df[, 3], k = window_size, fill = NA, align = "center")

#cut down the dataset to 16 to 19,202,098
  df_18to19 <- df[which.max(df$POS == 18000000):nrow(df), ]

#plot
  plot <- ggplot(df_18to19, aes(x = POS)) +
    geom_line(aes(y = smooth), color = "blue", linewidth = 1) + # Smoothed line
    labs(
	title = paste("Read depth from 18M to 19M", i, "Normalized 20,000 bp windows"),
	x = "Position (bp)",
	y = "Depth")

  pdf(paste0("Norm10000-Shunda-chr20-11to19-", i, ".pdf"))
  print(plot)
  dev.off()
}
