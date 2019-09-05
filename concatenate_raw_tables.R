library(data.table)
library(plyr)

setwd("/Users/av1936/Desktop/Google_Drive/collaborative_projects/blaser_alzheimers/feature_counts")

files <- list.files(pattern = "feature_counts")
print(files)
temp <- lapply(files, read.table, header=T)

for (i in 1:length(temp)){
  df <- temp[[i]]
  df <- df[,c(1,7)]
  names(df) <- c("genenames", "counts")
  temp[[i]] <- df
}

first <- temp[[1]]
for (i in 2:length(temp)){
  first <- merge(first, temp[[i]], by="genenames")
}

files = substr(files,1,nchar(files)-59)
names(first) <- c("genenames",files)

write.table(first,"feature_counts.txt",col.names=T, row.names=F,
            sep="\t", quote=F)


