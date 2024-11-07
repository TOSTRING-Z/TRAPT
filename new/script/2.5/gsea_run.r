# example
input.rank <- as.data.frame(list("name" = c("g1", "g2", "g3", "g4", "g5"), "value" = c(1:5)))
gs.db <- list("ts1\tg1\tg2\tg3", "ts2\tg1\tg2\tg3")

source("new/script/2.5/GSEA.1.0.11.R")

result <- GSEA(
    input.rank = input.rank,
    gs.db = gs.db,
    output.directory = "new/result/2.5/",
    gs.size.threshold.min = 1
)

# run
data <- read.table("new/result/2.5/decay_rate[aggregation].txt",header = T)
head(data)
input.rank <- data[,c(1,2)]
head(input.rank)

tr_type <-read.csv("new/result/2.5/41467_2020_16106_MOESM4_ESM.csv")[c("X","assignment")]
head(tr_type)

library(stringr)
library(dplyr)

tr_type.group <- tr_type %>%
  group_by(assignment) %>%
  summarise(paste(X, collapse = "\t"))

gs.db <- sapply(1:nrow(tr_type.group), function(row) {
    paste(tr_type.group[row, 1], "", tr_type.group[row, 2], sep = "\t")
})

source("new/script/2.5/GSEA.1.0.11.R")

result <- GSEA(
    input.rank = input.rank,
    gs.db = gs.db[1],
    output.directory = "new/result/2.5/",
    gs.size.threshold.min = 1
)
file.rename(from = "Rplots.pdf", to = "new/result/2.5/gsea-long_range_TF.pdf")

result <- GSEA(
    input.rank = input.rank,
    gs.db = gs.db[2],
    output.directory = "new/result/2.5/",
    gs.size.threshold.min = 1
)
file.rename(from = "Rplots.pdf", to = "new/result/2.5/gsea-short_range_TF.pdf")
