library(Biostrings)
library(iimi)
library(mltools)
library(data.table)
library(dplyr)
library(randomForest)
library(Rsamtools)
library(GenomicAlignments)

args <- commandArgs(trailingOnly=TRUE)

bam_path <- args[1]
unreliable_regions_path <- args[2]
model_path <- args[3]
sequence_info_path <- args[4]
output_path <- args[5]

unreliable_regions <- read.csv(unreliable_regions_path, header=TRUE, check.names=FALSE)
model <- readRDS(model_path)
nucleotide_info <- read.csv(sequence_info_path, header=TRUE)

# Convert BAM file.
rle_data <- convert_bam_to_rle(bam_file = bam_path)

newdata <- convert_rle_to_df(rle_data, unreliable_regions=unreliable_regions)

par(mar = c(1,2,1,1))
layout(matrix(c(1,1,2,5,5,6,3,3,4,7,7,8), nrow = 6))

# Run Model
prediction_virus <- predict_iimi(
  newdata = newdata,
  method = "xgb",
  trained_model = model,
  report_virus_level = TRUE
)

prediction_isolate <- predict_iimi(
  newdata = newdata,
  method = "xgb",
  trained_model = model,
  report_virus_level = FALSE
)

write.csv(prediction_virus, file.path(output_path, "prediction_virus.csv"), row.names = FALSE)
write.csv(prediction_isolate, file.path(output_path, "prediction_sequence.csv"), row.names = FALSE)


write_coverage <- function(cov, path) {
  cov_table <- data.frame(matrix(ncol=3, nrow=0))

  colnames(cov_table) <- c(
    "sequence_id",
    "lengths",
    "values"
  )

  for (sequence_id in names(cov)) {
    lengths <- cov[[sequence_id]]@lengths
    values <- cov[[sequence_id]]@values

    cs_lengths <- paste(lengths, collapse = ",")
    cs_values <- paste(values, collapse = ",")

    new_row <- list(
      sequence_id=sequence_id,
      lengths = cs_lengths,
      values=cs_values
    )

    cov_table <- rbind(cov_table, new_row)
  }

  write.csv(cov_table, path, row.names = FALSE)
}


unfiltered_coverage <- coverage(readGAlignments(BamFile(bam_path)))

write_coverage(unfiltered_coverage, file.path(output_path, "coverage.csv"))

untrustworthy <- unreliable_regions %>%
   group_by(`Virus segment`) %>%
   summarise(ranges = paste(paste(Start, End, sep = "-"), collapse = ", "))
write.csv(as.data.frame(untrustworthy), file.path(output_path,"untrustworthy.csv"), row.names = FALSE)
