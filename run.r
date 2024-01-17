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
mappability_profile_path <- args[2]
model_path <- args[3]
nucleotide_info_path <- args[4]
output_path <- args[5]

mappability_profile <- readRDS(mappability_profile_path)
model <- readRDS(model_path)
nucleotide_info <- read.csv(nucleotide_info_path, sep="\t", header=TRUE)

# Convert BAM file.
data <- convert_bam_to_cov(
  bam_file = bam_path,
  mappability_profile = mappability_profile,
  nucleotide_info = nucleotide_info
)

newdata <- data$ML_df

par(mar = c(1,2,1,1))
layout(matrix(c(1,1,2,5,5,6,3,3,4,7,7,8), nrow = 6))

# Run Model
prediction_virus <- predict_iimi(
  newdata = newdata,
  method = "rf",
  trained_model = trained_rf,
  report_result_level = 2
)

prediction_isolate <- predict_iimi(
  newdata = newdata,
  method = "rf",
  trained_model = trained_rf,
  report_result_level = 1
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


write_untrustworthy <- function(mappability_profile, path) {
  unt_table <- data.frame(matrix(ncol=2, nrow=0))

  colnames(unt_table) <- c(
    "sequence_id",
    "ranges"
  )

  for (sequence_id in names(mappability_profile)) {
    number_list <- mappability_profile[[sequence_id]]

    if (length(number_list) > 0) {
      breaks <- which(diff(number_list) > 1)

      # Create the pairs
      pairs <- matrix(nrow = length(breaks) + 1, ncol = 2)

      # Fill in the start and end points of each range
      pairs[, 1] <- c(number_list[1], number_list[breaks + 1])
      pairs[, 2] <- c(number_list[breaks], number_list[length(number_list)])

      # Convert to a list of coordinate pairs
      coordinate_pairs <- split(pairs, row(pairs))

      formatted_pairs <- vector("character", length(coordinate_pairs))

      # Iterate over the coordinate pairs and format them
      for (i in seq_along(coordinate_pairs)) {
        pair <- coordinate_pairs[[i]]
        formatted_pairs[i] <- paste(pair[1], pair[2], sep = "-")
      }

      # Combine all formatted pairs into a single string, separated by commas
      stringified <- paste(formatted_pairs, collapse = ",")

      unt_table <- rbind(unt_table, list(sequence_id, stringified))
    } else {
      unt_table <- rbind(unt_table, list(sequence_id, ""))
    }
  }

  colnames(unt_table) <- c(
    "sequence_id",
    "ranges"
  )

  write.csv(unt_table, path, row.names = FALSE)

}

unfiltered_coverage <- coverage(readGAlignments(BamFile(bam_path)))

write_coverage(unfiltered_coverage, file.path(output_path, "coverage.csv"))
write_untrustworthy(mappability_profile, file.path(output_path,"untrustworthy.csv"))
