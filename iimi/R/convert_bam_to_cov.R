#' @title convert_bam_to_cov
#'
#' @export
#' @importFrom log_info logger
#' @importFrom Rsamtools BamFile
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom IRanges coverage
#' @importFrom stats setNames
#'
#' @description Converts one or more indexed and sorted BAM files into a
#'     feature-extracted data frame after mappability profiling
#'
#' @param bam_file path to BAM file(s)
#' @param covs a list of Coverage profile(s) in RLE format. Can be one or more samples.
#' @param mappability_profile The mappability profile from a host genome (we only have
#'     Arabidopsis thaliana right now)
#' @param paired Indicate if the sequencing paired is single-end or paired-end
#'     reads. True if paired-end. False if single-end.
#'
#' @return A data frame object that contains the mapping result for each virus
#'     segment that the plant sample reads are aligned to and a RLE list of coverage
#'     information
convert_bam_to_cov <- function(
  bam_file,
  covs,
  mappability_profile = default_mappability_profile,
  nucleotide_info = default_nucleotide_info,
  paired = F
) {
  if (missing(covs)) {
    ## first, convert BAM files to coverages
    covs <- setNames(lapply(
      bam_file,
      function(x) {
        bams <- BamFile(x)

        if (!paired) {
          xread <- readGAlignments(bams)
        } else if (paired) {
          xread <- readGAlignmentPairs(bams)
        }

        cov <- coverage(xread)

        return(cov)
      }), sub(pattern = "(.*)\\.sorted.*$", replacement = "\\1", basename(bam_file)))
  }

  ### only keep reads that are mapped to the virus segment
  covs_only_mapped <- lapply(covs, function(x) {
    x[sapply(x, function(y){!all(y==0)})]
  })

  ## MP step
  for (sample in names(covs_only_mapped)) {
    for (seg in names(covs_only_mapped[[sample]])) {
      if (length(mappability_profile[[seg]]) > 0 &&
          (length(mappability_profile[[seg]])/length(covs_only_mapped[[sample]][[seg]])) < 0.2) {
        for (ii in mappability_profile[[seg]]) {
          covs_only_mapped[[sample]][[seg]][ii] = 0
        }
      }
    }
  }

  ## convert from cov to df
  model_data <- data.frame(matrix(ncol = 20, nrow = 0))

  colnames(model_data) <- c(
    "seg_id",
    "iso_id",
    "virus_name",
    "sample_id",
    "A_percent",
    "C_percent",
    "T_percent",
    "GC_percent",
    "avg_cov",
    "max_cov",
    "seg_len",
    "cov_2_percent",
    "cov_3_percent",
    "cov_4_percent",
    "cov_5_percent",
    "cov_6_percent",
    "cov_7_percent",
    "cov_8_percent",
    "cov_9_percent",
    "cov_10_percent"
  )
  
  max_cov <- vector()
  mean_cov <- vector()

  #use for loop to get the data frame
  for (sample in names(covs_only_mapped)) {
    for (seg in names(covs_only_mapped[[sample]])) {
      iso_id = nucleotide_info[nucleotide_info$seg_id == seg, 2]
      v_name <- nucleotide_info[nucleotide_info$seg_id == seg, 1]
      
      sample_id = sample
      seg_id = seg
      seg_length = nucleotide_info[nucleotide_info$seg_id == seg, 8]
      a_content = nucleotide_info[nucleotide_info$seg_id == seg, 4]
      c_content = nucleotide_info[nucleotide_info$seg_id == seg, 5]
      t_content = nucleotide_info[nucleotide_info$seg_id == seg, 6]
      gc_content = nucleotide_info[nucleotide_info$seg_id == seg, 7]
      max_val = max(covs_only_mapped[[sample]][[seg]])
      mean_val = mean(covs_only_mapped[[sample]][[seg]])
      max_cov <- append(max_cov,max_val)
      mean_cov <- append(mean_cov, mean_val)
      idx2 <- covs_only_mapped[[sample]][[seg]]@values > 2
      idx3 <- covs_only_mapped[[sample]][[seg]]@values > 3
      idx4 <- covs_only_mapped[[sample]][[seg]]@values > 4
      idx5 <- covs_only_mapped[[sample]][[seg]]@values > 5
      idx6 <- covs_only_mapped[[sample]][[seg]]@values > 6
      idx7 <- covs_only_mapped[[sample]][[seg]]@values > 7
      idx8 <- covs_only_mapped[[sample]][[seg]]@values > 8
      idx9 <- covs_only_mapped[[sample]][[seg]]@values > 9
      idx10 <- covs_only_mapped[[sample]][[seg]]@values > 10
      percent_2 <- sum(covs_only_mapped[[sample]][[seg]]@lengths[idx2])/sum(covs_only_mapped[[sample]][[seg]]@lengths)
      percent_3 <- sum(covs_only_mapped[[sample]][[seg]]@lengths[idx3])/sum(covs_only_mapped[[sample]][[seg]]@lengths)
      percent_4 <- sum(covs_only_mapped[[sample]][[seg]]@lengths[idx4])/sum(covs_only_mapped[[sample]][[seg]]@lengths)
      percent_5 <- sum(covs_only_mapped[[sample]][[seg]]@lengths[idx5])/sum(covs_only_mapped[[sample]][[seg]]@lengths)
      percent_6 <- sum(covs_only_mapped[[sample]][[seg]]@lengths[idx6])/sum(covs_only_mapped[[sample]][[seg]]@lengths)
      percent_7 <- sum(covs_only_mapped[[sample]][[seg]]@lengths[idx7])/sum(covs_only_mapped[[sample]][[seg]]@lengths)
      percent_8 <- sum(covs_only_mapped[[sample]][[seg]]@lengths[idx8])/sum(covs_only_mapped[[sample]][[seg]]@lengths)
      percent_9 <- sum(covs_only_mapped[[sample]][[seg]]@lengths[idx9])/sum(covs_only_mapped[[sample]][[seg]]@lengths)
      percent_10 <- sum(covs_only_mapped[[sample]][[seg]]@lengths[idx10])/sum(covs_only_mapped[[sample]][[seg]]@lengths)

      # Create a new row as a list
      new_row <- list(
        seg_id = seg_id,
        iso_id = iso_id,
        virus_name = v_name,
        sample_id = sample_id,
        A_percent = a_content,
        C_percent = c_content,
        T_percent = t_content,
        GC_percent = gc_content,
        avg_cov = mean_val,
        max_cov = max_val,
        seg_len = seg_length,
        cov_2_percent = percent_2,
        cov_3_percent = percent_3,
        cov_4_percent = percent_4,
        cov_5_percent = percent_5,
        cov_6_percent = percent_6,
        cov_7_percent = percent_7,
        cov_8_percent = percent_8,
        cov_9_percent = percent_9,
        cov_10_percent = percent_10
      )

      # Append the new row to the data frame
      model_data <- rbind(model_data, new_row)
    }
  }
  
  for (i in colnames(model_data[, -c(1:4)])) {
    model_data[[i]] = as.numeric(model_data[[i]])
  }

  list(ML_df = model_data, cov = covs_only_mapped)
}
