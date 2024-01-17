#' Nucleotide information of 3,138 virus segments
#'
#' A data set containing the GC content and other information about the virus
#' segments from the official Virtool virus data base. The variables are as
#' follows:
#'
#' @format A data frame with 3,138 rows and 7 variables:
#' \describe{
#'   \item{virus_name}{The virus name}
#'   \item{iso_id}{The virus isolate ID}
#'   \item{seg_id}{The virus segment ID}
#'   \item{A_percent}{The percentage of A nucleotides in the virus segment}
#'   \item{C_percent}{The percentage of C nucleotides in the virus segment}
#'   \item{T_percent}{The percentage of T nucleotides in the virus segment}
#'   \item{GC_percent}{The percentage of G and C nucleotides in the virus
#'   segment}
#'   \item{seg_len}{The length of the virus segment}
#'   }
"default_nucleotide_info"

#' A mappability profile for Arabidopsis.
#'
#' @format A data frame...
"default_mappability_profile"

#' The virus segments from the official Virtool virus data base
#'
#' A DNAStringSet object of the 3,138 virus segments.
#'
#'
"virus_segments"

#' A trained model using the default Random Forest settings
#'
#'

"trained_rf"



#' A trained model using the default XGBoost settings
#'
#'
#'

"trained_xgb"


#' Known diagnostics result of over 3,000 virus segments
#'
#' A data set containing the known truth about the diagnostics result for each
#' plant sample. It records whether the sample is infected with a virus segment.
#'
#' @format A data frame with 3,138 rows and 21 columns:
#' \describe{
#'   \item{row}{Each row is the name of a virus segment}
#'   \item{column}{Each column is the name of a plant sample}
#'   }
"example_diag"


#' Coverage profile of 21 plant samples.
#'
#' A list of coverage profiles for 21 plant samples. This is only a toy sample.
#' You can use it for running the examples in the vignette. We recommend using
#' more data to train the model, the more the better.
#'
#' @format A data frame with 3,138 rows and 21 columns:
#' \describe{
#'   \item{row}{Each row is the name of a virus segment}
#'   \item{column}{Each column is the name of a plant sample}
#'   }
"example_cov"
