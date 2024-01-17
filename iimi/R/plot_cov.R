#' @title plot_cov()
#'
#' @export
#' @importFrom graphics lines plot legend
#'
#'
#'
#'
#' @description Plots the coverage profile of the mapped plant sample.
#'
#' @param covs An RLE list of coverage information of one or more plant samples.
#' @param cex The size of the main title.
#' @param nucleotide_status Whether display the percentages of A nucleotide and GC
#'     content. Default is TRUE.
#' @param legend_status Whether display the legend of the plot. Default is TRUE.
#' @param window The size of the sliding window. Default is 75.
#'
#' @return The coverage profile of the mapped plant sample.

plot_cov <- function(covs, cex = 0.7, nucleotide_status = T, legend_status = T, window = 75) {
  for (sample in names(covs)) {
    counter = 1

    for (seg in names(covs[[sample]])) {
      plot(
        covs[[sample]][[seg]],
        type = "l",
        main = paste0(
          counter,
          "). ",
          sample,
          ": ",
          nucleotide_info[nucleotide_info$seg_id == seg,1],
          ",",
          nucleotide_info[nucleotide_info$seg_id == seg,2],
          ",",
          seg
        ),
        ylab = "",
        cex.main = cex,
        xlab = ""
      )

      if (length(mappability_profile[[seg]]) > 0) {
        lines(
          x = mappability_profile[[seg]],
          y = rep(0, length(mappability_profile[[seg]])),
          col = "red",
          lwd = 2.5
        )
      }

      if (nucleotide_status) {
        gc = rowSums(
          letterFrequencyInSlidingView(
            virus_segments[[seg]],
            window,
            c("G", "C")
          )
        ) / window

        a = rowSums(
          letterFrequencyInSlidingView(
            virus_segments[[seg]],
            window,
            c("A")
          )
        ) / window

        plot(gc, type = "l", ylim = c(0,1), xaxt = 'n', xlab = "", ylab = "")

        lines(a, lty = 3, col = "dimgrey")

        if (legend_status) {
          legend(
            x = "top",
            legend = c("Untrustworthy region", "GC content", "A percentage"),
            col = c("red", "black", "dimgrey"),
            lty = c(1,1,3),
            lwd = 2,
            horiz = TRUE,
            xpd = TRUE,
            bty = "n",
            inset = c(1,0),
            cex = cex
          )
        }
      }

      counter = counter + 1
    }
  }
}
