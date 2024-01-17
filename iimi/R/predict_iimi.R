#' @title predict_iimi()
#'
#' @export
#' @importFrom dplyr %>% group_by summarise first
#' @importFrom stats predict
#' @importFrom mltools sparsify
#' @importFrom data.table data.table
#'
#'
#'
#'
#' @description Uses a machine learning model to predict the infection status
#'    for the plant sample(s). User can use their own model if needed.
#'
#' @param newdata A matrix or data frame that contains the features extracted
#'    from the coverage profile using `convert_bam_to_cov()`.
#' @param method The machine learning method of choice, either Random Forest or
#'    XGBoost. Default is the XGBoost model
#' @param trained_model The trained model. If not provided, default model is used
#' @param report_result_level Indicate what level of result the result should
#'    report. "1" reports the segment, isolate, and virus name along with the
#'    diagnostics result. "2" reports the virus names along with the diagnostics
#'    result. Default is 2.
#'
#' @return A data frame of diagnostics result for each sample


predict_iimi <- function(newdata, method, trained_model, report_result_level = 2) {
  if (method == "rf") {
    if (missing(trained_model)) {
      model = trained_rf
    } else {
      model = trained_model
    }
    
    prediction <- predict(newdata = newdata, model)
    
    pred <- data.frame(
      prediction = prediction,
      seg_id = newdata$seg_id,
      iso_id = newdata$iso_id,
      virus_name = newdata$virus_name,
      sample_id = newdata$sample_id
    )
    

    if (report_result_level == 1) {
      result_df <- pred
    } else if (report_result_level == 2) {
      result_df <- pred %>% group_by(sample_id, virus_name) %>%
        summarise(
          virus_name = dplyr::first(virus_name),
          prediction = any(as.logical(prediction))
        )
    }
  }

  if (method == "xgb") {
    test = sparsify(data.table(newdata))

    if (missing(trained_model)) {
      model = trained_xgb
    } else {
      model = trained_model
    }

    pred = predict(newdata = test, model)

    pred <- data.frame(
      prediction = pred > 0.5,
      seg_id = newdata$seg_id,
      iso_id = newdata$iso_id,
      virus_name = newdata$virus_name,
      sample_id = newdata$sample_id
    )

    if (report_result_level == 1) {
      result_df <- pred
    } else if (report_result_level == 2) {
      result_df <- pred %>% group_by(sample_id, virus_name) %>%
        summarise(
          virus_name = first(virus_name),
          prediction = any(as.logical(prediction))
        )
    }
  }

  if (method %in% c("rf", "xgb") == F) {
    stop("`method` must be `rf` or `xgb`.")
  }

  result_df
}
