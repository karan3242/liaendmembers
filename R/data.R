#' LIA points of Silver Hoards from Israel
#'
#' Lead Isotope ratios of a Phoenician Silver Hoard from Tel Dor, Akko, Arad, Ein Hofez and  Eshtemona.
#' Reference Data set of Isotope rations and a pretrained model using default parameters from [train_data()]
#' @name data
#' @format **For artefact sampels:** Data frame with 3 variables
#' \describe{
#'   \item{pb64}{Lead Isotope ratio of 206Pb/204Pb}
#'   \item{pb74}{Lead Isotope ratio of 207Pb/204Pb}
#'   \item{pb84}{Lead Isotope ratio of 208Pb/204Pb}
#' }
#' @format **For reference data:** Object of class `ref.data` with 4 variables.
#' \describe{
#'   \item{group}{Regions of Ore samples}
#'   \item{pb64}{Lead Isotope ratio of 206Pb/204Pb}
#'   \item{pb74}{Lead Isotope ratio of 207Pb/204Pb}
#'   \item{pb84}{Lead Isotope ratio of 208Pb/204Pb}
#' }
#' @format **For Xgboost model:** List with xgb.Booster objects tranied usen parameters following (Shnyr et al., 2026)
#'
#' @source Eshel, T., Erel, Y., Yahalom-Mack, N., Tirosh, O., & Gilboa, A.
#'   (2019). Lead isotopes in silver reveal earliest Phoenician quest for metals
#'   in the west Mediterranean. Proceedings of the National Academy of Sciences,
#'   116(13), 6007–6012.
#'   \href{https://doi.org/10.1073/pnas.1817951116}{10.1073/pnas.1817951116}
#'
#' @inherit endmembers references examples
"tel_dor"

#' @rdname data
"akko"

#' @rdname data
"arad"

#' @rdname data
"ein_hofez"

#' @rdname data
"eshtemoa"

#' @rdname data
"ref_data"

#' @rdname data
"model"
