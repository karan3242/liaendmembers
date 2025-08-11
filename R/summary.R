#' Summary of LIA Endmember groups
#'
#' Summary of Lead Isotope analysis 'liaendmemebr' object
#'
#' @param object An object of class "liaendmembers"
#' @param ... further arguments passed to or from other methods.
#'
#' @returns
#' A list containing:
#'
#' \item{Counts}{The number of observations in each group}
#' \item{Tolerance}{Tolerance values used in the calculation.}
#' \item{Clamp}{Clamping values used in the calculation.}
#' \item{Data}{A data frame of the lead isotope ratios that were used.}
#'
#' @export
#' @inherit endmembers examples
summary.liaendmembers <- function(object, ...) {
     cat("Summary of End memebers:\n\n")
     cat("Tolarance:", object$tolarance, "\n")
     cat("Clamp:", object$clamp, "\n\n")
     cat("PCA Endmembers\n")
     print(object$pca_ends)
     cat("\n")
     count <- data.frame(
          "Group1" = nrow(object$group1),
          "Group2" = nrow(object$group2),
          "Miobjecting" = nrow(object$mixing),
          "Total" = sum(nrow(object$data))
     )
     row.names(count) <- "Counts"
     print(count)
     cat(rep("-", 18), "\n")
     print(summary(object$pca))
     invisible(list("Counts" = unlist(count),
                    "Tolarance"= object$tolarance,
                    "Clamp" = object$clamp,
                    "Data" = object$data))
}
