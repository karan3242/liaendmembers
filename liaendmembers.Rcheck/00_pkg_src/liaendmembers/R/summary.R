#' Summary of LIA Endmember groups
#'
#' Summary of Lead Isotope analysis 'liaendmemebr' object
#'
#' @param x An object of class "liaendmembers"
#' @param ... further arguments passed to or from other methods.
#'
#' @returns
#' A list containing:
#'
#' \item{Counts}{The number of observations in each group}
#' \item{Tolerance}{The tolerance value used in the calculation.}
#' \item{Data}{A data frame of the lead isotope ratios that were used.}
#'
#' @export
#' @inherit endmembers examples
summary.liaendmembers <- function(x, ...) {
     cat("Summary of End memebers:\n\n")
     cat("Tolarance:", x$tolarance, "\n\n")
     count <- data.frame(
          "Group1" = nrow(x$group1),
          "Group2" = nrow(x$group2),
          "Mixing" = nrow(x$mixing),
          "Total" = sum(nrow(x$group1), nrow(x$group2), nrow(x$mixing))
     )
     row.names(count) <- "Counts"
     print(count)
     cat(rep("-", 18), "\n")
     print(summary(x$pca))
     invisible(list("Counts" = unlist(count),
                    "Tolarance"= x$tolarance,
                    "Data" = rbind(x$group1, x$group2, x$mixing)))
}
