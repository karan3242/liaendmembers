summary.liaendmembers <- function(x) {
     cat("Summary of End memebers:\n\n")
     cat("Tolarance:", x$tolarance, "\n\n")
     count <- data.frame(
          "Group1" = nrow(x$group1),
          "Group2" = nrow(x$group2),
          "Mixing" = nrow(x$mixing)
     )
     row.names(count) <- "Counts"
     print(count)
     cat(rep("-", 14), "\n")
     print(summary(x$pca))
     invisible(list("Counts" = unlist(count),
                    "Tolarance"= x$tolarance,
                    "Data" = rbind(x$group1, x$group2, x$mixing)))
}
