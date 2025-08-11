test_that("summary workes", {
     expect_no_error(summary(endmembers(akko, names(akko))))
     expect_length(summary(endmembers(akko, names(akko))), 4)
})
