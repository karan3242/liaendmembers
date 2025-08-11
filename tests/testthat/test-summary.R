test_that("summary workes", {
     su <- summary(endmembers(akko, names(akko)))
     expect_no_error(su)
     expect_length(su, 4)
})
