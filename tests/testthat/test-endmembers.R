     test_that("endmembers() validation works as expected", {
     # Test for errors when required arguments are missing or invalid
     expect_error(endmembers(1:10))
     expect_error(endmembers(dor), "column names needed!")
     expect_error(endmembers(dor, col = c("pb64", "pb74")))


})

---

     test_that("endmembers() handles warnings and messages", {
          # This warning is expected when the threshold is too high
          expect_warning(endmembers(dor, names(dor), 0.5), "Overlap in endmembers between gorups. Suggest lower tolerance value.")

          # This message is expected with normality issues
          expect_message(endmembers(dor, names(dor)), "PC2 or PC3 are not normally distributed. This may indicate that their variation may not be random noise.")
     })

---

     test_that("endmembers() returns the correct class", {
          # Assign the result of the function call to a variable to avoid re-running it
          result <- endmembers(dor, names(dor))

          # Check if the returned object has the correct S3 classes
          expect_s3_class(result, "liaendmembers")
          expect_s3_class(result, "list")
     })
