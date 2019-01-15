context("Test_Mumsa.Comp")

input = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
test = Mumsa.Comp(input = input, output = "Mumsa_output",
remove.empty.align = FALSE)

test2 = read.delim("Mumsa_output.txt", sep = "\t", h=T)

test_that("test the output exported in the WD", {
  expect_equal(dim(test2)[1], 4)
  expect_equal(dim(test2)[2], 7)
  if(class(test[,2]) == "numeric"){a = 1} else {a = 0}
  expect_equal(a, 1)
})

test_that("test the output exported in the R envir.", {
  expect_equal(dim(test)[1], 4)
  expect_equal(dim(test)[2], 7)
  if(class(test[,2]) == "numeric"){a = 1} else {a = 0}
  expect_equal(a, 1)
})

# To clean the file created while running the example do the following:
file.remove("Mumsa_output.txt")
rm(test, test2)

# The function report an error message if some alignment are empty.
# Run the function without removing empty alignments
input = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
test = Mumsa.Comp(input = input, output = "Mumsa_output",
remove.empty.align = TRUE)

test2 = read.delim("Mumsa_output.txt", sep = "\t", h=T)

test_that("test the output exported in the WD", {
  expect_equal(dim(test2)[1], 4)
  expect_equal(dim(test2)[2], 7)
  if(class(test[,2]) == "numeric"){a = 1} else {a = 0}
  expect_equal(a, 1)
})

test_that("test the output exported in the R envir.", {
  expect_equal(dim(test)[1], 4)
  expect_equal(dim(test)[2], 7)
  if(class(test[,2]) == "numeric"){a = 1} else {a = 0}
  expect_equal(a, 1)
})


# To clean the file created while running the example do the following:
file.remove("Mumsa_output.txt")
