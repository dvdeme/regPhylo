context("Test_GeoCodeName")

data(Seq.Diastocapen)
Seq.DF1 = Seq.Diastocapen$Seq.DF1
# Run the function to retrieve geographic coordinates for sequences without coordinates.
Seq.DF3=GeoCodeName(input = Seq.DF1, output = "Seq.DF3.txt")

test_that("Test if there is two additional columns", {
  expect_equal(dim(Seq.DF3)[2], 28)
  if(length(which(colnames(Seq.DF3) == "Latitude_Y")) == 1) {a =1} else {a = 0}
  expect_equal(a, 1)
  if(length(which(colnames(Seq.DF3) == "Longitude_X")) == 1) {b =1} else {b = 0}
  expect_equal(b, 1)
  if(which(colnames(Seq.DF3) == "Geo_accuracy") == 25) {c =1} else {c = 0}
  expect_equal(c, 1)
  if(which(colnames(Seq.DF3) == "Location_used") == 21) {d =1} else {d = 0}
  expect_equal(d, 1)
})

test_that("Test if Inferred is present and Latitude and Longitude are numeric values", {
  if(length(which(Seq.DF3$Geo_accuracy == "Inferred")) == 1){a = 1} else {a = 0}
  expect_equal(a, 1)
  if(class(Seq.DF3$Latitude_Y[which(Seq.DF3$Geo_accuracy == "Inferred")]) == "numeric") {b=1} else {b=0}
  expect_equal(b, 1)
  if(class(Seq.DF3$Longitude_X[which(Seq.DF3$Geo_accuracy == "Inferred")]) == "numeric") {c=1} else {c=0}
  expect_equal(c, 1)
  expect_equal(Seq.DF3$Latitude_Y[which(Seq.DF3$Geo_accuracy == "Inferred")], -34.1573662)
  expect_equal(Seq.DF3$Longitude_X[which(Seq.DF3$Geo_accuracy == "Inferred")], 172.1349)
})

test_that("Test that the output table in the WD is identical to the output table in the R environment", {
  Seq.DF3.txt = read.delim("Seq.DF3.txt", sep = "\t", header = TRUE)
  expect_equal(dim(Seq.DF3.txt)[2], 28)
  if(length(which(colnames(Seq.DF3.txt) == "Latitude_Y")) == 1) {a =1} else {a = 0}
  expect_equal(a, 1)
  if(length(which(colnames(Seq.DF3.txt) == "Longitude_X")) == 1) {b =1} else {b = 0}
  expect_equal(b, 1)
  if(which(colnames(Seq.DF3.txt) == "Geo_accuracy") == 25) {c =1} else {c = 0}
  expect_equal(c, 1)
  if(which(colnames(Seq.DF3.txt) == "Location_used") == 21) {d =1} else {d = 0}
  expect_equal(d, 1)
  expect_equal(Seq.DF3.txt$Latitude_Y[which(Seq.DF3.txt$Geo_accuracy == "Inferred")], -34.1573662)
  expect_equal(Seq.DF3.txt$Longitude_X[which(Seq.DF3.txt$Geo_accuracy == "Inferred")], 172.1349)
})

file.remove("Seq.DF3.txt")
