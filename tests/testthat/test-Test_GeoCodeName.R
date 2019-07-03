context("Test_GeoCodeName")


data(Seq.Diastocapen)
Seq.DF1 = Seq.Diastocapen$Seq.DF1
# Run the function to retrieve geographic coordinates for sequences without coordinates, using AutoCorrection of place name for NZ.
Seq.DF3=GeoCodeName(input = Seq.DF1, output = "Seq.DF3.txt", AutoCorrNZ = TRUE)

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

rm(Seq.DF3)
file.remove("Seq.DF3.txt")


# Run the function to retrieve geographic coordinates for sequences without coordinates, using AutoCorrection of place name for NZ.
Seq.DF3b=GeoCodeName(input = Seq.DF1, output = "Seq.DF3.txt", AutoCorrNZ = FALSE)

test_that("Without AutoCorrNZ set to true the automatic correction of the place name related to NZ study are not performed", {
  if(which(Seq.DF3b[,25] == "NoLocationFound") == 3){e = 1} else {e = 0}
  expect_equal(e, 1)
  if(length(which(is.na(Seq.DF3b[,23]) == TRUE)) == 3){f = 1} else {f = 0}
  expect_equal(f, 1)
})

rm(Seq.DF3b)
file.remove("Seq.DF3.txt")


# Run the function to retrieve geographic coordinates for sequences without coordinates, without the AutoCorrNZ,
# but including the correction table.
correctionTab = matrix(c("Tasman Sea: , , NA, South Norfolk Ridge, off the Three Kings Islands",
                         "New Zealand: Three Kings Islands"), ncol=2)
colnames(correctionTab) = c("OriginalLocationName", "CorrectedLocationName")

Seq.DF3=GeoCodeName(input = Seq.DF1, output = "Seq.DF3.txt", AutoCorrNZ = FALSE, CorrTab = correctionTab)

test_that("Without AutoCorrNZ but with the correction table", {
  if(length(which(Seq.DF3$Geo_accuracy == "Inferred")) == 1){a = 1} else {a = 0}
  expect_equal(a, 1)
  if(class(Seq.DF3$Latitude_Y[which(Seq.DF3$Geo_accuracy == "Inferred")]) == "numeric") {b=1} else {b=0}
  expect_equal(b, 1)
  if(class(Seq.DF3$Longitude_X[which(Seq.DF3$Geo_accuracy == "Inferred")]) == "numeric") {c=1} else {c=0}
  expect_equal(c, 1)
  expect_equal(Seq.DF3$Latitude_Y[which(Seq.DF3$Geo_accuracy == "Inferred")], -34.1573662)
  expect_equal(Seq.DF3$Longitude_X[which(Seq.DF3$Geo_accuracy == "Inferred")], 172.1349)
})

rm(Seq.DF3)
file.remove("Seq.DF3.txt")

# Run the function with a table with all the sequences with geographic coordinates already.
Seq.DF1 = Seq.DF1[-c(1:3),]
Seq.DF3 = GeoCodeName(input = Seq.DF1, output = "Seq.DF3.txt", AutoCorrNZ = TRUE)

test_that("Test if there is two additional columns, When all sequences have geographic coordinates", {
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

test_that("Test if when all sequences have already some geographic coordinates the output is has expected", {
  expect_equal(dim(Seq.DF3)[1], 8)
  expect_equal(length(which(Seq.DF3[,25] == "From_DB")), 8)
  expect_warning(GeoCodeName(input = Seq.DF1, output = "Seq.DF3.txt", AutoCorrNZ = TRUE),
                 "All the sequences")
})

rm(Seq.DF3, Seq.DF1)
file.remove("Seq.DF3.txt")

# Run the function with a table with no place name for the sequences without geographic coordinates.
Seq.DF1 = Seq.Diastocapen$Seq.DF1
Seq.DF1[3,20] = NA
Seq.DF3 = GeoCodeName(input = Seq.DF1, output = "Seq.DF3.txt", AutoCorrNZ = TRUE)

test_that("Test if when all seqeunces without geographic coordinates neither place name the output is correct", {
  expect_equal(dim(Seq.DF3)[1], 11)
  expect_equal(dim(Seq.DF3)[2], 28)
  expect_warning(GeoCodeName(input = Seq.DF1, output = "Seq.DF3.txt", AutoCorrNZ = TRUE), "None of the DNA sequences")
})

rm(Seq.DF3)
file.remove("Seq.DF3.txt")
