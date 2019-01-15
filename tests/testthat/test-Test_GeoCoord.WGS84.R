context("test-Test_GeoCoord.WGS84")

test_that("Test GeoCoord.WGS84 split ", {
  data("Seq.Diastocapen")
  Seq.DF = Seq.Diastocapen$Seq.DF
  # Run the function
  Seq.DF1=GeoCoord.WGS84(input=Seq.DF, output="Seq.DF.txt")
  expect_equal(dim(Seq.DF1)[2], 26) # test there is one additional column
  # Test if the name of the new columns are Lat and Long
  if(length(which(colnames(Seq.DF1) == "Lat")) == 1) {a =1} else {a = 0}
  expect_equal(a, 1)
  if(length(which(colnames(Seq.DF1) == "Lat")) == 1) {b =1} else {b = 0}
  expect_equal(b, 1)
})

test_that("Test GeoCoord.WGS84 Lat and Long field are numeric", {
  data("Seq.Diastocapen")
  Seq.DF = Seq.Diastocapen$Seq.DF
  # Run the function
  Seq.DF1=GeoCoord.WGS84(input=Seq.DF, output="Seq.DF.txt")
  # Test if the colum Lat is numeric
  if(class(Seq.DF1$Lat)=="numeric") {a = 1} else {a = 0}
  expect_equal(a,1)
  if(class(Seq.DF1$Long)=="numeric") {b = 1} else {b = 0}
  expect_equal(b,1)
})

test_that("Test GeoCoord.WGS84 the mean Lat and mean long are numeric", {
  data("Seq.Diastocapen")
  Seq.DF = Seq.Diastocapen$Seq.DF
  # Run the function
  Seq.DF1=GeoCoord.WGS84(input=Seq.DF, output="Seq.DF.txt")
  # Test if the colum Lat is numeric
  if(class(mean(Seq.DF1$Lat, na.rm = TRUE)) == "numeric"){a = 1} else {a = 0}
  expect_equal(a,1)
  # Test if the colum Long is numeric
  if(class(mean(Seq.DF1$Long, na.rm = TRUE)) == "numeric"){b = 1} else {b = 0}
  expect_equal(b,1)
})

test_that("Test GeoCoord.WGS84 output is identical to the table in R envir", {
  data("Seq.Diastocapen")
  Seq.DF = Seq.Diastocapen$Seq.DF
  # Run the function
  Seq.DF1 = GeoCoord.WGS84(input=Seq.DF, output="Seq.DF.txt")
  Seq.DF.txt = read.delim("Seq.DF.txt", sep = "\t", header = TRUE)
  # test if the mean of the Lat and Long are equal between Seq.DF.txt and Seq.DF1
  if(mean(Seq.DF1$Lat, na.rm = TRUE) == mean(Seq.DF.txt$Lat, na.rm = TRUE)){a = 1} else {a = 0}
  expect_equal(a,1)
  # Test if the colum Long is numeric
  if(mean(Seq.DF1$Long, na.rm = TRUE) == mean(Seq.DF.txt$Long, na.rm = TRUE)){b = 1} else {b = 0}
  expect_equal(b,1)
})

file.remove("Seq.DF.txt")

