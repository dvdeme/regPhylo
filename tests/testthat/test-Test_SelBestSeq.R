context("Test_SelBestSeq")

Seq.DF5 = Seq.DF4$Seq.DF5

Seq.DF6 = SelBestSeq(input = Seq.DF5, output = "Alig_Seq.DF5.All",
                   RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
                   Alignment = T, MaxSeq = "ALL", gene.list = c("co1", "16srrna"),
                   SeqChoice = "Median")
dim(Seq.DF6)


test_that("Export all Sequences check output R envir", {
  expect_equal(dim(Seq.DF6)[1], 10)
  expect_equal(dim(Seq.DF6)[2], 33)
  # test that the NA for the distance to centroid of new zealand is in the same position (row)
  # that the NA in the longitude_X coordinates.
  expect_equal(which(is.na(Seq.DF6$Longitude_X)), which(is.na(Seq.DF6$Dist2RefPoint)))
  if(class(Seq.DF6$SeqLengthNoIndelsNoAmbig) == "integer") {a = 1} else {a = 0}
  expect_equal(a, 1)
  if(class(Seq.DF6$dist2med) == "integer") {b = 1} else {b = 0}
  expect_equal(b, 1)
})


Seq.DF6b= read.delim("Alig_Seq.DF5.All_InfoTabSelSeq.txt", sep = "\t", header = TRUE)

test_that("Export all Sequences check output in the WD", {
  expect_equal(dim(Seq.DF6b)[1], 10)
  expect_equal(dim(Seq.DF6b)[2], 33)
  # test that the NA for the distance to centroid of new zealand is in the same position (row)
  # that the NA in the longitude_X coordinates.
  expect_equal(which(is.na(Seq.DF6b$Longitude_X)), which(is.na(Seq.DF6b$Dist2RefPoint)))
  if(class(Seq.DF6b$SeqLengthNoIndelsNoAmbig) == "integer") {a = 1} else {a = 0}
  expect_equal(a, 1)
  if(class(Seq.DF6b$dist2med) == "integer") {b = 1} else {b = 0}
  expect_equal(b, 1)
})


ali16s = ape::read.dna("Alig_Seq.DF5.All_16srrna.fas", format = "fasta")
alico1 = ape::read.dna("Alig_Seq.DF5.All_co1.fas", format = "fasta")


test_that("Export all Sequences check seqeunces alignments", {
  if(class(alico1) == "DNAbin") {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(labels(ali16s)), 1)
  expect_equal(length(ali16s), 644)
  expect_equal(length(labels(alico1)), 9)
  expect_equal(length(alico1[[1]]), 654)
})

### selecting the best sequence and longest
Seq.DF6c = SelBestSeq(input = Seq.DF5, output = "Alig_Seq.DF5.Best",
                     RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
                     Alignment = T, MaxSeq = 1, gene.list = c("co1", "16srrna"),
                     SeqChoice = "Longest")


ali16sbest = ape::read.dna("Alig_Seq.DF5.Best_16srrna.fas", format = "fasta")
alico1best = ape::read.dna("Alig_Seq.DF5.Best_co1.fas", format = "fasta")


test_that("Export best sequence check sequences alignments", {
  if(class(alico1best) == "DNAbin") {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(labels(ali16sbest)), 1)
  expect_equal(length(ali16sbest), 644)
  expect_equal(length(labels(alico1best)), 1)
  expect_equal(length(alico1best), 654)
})

test_that("Export all Sequences check output R envir", {
  expect_equal(dim(Seq.DF6c)[1], 2)
  expect_equal(dim(Seq.DF6c)[2], 33)
  # test that the NA for the distance to centroid of new zealand is in the same position (row)
  # that the NA in the longitude_X coordinates.
  expect_equal(which(is.na(Seq.DF6c$Longitude_X)), which(is.na(Seq.DF6c$Dist2RefPoint)))

  if(class(Seq.DF6c$SeqLengthNoIndelsNoAmbig) == "integer") {a = 1} else {a = 0}
  expect_equal(a, 1)
  if(class(Seq.DF6c$dist2med) == "integer") {b = 1} else {b = 0}
  expect_equal(b, 1)
})



### selecting the best sequence and longest
Seq.DF6c = SelBestSeq(input = Seq.DF5, output = "Alig_Seq.DF5.Best",
                      RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
                      Alignment = T, MaxSeq = 1, gene.list = c("co1", "16srrna"),
                      SeqChoice = "Median")
ali16sbest = ape::read.dna("Alig_Seq.DF5.Best_16srrna.fas", format = "fasta")
alico1best = ape::read.dna("Alig_Seq.DF5.Best_co1.fas", format = "fasta")


test_that("Export best sequence (median) check sequences alignments", {
  if(class(alico1best) == "DNAbin") {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(labels(alico1best)), 1)
  expect_equal(length(alico1best), 654)
})


### Remove the geographic coordinates
colnames(Seq.DF5)
Seq.DF5b = cbind(Seq.DF5[, c(1:22)], Latitude_Y = rep(NA, dim(Seq.DF5)[1]), Longitude_X = rep(NA, dim(Seq.DF5)[1]),
                 Seq.DF5[, c(25:29)])

Seq.DF6c = SelBestSeq(input = Seq.DF5b, output = "Alig_Seq.DF5.Best",
                      RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
                      Alignment = T, MaxSeq = 1, gene.list = c("co1", "16srrna"),
                      SeqChoice = "Median")

ali16sbest = ape::read.dna("Alig_Seq.DF5.Best_16srrna.fas", format = "fasta")
alico1best = ape::read.dna("Alig_Seq.DF5.Best_co1.fas", format = "fasta")

test_that("Export best sequence (median) check sequences alignments, no geographic coordinates test the length", {
  if(class(alico1best) == "DNAbin") {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(labels(alico1best)), 1)
  expect_equal(length(alico1best), 652)
})

test_that("Test the R output in R envir, No geographic coordinates, median seq length", {
  expect_equal(sum(Seq.DF6c$dist2med), 0) # sum of the distance to the median length
  # of the distribution of the sequence length should be 0.
})

Seq.DF6c = SelBestSeq(input = Seq.DF5b, output = "Alig_Seq.DF5.Best",
                      RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
                      Alignment = T, MaxSeq = 1, gene.list = c("co1", "16srrna"),
                      SeqChoice = "Longest")

ali16sbest = ape::read.dna("Alig_Seq.DF5.Best_16srrna.fas", format = "fasta")
alico1best = ape::read.dna("Alig_Seq.DF5.Best_co1.fas", format = "fasta")

test_that("Export best sequence (Longest) check sequences alignments, no geographic coordinates test the length", {
  if(class(alico1best) == "DNAbin") {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(labels(alico1best)), 1)
  expect_equal(length(alico1best), 655)
})

Seq.DF6c$Dist2RefPoint

test_that("Test the R output in R envir, No geographic coordinates", {
  if(length(which(is.na(Seq.DF6c$Dist2RefPoint) == TRUE)) == 2){a = 1} else {a = 0}
  expect_equal(a, 1)
})


rm(Seq.DF6, Seq.DF6b, Seq.DF6c, ali16sbest, ali16s, alico1, alico1best)

file.remove("Alig_Seq.DF5.All_16srrna.fas", "Alig_Seq.DF5.Best_16srrna.fas",
            "Alig_Seq.DF5.Best_co1.fas", "Alig_Seq.DF5.All_co1.fas",
            "Alig_Seq.DF5.All_InfoTabSelSeq.txt", "Alig_Seq.DF5.Best_InfoTabSelSeq.txt")

