context("Test_Congr.NCBI.BOLD.perReposit")

test_that("Congr.NCBI.BOLD.perReposit assemble NCBI and BOLD data only using positive control", {
  data("Seq.Diastocapen")
  Seq.NCBI = Seq.Diastocapen$Seq.NCBI
  Seq.BOLD = Seq.Diastocapen$Seq.BOLD
  
  AllSeqDF=Congr.NCBI.BOLD.perReposit(input.NCBI=Seq.NCBI, input.BOLD=Seq.BOLD, output="AllSeq_NCBI_BOLD.txt")
  expect_equal(class(AllSeqDF), "data.frame")
  if(dim(AllSeqDF)[1]>0){a=1} else {a=0}
  expect_equal(a,1)
  OutPut1=read.delim("AllSeq_NCBI_BOLD.txt", sep="\t", h=T)
  expect_equal(dim(OutPut1)[1], dim(AllSeqDF)[1])
  expect_equal(dim(OutPut1)[2], dim(AllSeqDF)[2])
  expect_equal(dim(OutPut1)[2], 25) ### test if the column about the origin of the database is present
  expect_equal(length(unique(AllSeqDF[,25])), 3) ### test if the origines of the data base are properly
  # displayed especially when we use data information from the two Databases.
  file.remove("AllSeq_NCBI_BOLD.txt")
})

test_that("Congr.NCBI.BOLD.perReposit assemble NCBI, BOLD and the private repository using positive control", {
  data("Seq.Diastocapen")
  Seq.NCBI = Seq.Diastocapen$Seq.NCBI
  Seq.BOLD = Seq.Diastocapen$Seq.BOLD
  Seq.Perso = Seq.Diastocapen$Seq.Perso
  
  AllSeqDF2=Congr.NCBI.BOLD.perReposit(input.NCBI=Seq.NCBI, input.BOLD=Seq.BOLD,
                                       input.perReposit=Seq.Perso, perReposit="My.Rep",
                                       output="AllSeq_NCBI_BOLD_perRep.txt")
  expect_equal(class(AllSeqDF2), "data.frame")
  if(dim(AllSeqDF2)[1]>0){a=1} else {a=0}
  expect_equal(a,1)
  OutPut2=read.delim("AllSeq_NCBI_BOLD_perRep.txt", sep="\t", h=T)
  expect_equal(dim(OutPut2)[1], dim(AllSeqDF2)[1])
  expect_equal(dim(OutPut2)[2], dim(AllSeqDF2)[2])
  expect_equal(dim(OutPut2)[2], 25) ### test if the column about the origin of the database is present
  expect_equal(length(unique(AllSeqDF2[,25])), 4) ### test if the origines of the data base are properly
  # displayed especially when we use data information from the two or more databases.
  expect_equal(length(grep("My.Rep", AllSeqDF2[,25])), 1) ### 1 sequence is expecetd to come from My.Rep (the personal repository)
  file.remove("AllSeq_NCBI_BOLD_perRep.txt")
})

test_that("Congr.NCBI.BOLD.perReposit assemble NCBI and the private repository using positive control", {
  data("Seq.Diastocapen")
  Seq.NCBI = Seq.Diastocapen$Seq.NCBI
  Seq.Perso = Seq.Diastocapen$Seq.Perso

  AllSeqDF3=Congr.NCBI.BOLD.perReposit(input.NCBI=Seq.NCBI, input.perReposit=Seq.Perso,
                                       perReposit="My.Rep", output="AllSeq_NCBI_perRep.txt")
  expect_equal(class(AllSeqDF3), "data.frame")
  if(dim(AllSeqDF3)[1]>0){a=1} else {a=0}
  expect_equal(a,1)
  OutPut3=read.delim("AllSeq_NCBI_perRep.txt", sep="\t", h=T)
  expect_equal(dim(OutPut3)[1], dim(AllSeqDF3)[1])
  expect_equal(dim(OutPut3)[2], dim(AllSeqDF3)[2])
  expect_equal(dim(OutPut3)[2], 25) ### test if the column about the origin of the database is present
  expect_equal(length(unique(AllSeqDF3[,25])), 2) ### test if the origines of the data base are properly
  # displayed especially when we use data information from the two or more databases.
  expect_equal(length(grep("My.Rep", AllSeqDF3[,25])), 1) ### 1 sequence is expecetd to come from My.Rep (the personal repository)
  expect_equal(length(grep("NCBI", AllSeqDF3[,25])), 4) ### four sequences are expecetd to come from NCBI database
  file.remove("AllSeq_NCBI_perRep.txt")
})

test_that("Congr.NCBI.BOLD.perReposit assemble BOLD and the private repository using positive control", {
  data("Seq.Diastocapen")
  Seq.BOLD = Seq.Diastocapen$Seq.BOLD
  Seq.Perso = Seq.Diastocapen$Seq.Perso

  AllSeqDF4=Congr.NCBI.BOLD.perReposit(input.BOLD=Seq.BOLD, input.perReposit=Seq.Perso,
                                       perReposit="My.Rep", output="AllSeq_BOLD_perRep.txt")
  expect_equal(class(AllSeqDF4), "data.frame")
  if(dim(AllSeqDF4)[1]>0){a=1} else {a=0}
  expect_equal(a,1)
  OutPut4=read.delim("AllSeq_BOLD_perRep.txt", sep="\t", h=T)
  expect_equal(dim(OutPut4)[1], dim(AllSeqDF4)[1])
  expect_equal(dim(OutPut4)[2], dim(AllSeqDF4)[2])
  expect_equal(dim(OutPut4)[2], 25) ### test if the column about the origin of the database is present
  expect_equal(length(unique(AllSeqDF4[,25])), 2) ### test if the origines of the data base are properly
  # displayed especially when we use data information from the two or more databases.
  expect_equal(length(grep("My.Rep", AllSeqDF4[,25])), 1) ### 1 sequence is expecetd to come from My.Rep (the personal repository)
  expect_equal(length(grep("BOLD", AllSeqDF4[,25])), 8) ### 8 sequences are expected to come from BOLD.
  file.remove("AllSeq_BOLD_perRep.txt")
})

