
flowDir <- paste0(system.file(package = "PloidyPeaks"), "/gated_data/")
subsetDs <- c("A01-A01", "A04-D12", "A07-G12", "A09-A02", "T1-D08")
flowSet <- subset(
  list.files(flowDir),
  list.files(flowDir) %in% subsetDs
)
improperDataFormat <- c(
  "csv", "xls", "xlsx","html", "ppt",
  "pptx" ,"pdf", "doc", "docx"
)

test_that("Data is in proper flow format", {
  for(k in seq_len(length(flowSet))){
    expect_false(tools::file_ext(flowSet[k]) %in% improperDataFormat) 
  }
  
})

#Check xVariable is in flow sample
xVariable = "FITC-A"
test_that("xVariable is in flow sample", {
  for(k in seq_len(length(flowSet))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowSet[k]), transformation=FALSE
    )
    expect_false(!xVariable %in% flowName@parameters@data$name)
  }
})

#Complicated calculations in FlowPeakDetection()
test_that("Formula calculation is correct (1)", {
  x <- 5
  N1 <- 100
  g1SD <- 2
  g1Mean <- 4
  N2 <- 120
  g2SD <- 1
  g2Mean <- 10
  G1g2SD <- 2
  G1G2Mean <- 5
  G2g2SD <- 1.5
  G2G2Mean <- 15
  numDoublet1 <- 10
  numDoublet2 <- 8
  A1 <- 1
  B1 <- 0.5
  C1 <- -0.2
  A2 <- -1
  B2 <- 0.4
  C2 <- 0.3
  A3 <- 0.8
  B3 <- -0.6
  C3 <- 0.1
  expectedResult <- 19.45
  computedResult <- (N1/(sqrt(2*pi)*g1SD)*exp(-((x-g1Mean)^2)/(2 *g1SD^2))) +
    (A1 + B1*x + C1*(x^2))*(1/(sqrt(2 * pi)*g1SD*(x/g1Mean))*exp(-((x - g1Mean)^2)/(2 *(g1SD* (x/g1Mean))^2)))+
    (N2/(sqrt(2 * pi) * g2SD) * exp(-((x - g2Mean)^2)/(2 *g2SD^2)))+
    (A2 + B2*x + C2*(x^2))*(1/(sqrt(2 * pi)*g2SD*(x/g2Mean))*exp(-((x - g2Mean)^2)/(2 *(g2SD* (x/g2Mean))^2)))+
    (numDoublet1/(sqrt(2*pi)*G1g2SD)*exp(-((x-G1G2Mean)^2)/(2*G1g2SD^2)))+
    (A3 + B3*x + C3*(x^2))*(1/(sqrt(2*pi)*G1g2SD*(x/G1G2Mean))*exp(-((x-G1G2Mean)^2)/(2*(G1g2SD*(x/G1G2Mean))^2)))+
    (numDoublet2/(sqrt(2*pi)*G2g2SD)*exp(-((x-G2G2Mean)^2)/(2 *G2g2SD^2)))
  expect_equal(computedResult, expectedResult, tolerance = 1e-2)
})

test_that("Formula calculation is correct (2)", {
  x <- 5
  g1SD <- 1
  g1Mean <- 5
  g2SD <- 1.5
  g2Mean <- 10
  g1SD2 <-  2
  g1Mean2 <- 7
  g2SD2 <- 2.5
  g2Mean2 <- 12
  g2SD3 <- 3
  g2Mean3 <- 15
  A <- 1
  B <- 0.5
  C <- -0.2
  A2 <- -1
  B2 <- 0.3
  C2 <- 0.4
  A3 <- 0.8
  B3 <- -0.6
  C3 <- 0.2
  A4 <- 1.2
  B4 <- -0.8
  C4 <- 0.3
  pop1N1 <- 100
  pop1N2 <- 120
  pop2N1 <- 80
  pop2N2 <- 90
  pop3N2 <- 110
  expectedResult <- 49.45
  computedResult <- (pop1N1/(sqrt(2*pi)*g1SD)*exp(-((x-g1Mean)^2)/(2*g1SD^2)))+
    (pop1N2/(sqrt(2*pi)*g2SD)* exp(-((x-g2Mean)^2)/(2*g2SD^2)))+
    (pop2N1/(sqrt(2*pi)*g1SD2)*exp(-((x-g1Mean2)^2)/(2*g1SD2^2)))+
    (pop2N2/(sqrt(2*pi)*g2SD2)* exp(-((x-g2Mean2)^2)/(2*g2SD2^2)))+
    (pop3N2/(sqrt(2*pi)*g2SD3)* exp(-((x-g2Mean3)^2)/(2*g2SD3^2)))+
    (A + B*x + C*(x^2))*(1/(sqrt(2*pi)*g1SD*(x/g1Mean))*exp(-((x-g1Mean)^2)/(2*(g1SD*(x/g1Mean))^2)))+
    (A2 + B2*x + C2*(x^2))*(1/(sqrt(2*pi)*g2SD2*(x/g2Mean2))*exp(-((x-g2Mean2)^2)/(2*(g2SD2*(x/g2Mean2))^2)))+
    (A3 + B3*x + C3*(x^2))*(1/(sqrt(2*pi)*g2SD3*(x/g2Mean3))*exp(-((x-g2Mean3)^2)/(2*(g2SD3*(x/g2Mean3))^2)))+
    (A4 + B4*x + C4*(x^2))*(1/(sqrt(2*pi)*g2SD*(x/g2Mean))*exp(-((x-g2Mean)^2)/(2*(g2SD*(x/g2Mean))^2)))
  expect_equal(computedResult, expectedResult, tolerance = 1e-2)
})
