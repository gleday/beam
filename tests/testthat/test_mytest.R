context("mycontext")

data("TCPAprad")

test_that("output", {

  ##################################################
  ## beam estimation
  fit <- beam(X = TCPAprad, type="both", verbose=FALSE)

  ## Test beam-class object
  # Test table
  expect_equal(dim(fit@table), c(17766, 6))
  expect_equal(colnames(fit@table), c("m_cor", "m_logBF", "m_tail_prob", "p_cor", "p_logBF", "p_tail_prob"))
  expect_equal(fit@table[1,1], -0.3499239, tolerance=1e-5)
  expect_equal(fit@table[1,2], 10.04931, tolerance=1e-5)
  expect_equal(fit@table[1,3], 3.897943e-08, tolerance=1e-5)
  expect_equal(fit@table[1,4], 0.00415224, tolerance=1e-5)
  expect_equal(fit@table[1,5], -0.8375797, tolerance=1e-5)
  expect_equal(fit@table[1,6], 0.8700933, tolerance=1e-5)
    
  # Test deltaOpt
  expect_equal(length(fit@deltaOpt), 1)
  
  # Test alphaOpt
  expect_equal(length(fit@alphaOpt), 1)
  
  # Test dimX
  expect_equal(fit@dimX, c(164, 189))
  
  # Test type
  expect_equal(fit@type, "both")
  
  # Test varlabs
  expect_equal(length(fit@varlabs), 189)
  
  # Test gridAlpha
  expect_equal(dim(fit@gridAlpha), c(100, 3))
  
  # Test valOpt
  expect_equal(length(fit@valOpt), 1)
  
  # Test return.only
  expect_equal(fit@return.only, c("cor","BF","prob"))
  
  
  
  ## Test accessors beam-class object
  # mcor
  expect_equal(dim(mcor(fit)), c(189, 189))
  
  # mcor
  expect_equal(dim(pcor(fit)), c(189, 189))
  
  # bgraph
  expect_equal(class(bgraph(fit)), "igraph")
  
  # ugraph
  expect_equal(class(ugraph(fit)), "igraph")
  
  # summary
  expect_output(summary(fit))
  
  # print
  expect_output(print(fit))
  
  # show
  expect_output(show(fit))
  
  ##################################################
  ## beam selection
  sel <- beam.select(fit, thres=0.01, method = "BH") 
  
  ## Test beam.select-class object
  # Test marginal
  expect_equal(dim(sel@marginal), c(6101, 4))
  expect_equal(colnames(sel@marginal), c("m_cor", "m_logBF", "m_tail_prob", "m_tail_prob_BH"))
  
  # Test conditional
  expect_equal(dim(sel@conditional), c(154, 4))
  expect_equal(colnames(sel@conditional), c("p_cor", "p_logBF", "p_tail_prob", "p_tail_prob_BH"))
  
  # Test method
  expect_equal(sel@method, "BH")
  
  # Test thres
  expect_equal(sel@thres, 0.01)
  
  ## Test accessors beam.select-class object
  # Test marg
  expect_equal(is.data.frame(marg(sel)), TRUE)
  
  # Test cond
  expect_equal(is.data.frame(cond(sel)), TRUE)
  
  # Test mcor
  expect_equal(dim(mcor(sel)), c(189, 189))
  
  # Test pcor
  expect_equal(dim(pcor(sel)), c(189, 189))
  
  # summary
  expect_output(summary(sel))
  
  # print
  expect_output(print(sel))
  
  # show
  expect_output(show(sel))
  
  ## beam selection
  sel <- beam.select(fit, thres=0.01, method = "HC") 
  
  # Test method
  expect_equal(sel@method, "HC")
  
  # Test marg
  expect_equal(is.data.frame(marg(sel)), TRUE)
  
  # Test cond
  expect_equal(is.data.frame(cond(sel)), TRUE)
  
  # Test mcor
  expect_equal(dim(mcor(sel)), c(189, 189))
  
  # Test pcor
  expect_equal(dim(pcor(sel)), c(189, 189))
  
})