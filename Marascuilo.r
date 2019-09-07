# As you probably know, the Marascuilo procedure is used to analyze the difference between two proportions in a contingency table to determine if that difference in proportion is significant or not.

marascuilo = function(dataFrame,confidence=.95) {
  
  chiResult = chisq.test (dataFrame, correct=FALSE )
  xSquared = chiResult$statistic
  
  # Generate all possible pair-wise combinations of groups
  colNames = length(dataFrame)/2
  combos = combn(colNames , 2)
  numCombos = dim(combos)[2]  # combos is an array of pairs, we want the length
  
  
  # Allocate matrix (initially 0 rows) for results
  results = matrix(nrow=0, ncol=5, dimnames=getResultsColumNames() )
  
  chiSquaredConstant = calcChiSquaredConstant(dataFrame, confidence)
  for (i in 1: numCombos) { 
    newRow = testSignificanceOfAbsDiffVsCriticalRange(
      dataFrame, combos, i, chiSquaredConstant ) 
    results = rbind(results, newRow)        # append new row to results
  }
  
  
  # sort results so that the pair differences that most strikingly exceed 
  # the critical range appear toward the top.
  sortedResults = results[  order( results[,'abs.diff-critical.range'] ) , ]
  return (sortedResults )
}


calcChiSquaredConstant = function(dataFrame,confidence) {
  nRows = dim(dataFrame)[1]  
  nCols = dim(dataFrame)[2]  
  
  degreesFreedom =  (nRows-1) * (nCols-1) 
  chiSquaredConstant = sqrt( qchisq(confidence,degreesFreedom) )
  
  return (chiSquaredConstant)
}


getResultsColumNames =  function (numRows) {
  return ( list( c(), c('pair', 'abs.diff', 'critical.range', 'abs.diff-critical.range', 'significant')) )
}

# test significance for ith combination
#
testSignificanceOfAbsDiffVsCriticalRange = function(
  dataFrame, combos, i,  chiSquaredConstant) {
  
  results = matrix(nrow=1, ncol=5, dimnames=getResultsColumNames() )
  
  pair1=combos[1,i]
  pair2=combos[2,i]
  
  # sum column denoted by name 'pair1' into groupTotal1 
  groupTotal1 = sum( dataFrame[ , pair1])  
  groupTotal2 = sum( dataFrame[ , pair2])  # do same thing for pair2... 
  
  p1 = dataFrame[1, pair1] / groupTotal1 
  p2 = dataFrame[1, pair2] / groupTotal2
  p1Not = (1 - p1)
  p2Not = (1 - p2)
  
  absDiff = abs( p2  - p1 )
  
  criticalRange = chiSquaredConstant  * 
    sqrt(p1*p1Not/groupTotal1 + p2*p2Not/groupTotal2)
  results[1, 'pair'] = paste(pair1,"|",pair2) 
  results[1, 'abs.diff'] = round(absDiff,3)
  results[1, 'critical.range'] = round(criticalRange ,3)
  results[1, 'abs.diff-critical.range'] = round(absDiff - criticalRange ,3)
  
  
  if (absDiff > criticalRange) {
    results[1, 'significant'] = 'Y'
  } else {
    results[1, 'significant'] = 'N'
  } 
  
  return(results)
}
M <- t(as.matrix(rbind(c(213, 22), c(193, 42),c(189, 46),c(205, 30),
                       c(195, 40),c(165, 57),c(154, 80))))
dimnames(M) <- list(lethality = c("D", "N"),
                    algorithms = c("PROVEAN","SIFT", "I-Mutant", "PolyPhen", 
                              "PHD-SNP", "PANTHER", "SNP&GO"))
chisq.test(M)
marascuilo(M)


