\name{qVarSel-package}
\alias{qVarSel-package}
\alias{qVarSel}
\docType{package}
\title{
Selecting Variables for Clustering and Classification
}
\description{
For a given data matrix A and cluster centers/prototypes collected in the matrix P, the functions described here select a subset of variables Q that mostly explains/justifies P as prototipes. The functions are useful to reduce the dimension of the data for classification and to discard masking variables for clustering.
}
\details{
\tabular{ll}{
Package: \tab qVarSel\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-05-27\cr
License: \tab gpl-2\cr
}
The package is useful to reduce the variable dimension for clustering. The example below shows the sequence of the operations. First, k-means can applied to the whole data sets, to calculate prototypes P. Then, distances between units U and P are calculated ans stored in a matrix D. Then, apply package subroutine q-VarSelH to select the most important variables. Apply EM optimization on data D for full clustering parameters estimation.
}
\author{
Stefano Benati 

Maintainer: Stefano Benati <stefano.benati@unitn.it>
}
\references{
S. Benati, S. Garcia Quiles, J. Puerto "Optimization Methods to Select Variables for Clustering", Working Paper, Universidad de Sevilla, 2014
}
\keyword{ classif }
\keyword{ cluster }
\keyword{ optimize }
\examples{

  # Generate random cluster with masking variables
  
  require(clusterGeneration)
  tmp1 <- genRandomClust(numClust = 4, sepVal = 0.2, clustszind = 2,
                         rangeN = c( 100, 150 ),
                         numNonNoisy = 5, numNoisy = 5, numReplicate = 1, 
                         fileName = "chk1")
  a <- tmp1$datList$chk1_1
  ass <- tmp1$memList$chk1_1
  numunits <- length(ass)
  noiseindex <- tmp1$noisyList$chk1_1
  a <- scale(a)  #Standardzation for columns
  
  # calculate data prototypes using k-means

  sl2 <- kmeans(a, 4, iter.max = 200, 
                      nstart = 10, algorithm = "L")
  prototype = sl2$centers
  
  # calculate distances between observations and prototypes
  # Remark: d is a 3-dimensions matrix
  
  d = PrtDist(a, prototype)
  
  # Select 5 most representative variables, use 200 iterations
  
  lsH <- qVarSelH(d, 5, maxit = 200)
  
  # reduce the dimension of a
  
  sq = 1:(dim(a)[2])
  vrb = sq[lsH$x > 0.01]
  a_reduced = a[ ,vrb]
  
  # use the EM methodology for clustering on the reduced data

  require(mclust)
  sl1 <- Mclust(a_reduced, G = 4, modelName = "VVV") 
}
