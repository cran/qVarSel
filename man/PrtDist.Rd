\name{PrtDist}
\alias{PrtDist}
\title{
Calculation of the distances between units and centers
}
\description{
Given a data set A, with a[i,k] be the measure of variable k on unit i, and a prototype set P, with p[j,k] be the measure of variable k on prototype j,
the function calculate d[i,j,k], the i,j distance according variable k 
}
\usage{
PrtDist(a, 
        p)
}
\arguments{
  \item{a}{ Matrix with n rows (units) and m columns (variables)} 
  \item{p}{ Matrix with g rows (units) and m columns (variables)}
}

\details{d[i,j,k] is the squared distance: d[i,j,k] = (a[i,k] - p[j,k])^2
}
\value{d: The 3-dimensional matrix of i,j,k distances}

\references{
S. Benati, S. Garcia Quiles, J. Puerto "Optimization Methods to Select Variables for Clustering", Working Paper, Universidad de Sevilla, 2014
}
\author{
Stefano Benati
}
\note{
The function has been written to simplify the code examples. It has been written as a R script with minimal vectorization, therefore it is not computationally much efficient. Moreover d(i,j,k) is the square of differences and users may prefer to employ other dissimilarity measures. In that case, they should better write their own function.   
}
\examples{
## Generate random 2 cluster with 20 masking variables
## and 10 true variables
  
require(clusterGeneration)
tmp1 <- genRandomClust(numClust = 2, sepVal = 0.2, clustszind = 2,
                         rangeN = c( 100, 150 ),
                         numNonNoisy = 10, numNoisy = 20, numReplicate = 1, 
                         fileName = "chk1")
a <- tmp1$datList$chk1_1
a <- scale(a)  # Standardize for column


## Calculate two prototypes, using kmeans
y <- kmeans(a, 2, iter.max = 200, nstart = 10)
p = y$centers

## Calculate dist:
d <- PrtDist(a, p)
}
\keyword{ classif }
\keyword{ cluster }
\keyword{ optimize }
