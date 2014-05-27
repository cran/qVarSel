\name{qVarSelH}
\alias{qVarSelH}
\title{
Selection of Variables for Clustering or Data Dimension Reduction
}
\description{ The function implements the q-Vars heuristic described in the reference below. Given a 3-dimension matrix D, with d[i,j,k] being the distance between statistic unit i and prototype j measured through variable k, the function calculates the set of variables of cardinality q that mostly explains the prototypes.
}
\usage{
qVarSelH(d,
         q,
         maxit = 100)
}
\arguments{
  \item{d}{A numeric 3-dimensional matrix where elements d(i,j,k) are the     distances between observation i and cluster center/prototype j, that are measured through variable k.
  }
  \item{q}{A positive scalar, that is the number of variables to select}
  \item{maxit}{A positive scalar, that is the maximum number of iteration allowed}
}
\details{
The heuristic repeatedly selects a set of variables and then allocates units to prototypes, while a local optimum is reached. Random restart is used to continue the search until the maximum number of iteration is reached. 
}
\value{
\item{obj}{The value of the objective function}
\item{x}{A 0-1 vector describing wheter variable k is selected: If x[k] = 1 then k is selected}
\item{ass}{A vector of assignment of units to clusters: if ass[i] = j then unit i is assigned to the cluster represented by center/prototype j }
\item{bestit}{The iteration in which the optimal solution is found}
}
\references{
S. Benati, S. Garcia Quiles, J. Puerto "Optimization Methods to Select Variables for Clustering", Working Paper, Universidad de Sevilla, 2014
}
\author{
Stefano Benati
}
\note{
The methodology is heuristic and some steps are random. It may be the case that different runs provide different solutions. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
qVarSelLP
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


## Calculate two prototype, using kmeans
y <- kmeans(a, 2, iter.max = 200, nstart = 10)
p = y$centers

## Calculate dist:
d <- PrtDist(a, p)
           
## Calculate Best 10 variables:
lsH <- qVarSelH(d, 10, maxit = 200)

}
\keyword{ classif }
\keyword{ cluster }
\keyword{ optimize }
