\name{qVarSelLP}
\alias{qVarSelLP}
\title{
A Mixed Integer Linear Programming Formulation for the Variable Selection Problem  
}
\description{
The function solve the mixed integer linear programming formulation of the variable selection problem.}
\usage{
qVarSelLP(d,
          q,
          binary = FALSE,
          write = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ A 3-dimensional distance matrix, in which d[i,j,k] is the distance between unit i and prototype j, according to variable k.}
  \item{q}{The number of variables to selct.}
  \item{binary}{Set this value to TRUE if you wish to solve the problem with integer variables, or set it to FALSE if you just want to solve the continuous relaxation}
  \item{write}{Set this value to TRUE if you want that the optimization problem is exported in a file called qdistsel.lp}
}
\details{
The function solves the linear problem through the lpSolve solver available through the package lpSolveAPI. The linear programming formulation is the implementation of model F3 described in the paper below (without the median variables).
}
\value{
  \item{status}{The result of the optimization as output of the function solve() of library lpSolveAPI}
  \item{obj}{The value of the objective function}
  \item{x}{The value of the problem variables corresponding to the variable selection: x[j] = 1 means that variable j has been selected, 0 otherwise. If the continuous relaxation has been solved, the vector can contain fractional variables (most likely meaningless).}
}
\references{
S. Benati, S. Garcia Quiles, "A p-median model with distance selection", Working Paper, Universidad Carlos III de Madrid, 2012
}
\author{
Stefano Benati
}
\note{
The computational time to solve an integer programming problem can easily become exponential, therefore be carefull when set variable "binary" to TRUE, as you could wait days even to get the solution of a small scale problem.
Even though the continuos version can contain fractional variables, comparing the objective functions of subroutines qVarSelH and qVarSelLP is a certificate of the solution quality.
}


\seealso{
lpSolveAPI
}
\examples{
 ## Generate random cluster
  
  a <- rbind(cbind(rnorm(5, 0, 1), rnorm(5, 0, 3), rnorm(5, 0, 5)), 
            cbind(rnorm(5, 5, 1), rnorm(5, 5, 3), rnorm(5, 5, 5) ))
            
  ## calculate data prototypes using k-means

  sl2 <- kmeans(a, 2, iter.max = 100, nstart = 2)
  p = sl2$centers
  
  ## calculate distances between observations and prototypes
  ## Remark: d is a 3-dimensions matrix

  d = PrtDist(a, p)
  
  ## Select 2 most representative variables, use heuristic  
  
  lsH <- qVarSelH(d, 2, maxit = 200)

  ## Select 2 variables, use linear relaxation
  
  require(lpSolveAPI)
  lsC <- qVarSelLP(d, 2)
  
  
  ## check optimality
  
  if (abs(lsH$obj - lsC$obj) < 0.001)
      message = "Heuristic Solution is Optimal"
}
\keyword{ classif }
\keyword{ cluster }
\keyword{ optimize }
