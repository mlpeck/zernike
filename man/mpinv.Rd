\name{mpinv}
\alias{mpinv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Moore-Penrose generalized inverse
%%  ~~function to do ... ~~
}
\description{Computes the Moore-Penrose generalized inverse of a matrix using singular
    value decomposition.}
\usage{
mpinv(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A matrix}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Matrix containing the generalized inverse. 
    If X is an n x m matrix the return will have dimension m x n.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{M. L. Peck
}
\note{The threshold for determining if a matrix is rank deficient is
    \code{eps <- .Machine$double.eps * max(dim(X)) * S$d[1]}
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
    X <- matrix(rnorm(18), 6, 3) ## this should be full rank almost always
    mpinv(X) %*% X
    
    X <- matrix(1:18, 6, 3)  ## this is not
    mpinv(X) %*% X
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{array}
\keyword{math}% __ONLY ONE__ keyword per line
