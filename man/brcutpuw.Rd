\name{brcutpuw}
\alias{brcutpuw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Branch cut algorithm for phase unwrapping
}
\description{
Solves a modification of the assignment problem to minimize the total length of branch cuts.
}
\usage{
brcutpuw(phase, pen = 0, details = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phase}{
Matrix containing the wrapped phase map
}
  \item{pen}{
Optional penalty value for connecting a residue to an edge
}
  \item{details}{
boolean: if TRUE return some extra details for diagnostics
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix containing the unwrapped wavefront. If details==TRUE a named list starting with the unwrapped wavefront in puw.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
M.L. Peck
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qpuw}}
\code{\link{idiffpuw}}
}
\examples{
  set.seed(1234)
  PW <- wrap(matrix((0:100)*pi/10,101,101))
  ## need a border of NA's
  PW <- cbind(rep(NA,101), PW, rep(NA,101))
  PW <- rbind(rep(NA,103), PW, rep(NA,103))
  PW <- PW + rnorm(103^2)
  mtext(rmap(PW, plot=TRUE))
  PU <- brcutpuw(PW, details=TRUE)
  image(1:103, 1:103, PU$bcuts, col="blue", add=TRUE)
  X11()
  image(PU$puw, col=grey256, asp=1, useRaster=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}% use one of  RShowDoc("KEYWORDS")
\keyword{misc}% __ONLY ONE__ keyword per line
