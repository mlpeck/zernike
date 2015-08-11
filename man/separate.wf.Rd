\name{separate.wf}
\alias{separate.wf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Separate wavefronts
%%  ~~function to do ... ~~
}
\description{Separate ``polished in'' from ``instrumental'' aberrations if possible
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
separate.wf(zcm, theta, maxorder = 14)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zcm}{Matrix of observed Zernike coefficients
%%     ~~Describe \code{zcm} here~~
}
  \item{theta}{Vector of rotation angles (in radians)
%%     ~~Describe \code{theta} here~~
}
  \item{maxorder}{Maximum Zernike order to extract
%%     ~~Describe \code{maxorder} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{zcb}{Table of extracted coefficients and standard errors}
%%  \item{sumstats}{Zernike moments and some summary stats from least squares fits}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{M.L. Peck
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Utility}
