\name{gpcapsi}
\alias{gpcapsi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized Principal components algorithm for phase shifting interferometry
}
\description{A generalized principal components algorithm for phase shifting interferometry
    developed by the author. This is the ``low level'' implementation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{gpcapsi(im.mat, trace = 1)
}
\arguments{
  \item{im.mat}{Matrix containing the unmasked pixels from a set of interferograms.
}
  \item{trace}{Print progress of nonlinear solver every trace iterations.
      Use trace=0 for silent operation.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list with the following items:
    \item{phi}{Estimated wrapped phase.}
    \item{mod}{Estimated modulation.}
    \item{phases}{Estimated phase shifts.}
    \item{nlmin}{Return value from Rsolnp::solnp.}
    \item{snr}{An estimate of the S/N of the interferograms.}
    \item{eigen}{Eigenvalues of the crossproduct matrix}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{M. L. Peck}
\note{This is the low level interface to the algorithm. The matrix \code{im.mat} should
    contain the unmasked pixel values from the input interferogram array. No checks
    are made for valid data. This should normally be called through the high level
    function \code{\link{pcafit}}.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{pcapsi}}
        \code{\link{pcafit}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimization}
\keyword{models}% __ONLY ONE__ keyword per line
