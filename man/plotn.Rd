\name{plotn}
\alias{plotn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wavefront comparison plots
}
\description{Plot an arbitrary number of wavefronts and all differences.
}
\usage{
plotn(..., labels = NULL, addContours=FALSE, wftype = "net",
	col = rygcb(400), qt = c(0.01, 0.99))
}
\arguments{
  \item{\dots}{List of wavefront estimates as returned by \code{\link{wf_net}}.
}
  \item{labels}{Labels to identify the wavefronts.
}
  \item{addContours}{Boolean to add contours to top row plots}
  \item{wftype}{If the inputs are from \code{wf_net}, one of "net", "smooth",
  "residual".
}
  \item{col}{Color palette for top row of plot
}
  \item{qt}{Quantiles of differences to plot in comparisons.
}
}
\details{
  \dots can be any number of objects containing wavefront estimates as
  returned for example by \code{\link{wf_net}}.

  Wavefronts are displayed on the top row, and differences of all
  pairs on subsequent rows. Grayscale is used to render the
  difference plots, and the color palette given in \code{col}
  is used for the wavefronts.
}
\value{none
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}

\seealso{
  \code{\link{plot.pupil}}
  \code{\link{wf_net}}
}

\keyword{Graphics}

