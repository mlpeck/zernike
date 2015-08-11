\name{plotxs}
\alias{plotxs}
\title{Plot cross-sections (profiles) through a wavefront map.
}
\description{Plots an arbitrary number of cross-sections through
  a wavefront map, with one highlighted.
}
\usage{
plotxs(wf, cp, theta0 = 0, ylim = NULL, N = 4, n = 101,
	col0 = "black", col = "gray", lty = 2)
}
\arguments{
  \item{wf}{A matrix of wavefront values.
}
  \item{cp}{List of pupil parameters as returned by \code{\link{pupil.pars}}.
}
  \item{theta0}{Angle of highlighted profile, in degrees.
}
  \item{ylim}{range of heights to plot.
}
  \item{N}{Number of cross sections.
}
  \item{n}{Number of points for each cross section.
}
  \item{col0}{Highlight color.
}
  \item{col}{Cross section color.
}
  \item{lty}{Line type for plots.
}
}
\details{The cross sections are equally spaced in angle from 0 to \code{pi*(N-1)/N}.
Any angle can be specified for the highlighted profile at \code{theta0}.
}
\value{
none
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{
\code{\link{plot.pupil}} is the main wavefront plotting routine.
}
\keyword{Graphics}
