\name{rygcb}
\alias{rygcb}
\alias{rygcbm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A better rainbow.
}
\description{
Produces a rainbow color palette with colors ranging from
\code{"red"} to \code{"blue"} or \code{"magenta"}. 
Perceptual uniformity should be superior to \R's \code{rainbow}.
}
\usage{
rygcb(n)
rygcbm(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of color levels}
}
\details{The palette is created using \code{colorRampPalette}.}
\value{A vector of colors.}
\author{
M.L. Peck
}
\note{The call to \code{colorRampPalette} sets \code{space="Lab"}
and \code{interpolate="spline"} with the intent of creating
a more perceptually uniform rainbow.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{grey256}}
}
\examples{
plotsp <- function(spectrum) {
	sl <- length(spectrum)
	rgbv <- col2rgb(spectrum)
	plot((0:(sl-1))+0.5, rgbv[1,], type="l", col="red", xlim=c(0,sl),ylim=c(0,300),xlab="Index",ylab="Channel value")
	points((0:(sl-1))+0.5, rgbv[2,], type="l", col="green")
	points((0:(sl-1))+0.5, rgbv[3,], type="l", col="blue")
	grid()
	rect(0:(sl-1), 260, 1:sl, 300, col=spectrum, density=NA)
}
plotsp(rygcb(400))
X11()
plotsp(rygcbm(500))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Graphics}

