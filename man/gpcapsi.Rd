\name{gpcapsi}
\alias{gpcapsi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized Principal components algorithm for phase shifting interferometry
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{gpcapsi(im.mat, trace = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{im.mat}{
%%     ~~Describe \code{im.mat} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{pcapsi}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (im.mat, trace = 1) 
{
    M <- nrow(im.mat)
    nf <- ncol(im.mat)
    svd.im <- svd(im.mat)
    V <- svd.im$v[, 1:3]
    Phi <- svd.im$u[, 1:3] \%*\% diag(svd.im$d[1:3])
    phases <- atan2(V[, 3], V[, 2])
    phases <- wrap(phases - phases[1])
    S <- cbind(rep(1, nf), cos(phases), sin(phases))
    P <- crossprod(V, S)
    fn <- function(par, V, nf) {
        phases <- par[1:nf]
        phases <- wrap(phases - phases[1])
        P <- matrix(par[(nf + 1):(nf + 9)], 3, 3)
        S <- cbind(rep(1, nf), cos(phases), sin(phases))
        norm(V \%*\% P - S, "F")
    }
    sdmin <- Rsolnp::solnp(pars = c(phases, as.vector(P)), fun = fn, 
        LB = c(rep(-pi, nf), rep(-Inf, 9)), UB = c(rep(pi, nf), 
            rep(Inf, 9)), control = list(trace = trace, tol = 1e-10), 
        V = V, nf = nf)
    par <- sdmin$par
    phases <- par[1:nf]
    phases <- wrap(phases - phases[1])
    P <- matrix(par[(nf + 1):(nf + 9)], 3, 3)
    Phi <- Phi \%*\% mpinv(t(P))
    V <- V \%*\% P
    phi <- atan2(-Phi[, 3], Phi[, 2])
    mod <- sqrt(Phi[, 2]^2 + Phi[, 3]^2)
    r2 <- sum(svd.im$d[1:3]^2)/sum(svd.im$d^2)
    list(Phi = Phi, phi = phi, mod = mod/max(mod), phases = phases, 
        V = V, P = P, nlmin = sdmin, snr = sqrt(r2/(1 - r2)), 
        eigen = svd.im$d^2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimization}
\keyword{misc}% __ONLY ONE__ keyword per line
