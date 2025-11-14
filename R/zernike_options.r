## options for fringe analysis routines, wavefront fitting and wavefront display

zernike_options <- function(colors=topo.colors(256), refine=TRUE, puw_alg = "qual", fringescale=1,
                          wt=NULL, bgsub=TRUE,
                          maxiter=20, ptol=1.e-4, trace=1, nzcs = 2,
                          zc0=6:7,
                          satarget=c(0,0), astig.bath=c(0,0),
                          maxorder=14, isoseq=FALSE,
                          usecirc=FALSE, ext_prec=FALSE,
                          nthreads=parallel::detectCores()/2,
                          plots=TRUE, crop=FALSE,
                          sigma_ed=3., qt_ed=0.995, use_fftw=FALSE) {
  list(colors=colors,        ## color palette for wavefront plots
  plots=plots,               ## plot results in wavefront summaries?
  refine=refine,             ## do 2nd pass through psi algorithm?
  puw_alg=puw_alg,           ## phase unwrapping algorithm
  fringescale=fringescale,   ## waves per fringe in interferograms (usually 1 or 1/2)
  wt=wt,                     ## per frame weights in least squares PSI algorithm
  bgsub=bgsub,               ## subtract a background estimate in PC based PSI algorithms?
  maxiter=maxiter,           ## maximum no iterations for iterative PSI algorithms
  ptol=ptol,                 ## convergence tolerance for iterative PSI algorithms
  trace=trace,               ## some iterative algorithms can return info while working
  nzcs=nzcs,                 ## number of zernike coefficients to treat as variable in tiltpsiC
  zc0=zc0,                   ## coefficients to remove in net wavefront
  satarget=satarget,         ## target SA for numerical nulling
  astig.bath=astig.bath,     ## amount of astigmatism from Bath interferometer geometry
  maxorder=maxorder,         ## maximum order for Zernike polynomial fitting
  isoseq=isoseq,             ## use ISO sequenced Zernikes?
  usecirc=usecirc,           ## use circular Zernikes even for obstructed apertures?
  ext_prec=ext_prec,         ## extended precision annular Zernikes?
  nthreads=nthreads,         ## no threads to use with zpmCP
  crop=crop,                 ## crop wavefront related matrixes?
  sigma_ed=sigma_ed,         ## smoothing parameter for circle.auto
  qt_ed=qt_ed,               ## threshold parameter for circle.auto
  use_fftw=use_fftw)         ## use fftw based 2D FFT routines?
}
