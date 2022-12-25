## Version 3.7.5

* Formal writeup of the procedure for recursive generation of Zernike annular polynomials
  is included with the documentation.
* Added functions `rzernike_ann_direct`, `zapm_direct`, `zapm_iso_direct` to calculate radial and Zernike annular polynomials
  directly from formulas published by Mahajan.

## Version 3.7.4

* Removed functions `zapmC` and `zapm_cart` for approximate Zernike Annular polynomials.
* Added C++ functions `zapm`, `zapm_iso`, `zapm_128`, `zapm_iso_128` for nominally exact Zernike Annular polynomials, and their R wrappers.
  These return matrixes of annular Zernikes in extended Fringe or ISO/ANSI ordering. The functions `*_128` use 128 bit floats
  for somewhat higher numerical precision.
* Example code is included for the above functions.
* Also new are `rzernike_ann` and `rzernike_ann_128` to calculate the radial annular Zernikes. Finally, function `gol_welsch` calculates
  nodes and weights for Gauss-Legendre quadrature in the interval (eps^2, 1). This would normally not be called by the user, but it
  is in the namespace and documented. All of these are mostly for testing purposes.
* added boolean variable `usecirc` to the list of options returned by `psfit_options()`. If set to `TRUE` functions that use Zernikes for
  wavefront representations will use Zernike circle polynomials even for obstructed apertures. For the default value `FALSE` functions like
  `wf_net()` and `pupil` will toggle between Zernikes and annular Zernikes depending on whether the obstruction fraction is > 0.
* The boost library in package `BH` is now linked to and therefore must be installed (?? IS THIS TRUE FOR USERS ??)

## Version 3.7.2

* `crop` will return a square matrix with an even number of pixels per side. Side length is optionally specifiable. Interferograms should be completely contained within frames with some padding on all sides.  This function does absolutely no error checking.
* added option `lpbrcut` to list of available phase unwrapping routines set in `psfit_options()$puw_alg`. This will result in `lppuw::brcutpuw` being called as the phase unwrapper. Option `brcut` will result in calling `zernike::brcutpuw`.
* Consistent calls to phase unwrapping routine in `fftfit` and `vortexfit`.


## Version 3.7.1

New in this version: support for Zernike Annular polynomials approximated by orthogonalizing a matrix of Zernike polynomials:

* Added C++ functions `zapmC` and `zapm_cart` to create matrixes of Annular Zernikes in extended Fringe or ISO/ANSI sequence.
* Modified the functions `wf_net`, `fitzernikes`, and `pupil` to use annular Zernikes for annular (or obstructed) wavefronts. This is done automatically and involves no user visible changes.
* Expanded help for some topics and added examples for `wf_net` and `zpm_cart`.

## Version 3.7.0

A number of additions and changes have been made for this release, and some of the changes may break previously working code. This list _may_ not be comprehensive:

* Added `gradzpm_cart`, `zpm_cart`, and `norm_zpm` to create matrixes of Zernike polynomials and their derivatives in ISO/ANSI sequence and Cartesian coordinates. 
  See the documentation for further details.
* Added the item `isoseq` to the list in `psfit_options()` to tell other routines when to use `zpm_cart` instead of `zpm`. This is checked in the function `wf_net()` and if true `zpm_cart` is used in the Zernike fitting.
* Added function `makezlist.iso` to make a list of Zernike polynomial indexes in ISO/ANSI sequence.
* Added print, summary, and plot methods for the return lists from `psifit()`, `vortexfit()`, and `fftfit()`. These have been assigned the S3 class `wf_fitted`.
* **Potentially important!** Arguments to `pupil` have changed: instead of passing `zlist` as argument I now use `maxorder` which is in turn passed to `zpm` or `zpm_cart`. If the first element of zcoef is the piston term the piston entry in the argument list should be left null, which is the default. If zcoef starts with tilts, which it will be if returned by `wf_net` or one of the routines that call _it_, piston _must_ be non-null (usually should be 0). The function `pupil.arb()` continues to use `zlist` to provide a list of Zernike indexes since the usual use case would be to create a synthetic wavefront with just a few aberrations.
* The bonus functions `startest()` and `synth.interferogram()` have modified arguments and the logic of the code has been improved (maybe).
* I'm continuing the process of splitting out source files organized into topically related functions. This only applies to the source distribution.

## Release notes 3.6.0

Several changes were made for the current release. Most of these are "behind the scenes" and mostly made for performance or maintainability reasons. These include:

* C code for phase unwrapping routines has been converted to C++ for use with Rcpp. This mostly involved changing pointers to references in function calls.
* A threaded version of the code for filling matrixes of Zernike polynomials has been added. The number of threads used is set by the argument `nthreads` to `psfit_options()`. This is set by default to half the number of cores detected because all of the PCs in my collection have multithreading CPUs that report double the number of physical cores actually present. In practice threading might or might not be beneficial. Speed gains are more likely if very large numbers of ZP values are needed.

Some user visible changes:

* Added function `vortexfit()` implementing the vortex aka spiral quadrature transform algorithm of Larkin et al. This is the same algorithm with possibly slightly different implementation as used in DFTFringe.
* Added function `circle.hough()` to detect circular interferogram edges using the Hough circle transform. This is experimental and is not curently used by any other routine.
* The function `circle.pars()` for automated edge detection now uses the function `nlsrob()` from package robustbase for circle parameter estimation. An optional but enabled by default refinement step used in previous versions of the algorithm has been removed.

