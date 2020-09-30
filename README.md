Welcome to my R package for interferogram analysis and manipulation of Zernike polynomials. If you are new to R this README gives installation and basic usage instructions for Windows users. There are a number of required and optional packages to be installed, so please follow these instructions carefully. You need:

* A recent version of R itself (v3.6.3 or later is required). The current version at this writing is 4.0.2. If you're making a fresh installation just go with the current release. The windows binary installer can be downloaded from http://cran.r-project.org. Installation works like any Windows program -- I like to have the installer create desktop shortcuts and edit its properties to have it start up in the directory where I've stored data. If you are running a 64 bit version of Windows both 32 and 64 bit binaries will be installed. The 64 bit version's performance is slightly better and it can use all available memory, which may be useful if you work with large interferogram images. The only customization I make at install time is to select "SDI" windowing mode. This will display graphs in separate windows on the desktop

* Download the Windows binary of this package from the [releases](../releases) section. Do not unzip it!

* The following steps are done within R itself. It may be useful to run R as an administrator for these initial steps, but if you don't packages should still install somewhere accessible to you as an ordinary user. When you start R a console window opens with a small number of menu items and a text entry area. Commands are typed at the ">" prompt. You should have an active internet connection for the following steps:

* Install the following packages. This can be done with the menu item Packages/Install package(s)... The first time you select this in a session it will prompt you for a mirrored download site. Next it will display a selection box with all available packages. The ones you need are:

    + Rcpp. This package also has a large number of dependencies that should be automatically downloaded and installed.
    
    + RcppParallel.
    
    + RcppArmadillo.
    
    + rgl. This package provides interactive 3D graphics and is strictly speaking optional.
    
    + robustbase. The function `nlsrob` is used as the new default parameter estimation routine by the automated edgefinding function `circle.pars` and is also strictly speaking optional.
    
* You need not install these initially but they may be useful for certain tasks.
    
    + clue. The function `solve_LSAP` is used by the branch cut algorithm for phase unwrapping implemented in the function `brcutpuw`. This may outperform the default phase unwrapping routine in some situations.
    
    + data.table, dplyr, pixmap, mvtnorm.

* Now install package "zernike". This can be done with the menu item Packages/Install from local zip file... At the prompt just navigate to wherever you saved the zip file and select it. If it installs successfully a brief message will be sent to the console. You may get a warning message if you're running a different version of R from the one the package was built in. You can probably ignore this -- if the demos run you're in good shape.

* Three demo programs are included with the program. To run them enter the following commands:
```
demo(psiest, package="zernike", ask=FALSE)
```
```
demo(winfit, package="zernike", ask=FALSE)
```
```
demo(pcaest, package="zernike", ask=FALSE)
```
These demos illustrate basic PSI analysis, a sliding window analysis of a multiple cycle PSI sequence, and some alternative algorithms for generalized phase shifting interferometry.

* The main PSI analysis routine is `psifit()`. Type `example(psifit, package="zernike", ask=FALSE)` to see an example of its use. This uses the same data set as the demo programs.

* To access R's help in html format type `help.start()`. This will open up a window in your default browser with a help page containing links to the help files for all installed packages on your system as well as manuals, FAQs, and help files for R itself.

* You can exit R by typing q() or just closing the console window. You will be prompted whether to save the workspace. If you enter y any data in your workspace will be saved in binary format to a file named .RData, and commands entered in the current session will be saved to an ordinary text file named .Rhistory.

There is a vast quantity of documentation and literature about the R system. At a minimum you should read the [R for Windows FAQ](https://cran.r-project.org/bin/windows/base/rw-FAQ.html) and the introduction to R included as a PDF file in the software installation.

## Release notes

Several changes were made for the current release. Most of these are "behind the scenes" and mostly made for performance or maintainability reasons. These include:

* C code for phase unwrapping routines has been converted to C++ for use with Rcpp. This mostly involved changing pointers to references in function calls.
* A threaded version of the code for filling matrixes of Zernike polynomials has been added. The number of threads used is set by the argument `nthreads` to `psfit_options()`. This is set by default to half the number of cores detected because all of the PCs in my collection have multithreading CPUs that report double the number of physical cores actually present. In practice threading might or might not be beneficial. Speed gains are more likely if very large numbers of ZP values are needed.

Some user visible changes:

* Added function `vortexfit()` implementing the vortex aka spiral quadrature transform algorithm of Larkin et al. This is the same algorithm with possibly slightly different implementation as used in DFTFringe.
* Added function `circle.hough()` to detect circular interferogram edges using the Hough circle transform. This is experimental and is not curently used by any other routine.
* The function `circle.pars()` for automated edge detection now uses the function `nlsrob()` from package robustbase for circle parameter estimation. An optional but enabled by default refinement step used in previous versions of the algorithm has been removed.

***
PSI data are courtesy of Vladimir Galogaza. Steve Koehler provided valuable programming advice.

