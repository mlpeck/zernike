Welcome to my R package for interferogram analysis and manipulation of Zernike polynomials. If you are new to R this README gives installation and basic usage instructions for Windows users. There are a number of required and optional packages to be installed, so please follow these instructions carefully. You need:

* A recent version of R itself (v4.1.0 or later is required). The current version at this writing is 4.4.1. If you're making a fresh installation just go with the current release. The windows binary installer can be downloaded from http://cran.r-project.org. Installation works like any Windows program -- I like to have the installer create desktop shortcuts and edit its properties to have it start up in the directory where I've stored data. The only customization I make at install time is to select "SDI" windowing mode. This will display graphs in separate windows on the desktop.

* Download the Windows binary of this package from the [releases](../releases) section. Do not unzip it!

* The following steps are done within R itself. When you start R a console window opens with a small number of menu items and a text entry area. Commands are typed at the ">" prompt. You should have an active internet connection for the following steps:

* Install the following packages. This can be done with the menu item Packages/Install package(s)... The first time you select this in a session it will prompt you for a mirrored download site. Next it will display a selection box with all available packages. The ones you need are:

    + Rcpp. This package also has a large number of dependencies that should be automatically downloaded and installed.
    
    + RcppParallel.
    
    + RcppArmadillo.
    
    + BH.
    
* Optional but highly recommended packages
    
    + rgl. This package provides interactive 3D graphics.
    
    + robustbase. The function `nlsrob` is used as the new default parameter estimation routine by the automated edgefinding function `circle.pars`.
    
    + tinytable. Used to produce html formatted reports for the fringe analysis routines.
    
    + fftwtools. If installed the FFT based fringe analysis routines `vortexfit` and `fftfit` will use `fftw` for FFTs instead of the base R function. This may slightly improve the execution time of these functions.
    
* You need not install these initially but they may be useful for certain tasks.
    
    + clue. The function `solve_LSAP` is used by the branch cut algorithm for phase unwrapping implemented in the function `brcutpuw`. This may outperform the default phase unwrapping routine in some situations.
    
    + data.table, dplyr, pixmap, mvtnorm.
    
    + my package lppuw, available on github at [github.com/mlpeck/lppuw](https://github.com/mlpeck/lppuw), which contains some advanced methods of phase unwrapping. This also requires a package named [rcbc](https://dirkschumacher.github.io/rcbc/). A Windows binary is provided in the releases section of `lppuw` for anyone who needs it.
    

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

* The main PSI analysis routine is `psifit()`. Type `example(psifit, package="zernike", ask=FALSE, echo=FALSE)` to see an example of its use. This uses the same data set as the demo programs.

* I now have substantial examples for several high and not so high level functions in the package. Besides `psifit()` these include `vortexfit()`, `wf_net()`, `zpm_cart()`, `zapm()`, `zapm_iso()`, `zapm_128()`, and `zapm_iso_128()`. These are run in the same way as described for `psifit()` above.

* To access R's help in html format type `help.start()`. This will open up a window in your default browser with a help page containing links to the help files for all installed packages on your system as well as manuals, FAQs, and help files for R itself.

* You can exit R by typing q() or just closing the console window. You will be prompted whether to save the workspace. If you enter y any data in your workspace will be saved in binary format to a file named .RData, and commands entered in the current session will be saved to an ordinary text file named .Rhistory.

There is a vast quantity of documentation and literature about the R system. At a minimum you should read the [R for Windows FAQ](https://cran.r-project.org/bin/windows/base/rw-FAQ.html) and the introduction to R included as a PDF file in the software installation.

***
PSI data are courtesy of Vladimir Galogaza. Steve Koehler provided valuable programming advice.

