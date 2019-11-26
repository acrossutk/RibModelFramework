#!/bin/bash

#Instructions on how to build and install ribModel package from a git repository

#You may need to reinstall Rcpp
#If you want a subwindow in figures you need to install packages
##  survival
##  Matrix
##  Hmisc

##build the package outside the repository RibModelFramework
## The tar.gz name for package is automatically set
#
# Global install
#R CMD build RibModelFramework; sudo MAKE="make -j4" R CMD INSTALL $(ls -t ribModel_*.tar.gz| head -1)
#
#  Local install
R CMD build RibModelFramework; MAKE="make -j4" R CMD INSTALL -l ~/Rlibs $(ls -t ribModel_*.tar.gz| head -1)
