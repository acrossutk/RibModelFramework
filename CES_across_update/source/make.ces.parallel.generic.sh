#! /bin/bash
# script for compiling SEMPPR
#

g++  -static ./CES.c -g -lgsl -lgslcblas -lm -mtune=generic -fopenmp -O3 -o $1
