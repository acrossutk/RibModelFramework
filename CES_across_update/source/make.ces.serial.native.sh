#! /bin/bash
# script for compiling SEMPPR
#

g++  -static ./CES.c -g -lgsl -lgslcblas -lm -mtune=native -O3 -o $1
