#! /bin/bash
# script for compiling CES
#

g++  -static ./CES.c -g -lgsl -lgslcblas -lm -mtune=native -o $1
