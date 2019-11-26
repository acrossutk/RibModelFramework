#! /bin/bash
# script for compiling CES
#

./make.ces.parallel.generic.sh ../bin/CES_generic_parallel;
./make.ces.serial.generic.sh ../bin/CES_generic_serial;
./make.ces.serial.native.sh ../bin/CES_serial;
./make.ces.parallel.native.sh ../bin/CES;
