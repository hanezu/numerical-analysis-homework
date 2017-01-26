#!/bin/sh
# don't know why permission denied -> use python as script
make
./program.out > output.dat
gnuplot -e "splot 'output.dat'"
gnuplot -e 'splot "output.dat"'