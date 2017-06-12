#!/bin/bash

# aditional data processing commands here.

gnuplot << EOP

datafile = "file.txt"

set term jpeg size 1366,768 background "#eeeeee"
set output "600.jpg"

set grid x, y, mxtics, mytics

set mxtic 5
set mytic 5


set xlabel "X, Coordinate"
set ylabel "T, Temperature"

set title "T1 = 400, T2 = 900, T = 300 \n Helium"

set yrange [0:1000]
set xrange [0:1001]

plot datafile index 600 using 1:2 title "U, function" with lines


EOP

xdg-open 600.jpg
