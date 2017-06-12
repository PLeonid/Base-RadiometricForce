#!/bin/bash

# aditional data processing commands here.

gnuplot << EOP

datafile = "file.txt"

set term gif animate optimize delay 0 size 1366,768 background "#eeeeee"
set output "data.gif"

set grid x, y, mxtics, mytics

set mxtic 5
set mytic 5


set xlabel "X, Coordinate"
set ylabel "T, Temperature"

set title "Название \n Твоего Графика"

set yrange [0:1000]
set xrange [0:1001]

do for [i=1:600] {
plot datafile index i using 1:2 title "U, function" with lines
}

EOP

xdg-open data.gif
