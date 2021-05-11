#!/usr/bin/gnuplot
#set xrange [-0.5:1.5]
#set yrange [-0.5:1.5]
#set zrange[-0.5:1.5]
unset colorbox
set terminal png font "Droid, 14" size 750,750
set key off
unset xlabel
unset ylabel 
unset zlabel
unset title
unset border 
unset tics 
set size square
load "animate.gp"
