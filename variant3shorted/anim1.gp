
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
unset output
set output "out.png"
splot "out.txt" with lines lw 1 lt rgb 'black', "draw_cor_edges.txt" with lines lw 1 lt rgb  'red'
