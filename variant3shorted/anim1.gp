
unset colorbox

set key off
unset xlabel
unset ylabel 
unset zlabel
unset title
unset border 
unset tics 
set size square
splot "out.txt" with lines lw 1 lt rgb 'black', "draw_cor_edges.txt" with lines lw 1 lt rgb  'red'
