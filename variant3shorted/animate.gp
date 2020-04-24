unset output
set output "out.png"
splot "out.txt" with lines lw 1 lt rgb 'black', "draw_cor_edges.txt" with lines lw 1 lt rgb  'red'
unset output
set output "input.png"
splot "cub.txt" with lines lw 1 lt rgb 'black', "draw_cor_edges.txt" with lines lw 1 lt rgb  'red'
unset output

