set terminal png
set output "error.png"
set grid
plot "error.txt" with boxes fs solid 0.5 lt rgb 'red'

