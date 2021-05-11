if ! test -d gifs; then
	mkdir gifs
fi

for im_num in {1..10};
do
	./prog
	if ! test -d gifs/$im_num; then
		mkdir gifs/$im_num;
	fi

	for i in {1..50};
	do 
		./bfgs start_model.txt input.txt $i;
       		gnuplot load 'animator.gp';
       		cp out.png images_for_gif/out_$i.png;
       		cp hists.png images_for_gif/hists_$i.png;
	done
convert -delay 10 -loop 0 images_for_gif/hists_{1..50}.png gifs/$im_num/hist.gif;
convert -delay 10 -loop 0 images_for_gif/out_{1..50}.png gifs/$im_num/cub.gif;
rm images_for_gif/*png
done
