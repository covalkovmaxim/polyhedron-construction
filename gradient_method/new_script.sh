for i in {1..10}; 
do 
	./prog; 
	./bfgs start_model.txt input.txt 10000; 
	gnuplot load 'animator.gp';
       	cp input.png images/$i/input.png; 
	cp out.png images/$i/out.png;
       	cp input.txt images/$i/input.txt; 
	cp out.txt images/$i/out.txt;
	cp hists.png images/$i/hists.png;
done
