run: gradient_method.o polyhedron.o
	g++ gradient_method.o polyhedron.o -o result

gradient_method.o: gradient_method.cpp polyhedron.hpp
	g++ -c gradient_method.cpp

polyhedron.o: polyhedron.cpp polyhedron.hpp gnuplot.hpp
	g++ -c polyhedron.cpp
		
bfgs: nlopt.o polyhedron.o 
	g++ -L /usr/local/lib -I /usr/local/include nlopt_test.o polyhedron.o -o bfgs -lnlopt

nlopt.o: nlopt_test.cpp polyhedron.hpp
	g++ -c -L /usr/local/lib -I /usr/local/include nlopt_test.cpp -lnlopt
