
run: gradient_method.o polyhedron.o
	g++ gradient_method.o polyhedron.o -o result

gradient_method.o: gradient_method.cpp polyhedron.hpp
	g++ -c gradient_method.cpp

polyhedron.o: polyhedron.cpp polyhedron.hpp gnuplot.hpp
	g++ -c polyhedron.cpp
