
// gnuplot.hpp
#ifndef _GNUPLOT_H_
#define _GNUPLOT_H_

#include <cstdio>
#include <string>
#include <iostream>

#ifdef WIN32
    #define GNUPLOT_NAME "pgnuplot -persist"
#else
    #define GNUPLOT_NAME "gnuplot -persist"
#endif

using std::string;
using std::cerr;

class Gnuplot
{
public:
    Gnuplot() ;
    ~Gnuplot();
    void operator ()(const string & command); // отправить команду gnuplot
    void operator ()(const char* command);

protected:
    FILE *gnuplotpipe;
};

Gnuplot::Gnuplot()
{
    #ifdef WIN32
        gnuplotpipe = _popen(GNUPLOT_NAME, "w");
    #else
        gnuplotpipe  = popen(GNUPLOT_NAME, "w");
    #endif

    if (!gnuplotpipe)
    {
        cerr << ("Gnuplot not found !");
    }
}
Gnuplot::~Gnuplot()
{
    fprintf(gnuplotpipe,"exit\n");

    #ifdef WIN32
       _pclose(gnuplotpipe);
    #else
        pclose(gnuplotpipe);
    #endif
}
void Gnuplot::operator()(const string & command)
{
    fprintf(gnuplotpipe,"%s\n",command.c_str());
    fflush(gnuplotpipe); //без fflush ничего рисоваться не будет
}
void Gnuplot::operator()(const char* command)
{
    fprintf(gnuplotpipe,"%s\n",command);
    fflush(gnuplotpipe); //без fflush ничего рисоваться не будет
}

void plot_histogramm(double width)
{
    Gnuplot plot;

    plot("width="+std::to_string(width));
    plot("set terminal png size 1500, 800");
    plot("bin(x, s) = s*int(x/s) + width/2");
    plot("set boxwidth width");
    plot("set output 'histogramm_"+std::to_string(width)+".png'");
    plot("set multiplot layout 1,2");
    plot("plot 'data_for_gistogramm.dat' u (bin($1,width)):(1.0) s f w boxes fs solid 0.5 title 'гистограмма по x с шагом "+std::to_string(width)+"' lt rgb 'red'");
    //plot("set output 'y_histogramm_"+std::to_string(width)+".png'");
    plot("plot 'data_for_gistogramm.dat' u (bin($2,width)):(1.0) s f w boxes fs solid 0.5 title 'гистограмма по y с шагом "+std::to_string(width)+"' lt rgb 'red'");
    plot("unset multiplot");
}


#endif // #ifndef _GNUPLOT_H_
