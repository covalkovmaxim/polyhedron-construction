#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<time.h>
#include <unistd.h>
using namespace std;

int main()
{
  srand(time(NULL));
  FILE*fp=fopen("input.txt","w");
  fprintf(fp,"8\n2 5\n%f %f %f\n",
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);
  fprintf(fp,"%f %f %f\n   \n    \n",
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);
  

  fprintf(fp,"3 5\n%f %f %f\n",
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);
  fprintf(fp,"%f %f %f\n   \n    \n",
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);
  
  fprintf(fp,"4 5\n%f %f %f\n",
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);
  fprintf(fp,"%f %f %f\n   \n    \n",
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);
 
  fprintf(fp,"1 5\n%f %f %f\n",
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);

  fprintf(fp,"%f %f %f\n   \n    \n",
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);

   ////////////////////////////////////////////////////  
  
  fprintf(fp,"0 2\n%f %f %f\n",
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);
  fprintf(fp,"%f %f %f\n   \n    \n",
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);


  fprintf(fp,"0 3\n%f %f %f\n",
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);
  fprintf(fp,"%f %f %f\n   \n    \n",
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);
 
  fprintf(fp,"0 4\n%f %f %f\n",
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);
  fprintf(fp,"%f %f %f\n   \n    \n",
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);

  fprintf(fp,"0 1\n%f %f %f\n",
  1.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);

  fprintf(fp,"%f %f %f\n   \n    \n",
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.,
  0.+(2.*(rand()%2)-1.)*(rand()%1000)/2000.);

  fclose(fp);
  sleep(1);
  return 0;
}
