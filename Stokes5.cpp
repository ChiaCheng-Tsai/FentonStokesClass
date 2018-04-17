
// compared with Stokes, only the Output filenames are different

#include<iostream>

#include "FentonInout.h"
#include "Stokes.h"

int main()
{
 FentonInout<StokesClass> Inout(StokesSolution,(char *)"Data.dat",(char *)"Points.dat");
 Inout.Output((char *)"Surface(Stokes5).res",(char *)"Flowfield(Stokes5).res",(char *)"Solution(Stokes5).res",(char *)"Accuracy.txt");
 printf("\nTouch key to continue "); getch();
 printf("\n\nFinished\n");
 return 0;
}
