
#include "FentonInout.h"
#include "Stokes.h"

int main()
{
 FentonInout<StokesClass> Inout(StokesSolution,(char *)"Data.dat",(char *)"Points.dat");
 Inout.Solve();
 Inout.Output((char *)"Surface(Stokes).res",(char *)"Flowfield(Stokes).res",(char *)"Solution(Stokes).res",(char *)"Accuracy.txt");

 printf("\nTouch key to continue "); getch();
 printf("\n\nFinished\n");
 return 0;
}
