
#include "FentonInout.h"
#include "Fourier.h"

int main()
{
 FentonInout<FourierClass> Inout(FourierSolution,(char *)"Data.dat",(char *)"Convergence.dat",(char *)"Points.dat");
 Inout.Solve();
 Inout.Output((char *)"Surface(Fourier).res",(char *)"Flowfield(Fourier).res",(char *)"Solution(Fourier).res",(char *)"Accuracy.txt");
 printf("\nTouch key to continue "); getch();
 printf("\n\nFinished\n");
 return 0;
}
