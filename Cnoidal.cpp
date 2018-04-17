
#include "FentonInout.h"
#include "Cnoidal.h"

int main()
{
 FentonInout<CnoidalClass> Inout(CnoidalSolution,(char *)"Data.dat",(char *)"Points.dat");
 Inout.Output((char *)"Surface(Cnoidal).res",(char *)"Flowfield(Cnoidal).res",(char *)"Solution(Cnoidal).res",(char *)"Accuracy.txt");

 printf("\nTouch key to continue "); getch();
 printf("\n\nFinished\n");
 return 0;
}
