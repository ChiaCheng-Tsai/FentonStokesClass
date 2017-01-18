
#ifndef Fourier_H
#define Fourier_H
// Stokes theory calculations

#include <math.h>
#include <stdio.h>
#include <process.h>
#include <string.h>
#include <conio.h>
#include <stdlib.h>
#include <assert.h>
#define	ANSI

#define	pi 3.14159265358979324

#include "FentonInout.h"
#include "FentonSolve.h"

class FourierClass
{
 private:
  int num;
  double height,Hoverd;
  double criter;
  double *coeff,*cosa,*rhs1,*rhs2,*sina,*FourierTanh;
  double **sol;
 public:
  FourierClass(FentonInoutClass<FourierClass> & );
  void init(int Inoutn,char* InoutCase,double* Inoutz,double InoutCurrent,int InoutCurrent_criterion);
  double Eqns(double *rhs,int Inoutn,char* InoutCase,double* Inoutz,double InoutCurrent,int InoutCurrent_criterion);
  double Newton(int Inoutn,char* InoutCase,double* Inoutz,double InoutCurrent,int InoutCurrent_criterion);
};

// read data (n, Case, T, L, MaxH, Current_criterion, Current) from Inout
// and set Method, z[], B[], Y[], H, kd, kH, and SU to Inout

FourierClass::FourierClass(FentonInoutClass<FourierClass> & Inout)
{
 int	i, j, iter, m, ns;
 double	dhe, dho, error, **CC;
 double sum;

 if(Inout.n==1)	assert(false); // This is for testing first-order theories

 sprintf(Inout.Method, "Fourier method with %d terms in series", Inout.n);
 printf( "\n# Solution by Fourier method with %d terms in series", Inout.n);

 num=2*Inout.n+10;
 dhe=Inout.Height/Inout.nstep;
 dho=Inout.MaxH/Inout.nstep;

 CC = new double* [num+1];

 for(i=0;i<num+1;i++)
  CC[i] = new double [num+1];

 for ( j=1; j <=num ; ++j)
 {
  for ( i=1; i <=num ; ++i)CC[j][i] = 0.;
  CC[j][j] = 1.;
 }

 Inout.z = new double[num+1];
 Inout.Y = new double[Inout.n+1]; // bug fixed by Tsai
 Inout.B = new double[Inout.n+1];

// see page 359 in Fenton (1988) for the definition of z and B
// see Eq. (3.5) in Fenton (1999) for the definition of B
// see page 363 in Fenton (1988) for the definition of Y

 rhs1 = new double[num+1];
 rhs2 = new double[num+1];
 coeff = new double[num+1];
 cosa = new double[2*Inout.n+1];
 sina = new double[2*Inout.n+1];
 FourierTanh = new double[Inout.n+1];

 sol = new double* [num+1];
 for(i=0;i<num+1;i++)
  sol[i] = new double [3];

//	Commence stepping through steps in wave height
 for ( ns = 1 ; ns <= Inout.nstep ; ns++ )
 {
  height=ns*dhe;
  Hoverd=ns*dho;
  fprintf(Inout.monitor,"\n\nHeight step %2d of %2d", ns, Inout.nstep);

//	Calculate initial linear solution
  if(ns <= 1)
    init(Inout.n,Inout.Case,Inout.z,Inout.Current,Inout.Current_criterion);
//	Or, extrapolate for next wave height, if necessary

  else
   for ( i=1 ; i <= num ; i++ )
    Inout.z[i]=2.*sol[i][2]-sol[i][1];

//	Commence iterative solution
  for (iter=1 ; iter <= Inout.number ; iter++ )
  {
   fprintf(Inout.monitor,"\nIteration%3d:", iter);

//	Calculate right sides of equations and differentiate numerically
//	to obtain Jacobian matrix, then solve matrix equation

   error = Newton(Inout.n,Inout.Case,Inout.z,Inout.Current,Inout.Current_criterion);

//	Convergence criterion satisfied?

   fprintf(stdout," Mean of corrections to free surface: %8.1e", error);

   if(ns == Inout.nstep)	criter = 1.e-10 ;
   else			criter = Inout.crit;

   if((error < criter * fabs(Inout.z[1]))  && iter > 1 ) break;
   if(iter == Inout.number)
   {
    fprintf(stdout,"\nNote that the program still had not converged to the degree specified\n");
   }

//	Operations for extrapolations if more than one height step used

   if(ns == 1)
    for ( i=1 ; i<=num ; i++ )
     sol[i][2] = Inout.z[i];
   else
    for ( i=1 ; i<=num ; i++ )
    {
     sol[i][1] = sol[i][2];
     sol[i][2] = Inout.z[i];
    }
  }

//	Fourier coefficients (for surface elevation by slow Fourier transform)

  for ( Inout.Y[0] = 0., j = 1 ; j <= Inout.n ; j++ )
  {
   Inout.B[j]=Inout.z[j+Inout.n+10];
   sum = 0.5*(Inout.z[10]+Inout.z[Inout.n+10]*pow(-1.,(double)j));
   for ( m = 1 ; m <= Inout.n-1 ; m++ )
    sum += Inout.z[10+m]*cosa[(m*j)%(Inout.n+Inout.n)];
   Inout.Y[j] = 2. * sum / Inout.n;
  }
 } // End stepping through wave heights

 delete [] rhs1;
 delete [] rhs2;
 delete [] coeff;
 delete [] cosa;
 delete [] sina;
 delete [] FourierTanh;

 for(i=0;i<num+1;i++)
  delete [] sol[i];
 delete [] sol;

 for(i=0;i<num+1;i++)
  delete [] CC[i];
 delete [] CC;

}

// **************************************************
// CALCULATE INITIAL SOLUTION FROM LINEAR WAVE THEORY
// **************************************************

void FourierClass::init(int Inoutn,char* InoutCase,double* Inoutz,double InoutCurrent,int InoutCurrent_criterion)
{
 int i;
 double a, b, t;

 iff(InoutCase,Period)
 {
  a=4.*pi*pi*height/Hoverd;
  b=a/sqrt(tanh(a));
  t=tanh(b);
  Inoutz[1]=b+(a-b*t)/(t+b*(1.-t*t));
 }
 else
  Inoutz[1]=2.*pi*height/Hoverd;

 Inoutz[2]=Inoutz[1]*Hoverd;
 Inoutz[4]=sqrt(tanh(Inoutz[1]));
 Inoutz[3]=2.*pi/Inoutz[4];

 if(InoutCurrent_criterion==1)
 {
  Inoutz[5]=InoutCurrent*sqrt(Inoutz[2]);
  Inoutz[6]=0.;
 }
 else
 {
  Inoutz[6]=InoutCurrent*sqrt(Inoutz[2]);
  Inoutz[5]=0.;
 }

 Inoutz[7]=Inoutz[4];
 Inoutz[8]=0.;
 Inoutz[9]=0.5*Inoutz[7]*Inoutz[7];
 cosa[0]=1.;
 sina[0]=0.;
 Inoutz[10]=0.5*Inoutz[2];

 for( i=1 ; i<=Inoutn ; i++ )
 {
  cosa[i]=cos(i*pi/Inoutn);
  cosa[i+Inoutn]=cos((i+Inoutn)*pi/Inoutn);
  sina[i]=sin(i*pi/Inoutn);
  sina[i+Inoutn]=sin((i+Inoutn)*pi/Inoutn);
  Inoutz[Inoutn+i+10]=0.;
  Inoutz[i+10]=0.5*Inoutz[2]*cosa[i];
 }
 Inoutz[Inoutn+11]=0.5*Inoutz[2]/Inoutz[7];

 for( i=1 ; i<=9 ; i++ )
  sol[i][1] = Inoutz[i];
 for( i=10 ; i<=num ; i++ )
  sol[i][1] = 0.;
}

//	EVALUATION OF EQUATIONS.

double FourierClass::Eqns(double *rhs,int Inoutn,char* InoutCase,double* Inoutz,double InoutCurrent,int InoutCurrent_criterion)
{
 int i, j, m, it, nm;
 double c, e, s, u, v;
 double psi; // added by Tsai

 rhs[1]=Inoutz[2]-Inoutz[1]*Hoverd;

 iff(InoutCase,Wavelength)
  rhs[2]=Inoutz[2]-2.*pi*height;
 else
  rhs[2]=Inoutz[2]-height*Inoutz[3]*Inoutz[3];

 rhs[3]=Inoutz[4]*Inoutz[3]-pi-pi;
 rhs[4]=Inoutz[5]+Inoutz[7]-Inoutz[4];
 rhs[5]=Inoutz[6]+Inoutz[7]-Inoutz[4];

 rhs[5]=rhs[5]-Inoutz[8]/Inoutz[1];

 for (i=1; i<=Inoutn; i++ )
 {
  coeff[i]=Inoutz[Inoutn+i+10];
  FourierTanh[i] = tanh(i*Inoutz[1]);
 }

 it=6;
 if(InoutCurrent_criterion==1)it=5;
 rhs[6]=Inoutz[it]-InoutCurrent*sqrt(Inoutz[1]); // Correction made 20.5.2013, z[2] changed to z[1]
 rhs[7]=Inoutz[10]+Inoutz[Inoutn+10];
 for (i=1 ; i<= Inoutn-1 ; i++ )
  rhs[7]=rhs[7]+Inoutz[10+i]+Inoutz[10+i];

 rhs[8]=Inoutz[10]-Inoutz[Inoutn+10]-Inoutz[2];

 for ( m=0 ; m <= Inoutn ; m++ )
 {
  psi=0.;
  u=0.;
  v=0.;

  for (j=1 ; j <= Inoutn ; j++ )
  {
   nm = (m*j) % (Inoutn+Inoutn);
   e=exp(j*(Inoutz[10+m]));
   s=0.5*(e-1./e);
   c=0.5*(e+1./e);
   psi=psi+coeff[j]*(s+c*FourierTanh[j])*cosa[nm];
   u=u+j*coeff[j]*(c+s*FourierTanh[j])*cosa[nm];
   v=v+j*coeff[j]*(s+c*FourierTanh[j])*sina[nm];
  }
  rhs[m+9]=psi-Inoutz[8]-Inoutz[7]*Inoutz[m+10];
  rhs[Inoutn+m+10]=0.5*(pow((-Inoutz[7]+u),2.)+v*v)+Inoutz[m+10]-Inoutz[9];
 }

 for (j=1, s=0. ; j <= num ; j++ ) s += rhs[j]*rhs[j];
 return s;
}

// **************************************************
//	SET UP JACOBIAN MATRIX AND SOLVE MATRIX EQUATION
// **************************************************

double FourierClass::Newton(int Inoutn,char* InoutCase,double* Inoutz,double InoutCurrent,int InoutCurrent_criterion)
{
 double	**a, *rhs, *x;
 double 	h, sum;

 int i, j;

 Eqns(rhs1,Inoutn,InoutCase,Inoutz,InoutCurrent,InoutCurrent_criterion);

 rhs = new double[num+1];
 x = new double[num+1];
 a = new double* [num+1];
 for(i=0;i<num+1;i++)
  a[i] = new double [num+1];

 for ( i=1 ; i<=num ; i++ )
 {
  h=0.01*Inoutz[i];
  if(fabs(Inoutz[i]) < 1.e-4) h = 1.e-5;
  Inoutz[i]=Inoutz[i]+h;
  Eqns(rhs2,Inoutn,InoutCase,Inoutz,InoutCurrent,InoutCurrent_criterion);
  Inoutz[i]=Inoutz[i]-h;
  rhs[i] = -rhs1[i];
  for ( j=1 ; j<=num ; j++ )
   a[j][i] = (rhs2[j] - rhs1[j])/h;
 }

// **************************************************
// SOLVE MATRIX EQUATION
// **************************************************

 FentonSolve(a, rhs, num, num, x);

 for ( i=1 ; i<=num ; i++ )
  Inoutz[i] += x[i];

 for ( sum = 0., i=10 ; i<= Inoutn+10 ; i++ )
  sum += fabs(x[i]);

 sum /= Inoutn;

 delete [] rhs;
 delete [] x;

 for(i=0;i<num+1;i++)
  delete [] a[i];
 delete [] a;

 return(sum);
}

#endif

