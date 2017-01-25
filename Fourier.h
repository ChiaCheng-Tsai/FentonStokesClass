
#ifndef Fourier_H
#define Fourier_H
// Fourier theory calculations

#include <math.h>
#include <stdio.h>
#include <process.h>
#include <string.h>
#include <conio.h>
#include <stdlib.h>
#include <assert.h>
#define	ANSI

#include "FentonHeaders.h"
#include "FentonInout.h"
#include "FentonSolve.h"
#include "FourierStokesData.h"

class FourierClass
{
 private:
  int num;
  double height,Hoverd;
  double criter;
  double *coeff,*cosa,*rhs1,*rhs2,*sina,*FourierTanh;
  double **sol;
  FourierStokesData FSD;
 public:
  FourierClass(){};
  FourierClass(FentonInout<FourierClass> & );
  ~FourierClass();
  void init(int Inoutn,char* InoutCase,double* Inoutz,double InoutCurrent,int InoutCurrent_criterion);
  double Eqns(double *rhs,int Inoutn,char* InoutCase,double* Inoutz,double InoutCurrent,int InoutCurrent_criterion);
  double Newton(int Inoutn,char* InoutCase,double* Inoutz,double InoutCurrent,int InoutCurrent_criterion);

  double Get_kd(){return FSD.kd;};
  double Get_kh(){return FSD.kh;};
  double Get_kH(){return FSD.kH;};
  double Get_SU(){return FSD.SU;};
  double Get_L(){return FSD.L;};
  double Get_H(){return FSD.H;};
  int Get_n(){return FSD.n;};
  double Get_z(int i){return FSD.z[i];};
  void Output(FILE *Solution){FSD.Output(Solution);};
  double Surface(double x){return FSD.Surface(x);};
  void Point(double X,double y,
             double &phi,double &u,double &v,double &dphidt,double &dudt,double &dvdt,double &ut,double &vt,
             double &ux,double &vx,double &uy,double &vy,double &Pressure,double &Bernoulli_check)
             {FSD.Point(X,y,phi,u,v,dphidt,dudt,dvdt,ut,vt,ux,vx,uy,vy,Pressure,Bernoulli_check);};

  // dummy member function
  double Get_h(){assert(false);return 0;};
  double Get_R_d(){assert(false);return 0;};
  double Get_c(){assert(false);return 0;};
  double eta_h(double x){assert(false);return 0;};
  double u_h(double x, double Y){assert(false);return 0;};
  double v_h(double x, double Y){assert(false);return 0;};

};

// read data (n, Case, T, L, MaxH, Current_criterion, Current) from Inout
// and set Method, z[], B[], Y[], H, kd, kH, and SU to Inout

// Height, number, crit, and nstep are additionally needed compared with StokesClass

FourierClass::FourierClass(FentonInout<FourierClass> & Inout)
{
 int	i, j, iter, m, ns;
 double	dhe, dho, error;
 double sum;

 FSD.Current=Inout.Current;
 FSD.Current_criterion=Inout.Current_criterion;
 FSD.H=Inout.H;
 FSD.T=Inout.T;
 FSD.L=Inout.L;
 FSD.n=Inout.n;

 strcpy(FSD.Case,Inout.Case);
 strcpy(FSD.Title,Inout.Title);

 sprintf(Inout.Method, "Fourier method with %d terms in series", FSD.n);
 sprintf(FSD.Method, "Fourier method with %d terms in series", FSD.n);
 printf( "\n# Solution by Fourier method with %d terms in series", FSD.n);

 num=2*FSD.n+10;
 dhe=Inout.Height/Inout.nstep;
 dho=Inout.MaxH/Inout.nstep;

 FSD.z = new double[num+1];
 FSD.Y = new double[FSD.n+1]; // bug fixed by Tsai
 FSD.B = new double[FSD.n+1];

// see page 359 in Fenton (1988) for the definition of z and B
// see Eq. (3.5) in Fenton (1999) for the definition of B
// see page 363 in Fenton (1988) for the definition of Y

 rhs1 = new double[num+1];
 rhs2 = new double[num+1];
 coeff = new double[num+1];
 cosa = new double[2*FSD.n+1];
 sina = new double[2*FSD.n+1];
 FourierTanh = new double[FSD.n+1];

 sol = new double* [num+1];
 for(i=0;i<num+1;i++)
  sol[i] = new double [3];

//	Commence stepping through steps in wave height
 for ( ns = 1 ; ns <= Inout.nstep ; ns++ )
 {
  height=ns*dhe;
  Hoverd=ns*dho;
  printf("\n\nHeight step %2d of %2d", ns, Inout.nstep);

//	Calculate initial linear solution
  if(ns <= 1)
    init(FSD.n,FSD.Case,FSD.z,FSD.Current,FSD.Current_criterion);
//	Or, extrapolate for next wave height, if necessary

  else
   for ( i=1 ; i <= num ; i++ )
    FSD.z[i]=2.*sol[i][2]-sol[i][1];

//	Commence iterative solution
  for (iter=1 ; iter <= Inout.number ; iter++ )
  {
   printf("\nIteration%3d:", iter);

//	Calculate right sides of equations and differentiate numerically
//	to obtain Jacobian matrix, then solve matrix equation

   error = Newton(FSD.n,FSD.Case,FSD.z,FSD.Current,FSD.Current_criterion);

//	Convergence criterion satisfied?

   printf(" Mean of corrections to free surface: %8.1e", error);

   if(ns == Inout.nstep)	criter = 1.e-10 ;
   else			criter = Inout.crit;

   if((error < criter * fabs(FSD.z[1]))  && iter > 1 ) break;
   if(iter == Inout.number)
   {
    fprintf(stdout,"\nNote that the program still had not converged to the degree specified\n");
   }

//	Operations for extrapolations if more than one height step used

   if(ns == 1)
    for ( i=1 ; i<=num ; i++ )
     sol[i][2] = FSD.z[i];
   else
    for ( i=1 ; i<=num ; i++ )
    {
     sol[i][1] = sol[i][2];
     sol[i][2] = FSD.z[i];
    }
  }

//	Fourier coefficients (for surface elevation by slow Fourier transform)

  for ( FSD.Y[0] = 0., j = 1 ; j <= FSD.n ; j++ )
  {
   FSD.B[j]=FSD.z[j+FSD.n+10];
   sum = 0.5*(FSD.z[10]+FSD.z[FSD.n+10]*pow(-1.,(double)j));
   for ( m = 1 ; m <= FSD.n-1 ; m++ )
    sum += FSD.z[10+m]*cosa[(m*j)%(FSD.n+FSD.n)];
   FSD.Y[j] = 2. * sum / FSD.n;
  }
 } // End stepping through wave heights

 FSD.L = 2*Fenton_pi/FSD.z[1];
 FSD.H = FSD.z[2]/FSD.z[1];

}

// **************************************************
// CALCULATE INITIAL SOLUTION FROM LINEAR WAVE THEORY
// **************************************************

FourierClass::~FourierClass()
{
 delete [] rhs1;
 delete [] rhs2;
 delete [] coeff;
 delete [] cosa;
 delete [] sina;
 delete [] FourierTanh;

 for(int i=0;i<num+1;i++)
  delete [] sol[i];
 delete [] sol;
}

void FourierClass::init(int Inoutn,char* InoutCase,double* Inoutz,double InoutCurrent,int InoutCurrent_criterion)
{
 int i;
 double a, b, t;

 iff(InoutCase,Period)
 {
  a=4.*Fenton_pi*Fenton_pi*height/Hoverd;
  b=a/sqrt(tanh(a));
  t=tanh(b);
  Inoutz[1]=b+(a-b*t)/(t+b*(1.-t*t));
 }
 else
  Inoutz[1]=2.*Fenton_pi*height/Hoverd;

 Inoutz[2]=Inoutz[1]*Hoverd;
 Inoutz[4]=sqrt(tanh(Inoutz[1]));
 Inoutz[3]=2.*Fenton_pi/Inoutz[4];

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
  cosa[i]=cos(i*Fenton_pi/Inoutn);
  cosa[i+Inoutn]=cos((i+Inoutn)*Fenton_pi/Inoutn);
  sina[i]=sin(i*Fenton_pi/Inoutn);
  sina[i+Inoutn]=sin((i+Inoutn)*Fenton_pi/Inoutn);
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
  rhs[2]=Inoutz[2]-2.*Fenton_pi*height;
 else
  rhs[2]=Inoutz[2]-height*Inoutz[3]*Inoutz[3];

 rhs[3]=Inoutz[4]*Inoutz[3]-Fenton_pi-Fenton_pi;
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

