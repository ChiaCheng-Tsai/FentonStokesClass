
#ifndef FourierStokesData_H
#define FourierStokesData_H
// Stokes theory calculations

#include <math.h>
#include <stdio.h>
#include <process.h>
#include <string.h>
#include <conio.h>
#include <stdlib.h>
#include <assert.h>
#define	ANSI

class FourierStokesData
{
 private:
  // FentonInout input variables
  double Current,H,T,L;
  int Current_criterion,n;
  char Title[100],Case[20],Method[100];

  // Class internal variables
  double kd,kh,kH,SU;

  double *B,*Y,*z;
// see page 359 in Fenton (1988) for the definition of z and B
// see Eq. (3.5) in Fenton (1999) for the definition of B
// see page 363 in Fenton (1988) for the definition of Y
 public:
  FourierStokesData(){};
  ~FourierStokesData();
  void Output(FILE *SolutionOutput);
  double Surface(double x);
  void Point(double X,double y,
             double &phi,double &u,double &v,double &dphidt,double &dudt,double &dvdt,double &ut,double &vt,
             double &ux,double &vx,double &uy,double &vy,double &Pressure,double &Bernoulli_check);

  friend class StokesClass;
  friend class FourierClass;
};

FourierStokesData::~FourierStokesData()
{
 delete [] z;
 delete [] Y;
 delete [] B;
}

void FourierStokesData::Output(FILE *SolutionOutput)
{
 int i;
 double c,ce,cs,ubar,Q,R,pulse,ke,pe,ub2,sxx,f,s;

 kd = z[1];
 L=2*Fenton_pi/z[1];
 H=z[2]/z[1];
 T=z[3]/sqrt(z[1]);
 c=z[4]/sqrt(z[1]);
 ce=z[5]/sqrt(z[1]);
 cs=z[6]/sqrt(z[1]);
 ubar=z[7]/sqrt(z[1]);
 Q=ubar-z[8]/pow(kd,1.5);
 R=1+z[9]/z[1];

 pulse=z[8]+z[1]*z[5];
 ke=0.5*(z[4]*pulse-z[5]*Q*pow(kd,1.5));

 // Calculate potential energy, not by computing the mean of 1/2 (eta-d)^2
 // but by exploiting orthogonality of the cosine functions to give the sum of 1/4 Y[i]^2
 pe = 0;
 for(i=1;i<=n;++i)
  pe += 0.25*pow(Y[i],2);

 ub2=2.*z[9]-z[4]*z[4];
 sxx=4.*ke-3.*pe+ub2*z[1]+2.*z[5]*(z[7]*z[1]-z[8]);
 f=z[4]*(3.*ke-2.*pe)+0.5*ub2*(pulse+z[4]*z[1])+z[4]*z[5]*(z[7]*z[1]-z[8]);
 //q=z[7]*z[1]-z[8];
 //r=z[9]+z[1];
 s=sxx-2.*z[4]*pulse+(z[4]*z[4]+0.5*z[1])*z[1];

 fprintf(SolutionOutput, "\n\n# Stokes-Ursell parameter (SU): %7.3f", 0.5*z[2]/pow(z[1],3));
 fprintf(SolutionOutput, "\n\n# Integral quantities - notation from Fenton (1988)");
 fprintf(SolutionOutput, "\n# (1) Quantity, (2) symbol, solution non-dimensionalised by (3) g & wavenumber, and (4) g & mean depth\n");
 fprintf(SolutionOutput, "\n# Water depth                        (d)" LO LO, z[1], 1.);
 fprintf(SolutionOutput, "\n# Wave length                   (lambda)" LO LO, 2*Fenton_pi, L);
 fprintf(SolutionOutput, "\n# Wave height                        (H)" LO LO, z[2], H);
 fprintf(SolutionOutput, "\n# Wave period                      (tau)" LO LO, z[3], T);
 fprintf(SolutionOutput, "\n# Wave speed                         (c)" LO LO, z[4], c);
 fprintf(SolutionOutput, "\n# Eulerian current                 (u1_)" LO LO, z[5], ce);
 fprintf(SolutionOutput, "\n# Stokes current                   (u2_)" LO LO, z[6], cs);
 fprintf(SolutionOutput, "\n# Mean fluid speed in frame of wave (U_)" LO LO, z[7], ubar);
 fprintf(SolutionOutput, "\n# Volume flux due to waves           (q)" LO LO, z[8], z[8]/pow(kd,1.5));
 fprintf(SolutionOutput, "\n# Bernoulli constant                 (r)" LO LO, z[9], z[9]/kd);
 fprintf(SolutionOutput, "\n# Volume flux                        (Q)" LO LO, Q*pow(kd,1.5), Q);
 fprintf(SolutionOutput, "\n# Bernoulli constant                 (R)" LO LO, R*kd, R);
 fprintf(SolutionOutput, "\n# Momentum flux                      (S)" LO LO, s, s/kd/kd );
 fprintf(SolutionOutput, "\n# Impulse                            (I)" LO LO, pulse, pulse/pow(kd,1.5));
 fprintf(SolutionOutput, "\n# Kinetic energy                     (T)" LO LO, ke, ke/kd/kd);
 fprintf(SolutionOutput, "\n# Potential energy                   (V)" LO LO, pe, pe/kd/kd);
 fprintf(SolutionOutput, "\n# Mean square of bed velocity     (ub2_)" LO LO, ub2, ub2/kd);
 fprintf(SolutionOutput, "\n# Radiation stress                 (Sxx)" LO LO, sxx, sxx/kd/kd);
 fprintf(SolutionOutput, "\n# Wave power                         (F)" LO LO, f, f/pow(kd,2.5));

 fprintf(SolutionOutput, "\n\n# Dimensionless coefficients in Fourier series" );
 fprintf(SolutionOutput, "\n# Potential/Streamfn\tSurface elevations" );
 fprintf(SolutionOutput, "\n# j, B[j], & E[j], j=1..N\n" );

 for ( i=1 ; i <= n ; i++ )
  fprintf(SolutionOutput, "\n%2d\t%15.7e\t%15.7e", i, B[i], Y[i]);
 fprintf(SolutionOutput, "\n\n" );

 return;
}

// Surface elevation

double FourierStokesData::Surface(double x)
{
 int j;
 double kEta;

 kEta = kd;
 for ( j = 1 ; j < n ; j++ )
	kEta += Y[j] * cos(j*x*kd);
 kEta += 0.5*Y[n] * cos(n*x*kd);
 return (kEta/kd);
}


// Velocities, accelerations, and pressure at a point

void FourierStokesData::Point(double X,double y,
  double &phi,double &u,double &v,double &dphidt,double &dudt,double &dvdt,double &ut,double &vt,
  double &ux,double &vx,double &uy,double &vy,double &Pressure,double &Bernoulli_check)
{

 int j;
 double Cosh, Sinh, Sin, Cos;
 double c,ce,R;

 c=z[4]/sqrt(z[1]);
 ce=z[5]/sqrt(z[1]);
 R=1+z[9]/z[1];

 u = v = ux = vx = phi = 0.;

 for ( j = 1 ; j <= n ; j++ )
 {
  Cos  = cos(j*X*kd);
  Sin  = sin(j*X*kd);

  if(y<0.5)
  {
   Cosh = cosh(j*kd*y)/cosh(j*kd);
   Sinh = sinh(j*kd*y)/cosh(j*kd);
  }
  else
  {
   Cosh = cosh(j*kd*(y-1.))+sinh(j*kd*(y-1.))*tanh(j*kd);
   Sinh = sinh(j*kd*(y-1.))+cosh(j*kd*(y-1.))*tanh(j*kd);
  }

  phi += B[j] * Cosh * Sin;
  u += j * B[j] * Cosh * Cos;
  v += j * B[j] * Sinh * Sin;
  ux += - j * j * B[j] * Cosh * Sin;
  vx += j * j * B[j] * Sinh * Cos;
 }

// All PHI, u, v, ux and vx are dimensionless w.r.t. g & k.
// Now convert to dimensionless w.r.t. d.

 phi /= pow(kd,1.5);
 u /= pow(kd,0.5);
 v /= pow(kd,0.5);
 ux *= pow(kd,0.5);
 vx *= pow(kd,0.5);

 u = ce + u;
 phi = ce * X + phi;
 dphidt = -c * u;

 ut = -c * ux;
 vt = -c * vx;
 uy = vx;
 vy = -ux;
 dudt = ut + u*ux + v*uy;
 dvdt = vt + u*vx + v*vy;
 Pressure = R - y - 0.5 * ((u-c)*(u-c)+v*v);
 Bernoulli_check = dphidt + Pressure + y + 0.5*(u*u+v*v)-(R-0.5*c*c);

 return;
}

#endif

