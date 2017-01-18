
#ifndef FentonInout_H
#define FentonInout_H

#include <math.h>
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define Skip		fgets(dummy,100,DataInput)
#define pi			3.14159265358979324
#define iff(x,y)	if(strcmp(x,#y)==0)
#define HI			"\t%10.6f"
#define LO			"\t%7.4f"
#define Readtext(x)		fgets(x,400,DataInput); x[strlen(x)-1] = '\0'
#define Read(x,y)		{fscanf(DataInput,"%"#y, &x);Skip;}

template<typename SolutionType>
class FentonInoutClass
{
 private:
  FILE *monitor,*DataInput,*ConvergenceInput,*PointsInput,*ElevationOutput,*FlowfieldOutput,*SolutionOutput;
  char Data_file[50],Convergence_file[50],Points_file[50]; // for input
  char Elevation_file[50],Flowfield_file[50],Solution_file[50]; // for output

  // Convergence_file for FourierClass only

  char Extra;

  FILE *DiagnosticOutput;
  char Diagnostic_file[50];

  char Title[100],dummy[100],Case[20],Currentname[10],Current1[10],Current2[10],Method[100];
  int Current_criterion,n,nstep,number,Nprofiles,Surface_points,Points;

  // nstep: Number of height steps to reach H/d
  // number: Maximum number of iterations for each height step
  // Surface_points: Number of points on free surface
  // Nprofiles: Number of velocity profiles over half a wavelength to print out
  // Points: Number of vertical points in each profile

  double *B,*Y,*z;
// see page 359 in Fenton (1988) for the definition of z and B
// see Eq. (3.5) in Fenton (1999) for the definition of B
// see page 363 in Fenton (1988) for the definition of Y

  double Bernoulli_check,c,ce,crit,criter,cs,Current,dphidt,dudt,dvdt,Eta,f,H,Height,Highest,SU,Hoverd,height,kd,kH,ke,L,
         MaxH,Mean_half_keta_squared,pe,phi,Pressure,pulse,Q,q,R,r,s,sum,sxx,T,U,u,ub2,ubar,ut,ux,uy,v,vt,vx,vy;

  double Surface(double x);
  void Point(double X, double y);
  void Title_block(FILE*);
  int Read_data();
 public:
  FentonInoutClass(char* data_file,char* convergence_file,char* points_file);
  FentonInoutClass(char* data_file,char* points_file);
  ~FentonInoutClass();
  void Solve();
// read data (n, Case, T, L, MaxH, Current_criterion, Current) from Inout
// and set Method, z[], B[], Y[], H, kd, kH, and SU to Inout
  void Output(char* elevation_file,char* flowfield_file,char* solution_file,char* diagnostic_file);

  friend class StokesClass;
  friend class FourierClass;

};

template<typename SolutionType>
FentonInoutClass<SolutionType>::FentonInoutClass(char* data_file,char* convergence_file,char* points_file)
{
 Extra='y';

 strcpy(Current1,"Euler");
 strcpy(Current2,"Stokes");
 strcpy(Data_file,data_file);
 strcpy(Convergence_file,convergence_file);
 strcpy(Points_file,points_file);

 monitor = stdout;
 assert(Read_data()==1);
}

template<typename SolutionType>
FentonInoutClass<SolutionType>::FentonInoutClass(char* data_file,char* points_file)
{
 Extra='n';

 strcpy(Current1,"Euler");
 strcpy(Current2,"Stokes");
 strcpy(Data_file,data_file);
 strcpy(Points_file,points_file);

 monitor = stdout;
 assert(Read_data()==1);
}

template<typename SolutionType>
FentonInoutClass<SolutionType>::~FentonInoutClass()
{
 delete [] z;
 delete [] Y;
 delete [] B;
}

template<typename SolutionType>
int FentonInoutClass<SolutionType>::Read_data()
{
 DataInput = fopen(Data_file,"r");
 Readtext(Title);
 iff(Title, FINISH) return(0);
 fprintf(monitor,"# %s", Title);

 Read(MaxH,lf);
 fprintf(monitor,"\n\n# Height/Depth:%6.3f", MaxH);
 fscanf(DataInput,"%s", Case); Skip;
 iff(Case,Wavelength)
 {
  Read(L,lf);
  fprintf(monitor,"  Length/Depth:%7.2f", L);
  Height = MaxH/L;
 }
 iff(Case,Period)
 {
  Read(T,lf);
  Height = MaxH/(T*T);
  fprintf(monitor,"  Dimensionless Period T*sqrt(g/d):%7.2f", T);
 }
 Read(Current_criterion,d);
 Read(Current,lf);
 if(Current_criterion == 1) strcpy(Currentname, Current1);
 if(Current_criterion == 2) strcpy(Currentname, Current2);
 fprintf(monitor,"\n# Current criterion: %s,  Dimensionless value:%6.3lf", Currentname, Current);

 Read(n,d);
 Read(nstep,d);

// If wavelength is known at this stage the program calculates the highest wave,
// sees how close this wave is and allocates the number of steps automatically.

 iff(Case,Wavelength)
 {
  Highest = (0.0077829*L*L*L+0.0095721*L*L+0.141063*L)
  /(0.0093407*L*L*L+0.0317567*L*L+0.078834*L+1);
  printf("\n# MaxH/Highest: %f of the maximum of H/d = %f",MaxH/Highest,Highest);
  //nstep = 1+20*pow(MaxH/Highest,2);
 }

 fclose(DataInput);

 if(Extra=='y')
 {

// Convergence criteria for FourierClass
  ConvergenceInput = fopen(Convergence_file,"r");
  fgets(dummy,400,ConvergenceInput);
  fscanf(ConvergenceInput,"%d", &number);fgets(dummy,400,ConvergenceInput);
  fscanf(ConvergenceInput,"%le", &crit);fgets(dummy,400,ConvergenceInput);
  fclose(ConvergenceInput);

 }

// Number of data points to present results for

 PointsInput = fopen(Points_file,"r");
 fgets(dummy,400,PointsInput);
// Number of points on surface profile (clustered quadratically near crest)
 fscanf(PointsInput,"%d", &Surface_points);fgets(dummy,400,PointsInput);
// Number of vertical profiles
 fscanf(PointsInput,"%d", &Nprofiles);fgets(dummy,400,PointsInput);
// Number of points in each profile
 fscanf(PointsInput,"%d", &Points);fgets(dummy,400,PointsInput);

 fclose(PointsInput);

 return(1);

}

//	PRINT OUT TITLE BLOCK

template<typename SolutionType>
void FentonInoutClass<SolutionType>::Title_block(FILE* file)
{
 fprintf(file,"# %s", Title);
 fprintf(file,"\n\n# Height/Depth:%6.3f", z[2]
         / z[1]);
 fprintf(file,"  Length/Depth:%7.2f", 2*pi/z[1]);
 fprintf(file,"  Dimensionless Period T*sqrt(g/d):%7.2f", z[3]/sqrt(z[1]));
 // Highest wave - eqn (32) of Fenton (1990)
 L = 2*pi/z[1];
 Highest = (0.0077829*L*L*L+0.0095721*L*L+0.141063*L)
   		  /(0.0093407*L*L*L+0.0317567*L*L+0.078834*L+1);
 fprintf(file,"\n# (A height of %3.0lf\%% of the maximum of H/d =%6.3f for this length)", z[2]/z[1]/Highest*100., Highest);
 fprintf(file,"\n# Current criterion: %s,  Dimensionless value:%6.3lf\n", Currentname, Current);
 fprintf(file,"\n# Solution by %s", Method);
}

// read data (n, Case, T, L, MaxH, Current_criterion, Current) from *this
// and set Method, z[], B[], Y[], H, kd, kH, and SU to *this
template<typename SolutionType>
void FentonInoutClass<SolutionType>::Solve()
{
 SolutionType SolutionObject(*this);
}

template<typename SolutionType>
void FentonInoutClass<SolutionType>::Output(char* elevation_file,char* flowfield_file,char* solution_file,char* diagnostic_file)
{

 strcpy(Elevation_file,elevation_file);
 strcpy(Flowfield_file,flowfield_file);
 strcpy(Solution_file,solution_file);
 strcpy(Diagnostic_file,diagnostic_file);

 int 		i, I;
 double 	X, eta, y;

 SolutionOutput = fopen(Solution_file,"w");
 Title_block(SolutionOutput);

 kd = z[1];
 L=2*pi/z[1];
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
 q=z[7]*z[1]-z[8];
 r=z[9]+z[1];
 s=sxx-2.*z[4]*pulse+(z[4]*z[4]+0.5*z[1])*z[1];

 fprintf(SolutionOutput, "\n\n# Stokes-Ursell parameter (SU): %7.3f", 0.5*z[2]/pow(z[1],3));
 fprintf(SolutionOutput, "\n\n# Integral quantities - notation from Fenton (1988)");
 fprintf(SolutionOutput, "\n# (1) Quantity, (2) symbol, solution non-dimensionalised by (3) g & wavenumber, and (4) g & mean depth\n");
 fprintf(SolutionOutput, "\n# Water depth                        (d)" LO LO, z[1], 1.);
 fprintf(SolutionOutput, "\n# Wave length                   (lambda)" LO LO, 2*pi, L);
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

 fclose(SolutionOutput);

//	Surface - print out coordinates of points on surface for plotting plus check of pressure on surface

 ElevationOutput = fopen(Elevation_file,"w");
 fprintf(ElevationOutput,  "# %s\n", Title);
 fprintf(ElevationOutput,  "\n# Solution by %s\n", Method);
 fprintf(ElevationOutput,  "\n# Surface of wave - trough-crest-trough,");
 fprintf(ElevationOutput,  " note quadratic point spacing clustered around crest");
 fprintf(ElevationOutput,  "\n# Non-dimensionalised with respect to depth");
 fprintf(ElevationOutput,  "\n# X/d, eta/d, & check of surface pressure\n");

 for ( i=-Surface_points/2 ; i <= Surface_points/2; i++)
 {
  X = 2 * L * (i * fabs(i)/Surface_points/Surface_points);	//NB Quadratic point spacing, clustered near crest
  eta = Surface(X);
  Point(X,eta);
  fprintf(ElevationOutput,  "\n%8.4lf\t%7.4f\t%7.0e", X, eta, Pressure);
 }
 fprintf(ElevationOutput,  "\n\n");
 fclose(ElevationOutput);

//	fprintf(out,   out Velocity and acceleration profiles plus check of Bernoulli


 FlowfieldOutput = fopen(Flowfield_file,"w");
 fprintf(FlowfieldOutput,  "# %s\n", Title);
 fprintf(FlowfieldOutput,  "\n# Solution by %s\n", Method);
 fprintf(FlowfieldOutput,  "\n# Velocity and acceleration profiles and Bernoulli checks\n");
 fprintf(FlowfieldOutput,  "\n# All quantities are dimensionless with respect to g and/or d\n");
 fprintf(FlowfieldOutput,  "\n#*******************************************************************************");
 fprintf(FlowfieldOutput,  "\n# y        u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check  ");
 fprintf(FlowfieldOutput,  "\n# -     -------------   -------  ------   -----  ------------- ---------------  ");
 fprintf(FlowfieldOutput,  "\n# d        sqrt(gd)       gd        g       g       sqrt(g/d)        gd         ");
 fprintf(FlowfieldOutput,  "\n#*******************************************************************************");

 for(I = 0; I <= Nprofiles ; ++I)
 {
  X = 0.5 * L * I/(Nprofiles);
  eta = Surface(X);
  fprintf(FlowfieldOutput,  "\n\n# X/d = %8.4f, Phase = %6.1f°\n", X, X/L*360);

  for(i=0 ; i <= Points; ++i)
  {
   y = (i)*eta/(Points);
   Point(X, y);
   fprintf(FlowfieldOutput,  "\n%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f",
        y, u, v, dphidt, ut, vt, ux, uy, Bernoulli_check);
  }
 }
 fprintf(FlowfieldOutput,  "\n\n");
 fclose(FlowfieldOutput);

 double Velo[Points+1], sum1, sum2, ucm;
 X = 0.;
 eta = Surface(X);

 for(i=0 ; i <= Points; ++i)
 {
  y = (i)*eta/(Points);
  Point(X, y);
  Velo[i] = u;
 }

 for(i=1, sum1=0; i <= Points-1; i+=2) sum1 += Velo[i];
 for(i=2, sum2=0; i <= Points-2; i+=2) sum2 += Velo[i];
 ucm = (Velo[0]+4*sum1+2*sum2+Velo[Points])/3./Points;

 DiagnosticOutput = fopen(Diagnostic_file,"a");
 fprintf(DiagnosticOutput,"\n%s %2d\t%7.4f\t%8.3f\t%7.3f\t%3.0f\t%8.5f",
			Method, n, H, L, 0.5*z[2]/pow(z[1],3), z[2]/z[1]/Highest*100., ucm);
 fclose(DiagnosticOutput);
}

// Surface elevation

template<typename SolutionType>
double FentonInoutClass<SolutionType>::Surface(double x)
{
 int j;
 static double kEta;

 kEta = kd;
 for ( j = 1 ; j < n ; j++ )
	kEta += Y[j] * cos(j*x*kd);
 kEta += 0.5*Y[n] * cos(n*x*kd);
 return (kEta/kd);
}


// Velocities, accelerations, and pressure at a point

template<typename SolutionType>
void FentonInoutClass<SolutionType>::Point(double X, double y)
{
 int j;
 double Cosh, Sinh, Sin, Cos;

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

// All PHI, PSI, u, v, ux and vx are dimensionless w.r.t. g & k.
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

