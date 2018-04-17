
#ifndef FentonInout_H
#define FentonInout_H

#include <math.h>
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "FentonHeaders.h"
#include "FentonElliptic.h"

enum { FourierSolution = 0, StokesSolution = 1, CnoidalSolution = 2};

template<typename SolutionClass>
class FentonInout
{
 private:
  int ModelType;
  FILE *monitor,*DataInput,*ConvergenceInput,*PointsInput,*SurfaceOutput,*FlowfieldOutput,*SolutionOutput;
  char Data_file[50],Convergence_file[50],Points_file[50]; // for input
  char Surface_file[50],Flowfield_file[50],Solution_file[50]; // for output

  // Convergence_file for FourierClass only

  char Extra;

  FILE *AccuracyOutput;
  char Accuracy_file[50];

  char dummy[100],Currentname[10],Current1[10],Current2[10];

  // for Read_data() from DataInput
  char Title[100],Case[20],Method[100];
  int Current_criterion,n,nstep;
  double H,L,T,Current;

  int number,Nprofiles,Surface_points,Points;

  // nstep: Number of height steps to reach H/d
  // number: Maximum number of iterations for each height step
  // Surface_points: Number of points on free surface
  // Nprofiles: Number of velocity profiles over half a wavelength to print out
  // Points: Number of vertical points in each profile

  double MaxH,Height,crit;

  SolutionClass SolutionObject;

  void Title_block(FILE*);
  int Read_data();
  double Get_Highest(double LL){return (0.0077829*LL*LL*LL+0.0095721*LL*LL+0.141063*LL)
    		  /(0.0093407*LL*LL*LL+0.0317567*LL*LL+0.078834*LL+1);};

 public:
  FentonInout(int modeltype,char* data_file,char* convergence_file,char* points_file);
  FentonInout(int modeltype,char* data_file,char* points_file);
  void Solve();
// read data (n, Case, T, L, MaxH, Current_criterion, Current) from Inout
  void Output(char* elevation_file,char* flowfield_file,char* solution_file,char* diagnostic_file);

  friend class StokesClass;
  friend class FourierClass;
  friend class CnoidalClass;

};

template<typename SolutionClass>
FentonInout<SolutionClass>::FentonInout(int modeltype,char* data_file,char* convergence_file,char* points_file):ModelType(modeltype)
{
 Extra='y';

 strcpy(Current1,"Euler");
 strcpy(Current2,"Stokes");
 strcpy(Data_file,data_file);
 strcpy(Convergence_file,convergence_file);
 strcpy(Points_file,points_file);

 monitor = stdout;
 assert(Read_data()==1);
 Solve();
}

template<typename SolutionClass>
FentonInout<SolutionClass>::FentonInout(int modeltype,char* data_file,char* points_file):ModelType(modeltype)
{
 Extra='n';

 strcpy(Current1,"Euler");
 strcpy(Current2,"Stokes");
 strcpy(Data_file,data_file);
 strcpy(Points_file,points_file);

 monitor = stdout;
 assert(Read_data()==1);
 Solve();
}

template<typename SolutionClass>
int FentonInout<SolutionClass>::Read_data()
{
 if ( ModelType==FourierSolution || ModelType==StokesSolution )
 {
  DataInput = fopen(Data_file,"r");
  Readtext(DataInput,Title);
  iff(Title, FINISH) return(0);
  fprintf(monitor,"# %s", Title);

  Read(DataInput,MaxH,lf);
  fprintf(monitor,"\n\n# Height/Depth:%6.3f", MaxH);
  fscanf(DataInput,"%s", Case); Skip(DataInput);
  iff(Case,Wavelength)
  {
   Read(DataInput,L,lf);
   fprintf(monitor,"  Length/Depth:%7.2f", L);
   Height = MaxH/L;
  }
  iff(Case,Period)
  {
   Read(DataInput,T,lf);
   Height = MaxH/(T*T);
   fprintf(monitor,"  Dimensionless Period T*sqrt(g/d):%7.2f", T);
  }
  Read(DataInput,Current_criterion,d);
  Read(DataInput,Current,lf);
  if(Current_criterion == 1) strcpy(Currentname, Current1);
  if(Current_criterion == 2) strcpy(Currentname, Current2);
  fprintf(monitor,"\n# Current criterion: %s,  Dimensionless value:%6.3lf", Currentname, Current);

  Read(DataInput,n,d);
  Read(DataInput,nstep,d);

// If wavelength is known at this stage the program calculates the highest wave,
// sees how close this wave is and allocates the number of steps automatically.
  double HHighest;

  iff(Case,Wavelength)
  {
   HHighest = Get_Highest(L);
   printf("\n# MaxH/Highest: %f of the maximum of H/d = %f",MaxH/HHighest,HHighest);
  //nstep = 1+20*pow(MaxH/Highest,2);
  }

  fclose(DataInput);
 }
 else if ( ModelType==CnoidalSolution )
 {
  DataInput = fopen(Data_file,"r");
  Readtext(DataInput,Title);
  Read(DataInput,H,lf);
  fscanf(DataInput,"%s", Case); Skip(DataInput);

  iff(Case,Wavelength)
	{
	 Read(DataInput,L,lf);
	 if(L<10.) {printf("The wavelength is less than 10. Cnoidal theory should not be applied");}
	}

  iff(Case,Period)
	{
	 Read(DataInput,T,lf);
	 if(T<10.) {printf("The period is less than 10. Cnoidal theory should not be applied"); }
	}
  Read(DataInput,Current_criterion,d);
  Read(DataInput,Current,lf);
  Read(DataInput,n,d);
  if(n>6) n=6; //This was 6, here for purposes of testing the straight theory.

  fclose(DataInput);
 }
 else
  assert(false);

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

 Readtext(PointsInput,dummy);
 // Number of points on surface profile (clustered quadratically near crest)
 Read(PointsInput,Surface_points,d);
 // Number of vertical profiles
 Read(PointsInput,Nprofiles,d);
 // Number of points in each profile
 Read(PointsInput,Points,d);
 fclose(PointsInput);

 return(1);

}

//	PRINT OUT TITLE BLOCK

template<typename SolutionClass>
void FentonInout<SolutionClass>::Title_block(FILE* file)
{
 if ( ModelType==FourierSolution || ModelType==StokesSolution )
 {
  fprintf(file,"# %s", Title);
  fprintf(file,"\n\n# Height/Depth:%6.3f", SolutionObject.Get_z(2)/ SolutionObject.Get_z(1));
  fprintf(file,"  Length/Depth:%7.2f", 2*Fenton_pi/SolutionObject.Get_z(1));
  fprintf(file,"  Dimensionless Period T*sqrt(g/d):%7.2f", SolutionObject.Get_z(3)/sqrt(SolutionObject.Get_z(1)));
  // Highest wave - eqn (32) of Fenton (1990)
  double LL = 2*Fenton_pi/SolutionObject.Get_z(1);
  double HHighest = Get_Highest(LL);
  fprintf(file,"\n# (A height of %3.0lf\%% of the maximum of H/d =%6.3f for this length)", SolutionObject.Get_z(2)/SolutionObject.Get_z(1)/HHighest*100., HHighest);
  fprintf(file,"\n# Current criterion: %s,  Dimensionless value:%6.3lf\n", Currentname, Current);
  fprintf(file,"\n# Solution by %s", Method);
 }
 else if ( ModelType==CnoidalSolution )
 {
  fprintf(file,"%s", Title);

  fprintf(file,"\n\nHeight/Depth:%6.3f", H);
  iff(Case,Wavelength)
	fprintf(file,"  Length/Depth:%7.2f", L);
  iff(Case,Period)
	fprintf(file,"  Dimensionless Period:%7.2f", T);

  fprintf(file,"\nCurrent criterion: %d  Magnitude:%6.3lf", Current_criterion,Current);
  fprintf(file,"\n\nCnoidal theory");
 }
 else
 {
  assert(false);
 }

}

// read data (n, Case, T, L, MaxH, Current_criterion, Current) from *this
template<typename SolutionClass>
void FentonInout<SolutionClass>::Solve()
{
 if ( ModelType==CnoidalSolution )
  Title_block(monitor);

 SolutionObject=SolutionClass(*this);

}

template<typename SolutionClass>
void FentonInout<SolutionClass>::Output(char* elevation_file,char* flowfield_file,char* solution_file,char* diagnostic_file)
{

 strcpy(Surface_file,elevation_file);
 strcpy(Flowfield_file,flowfield_file);
 strcpy(Solution_file,solution_file);
 strcpy(Accuracy_file,diagnostic_file);

 if ( ModelType==FourierSolution || ModelType==StokesSolution )
 {

  int 		i, I;
  double 	X, eta, y;

  double HHighest,phi,u,v,dphidt,dudt,dvdt,ut,vt,ux,vx,uy,vy,Pressure,Bernoulli_check;

  double LL=SolutionObject.Get_L(),HH=SolutionObject.Get_H();

  SolutionOutput = fopen(Solution_file,"w");
  Title_block(SolutionOutput);
  SolutionObject.Output(SolutionOutput);
  fclose(SolutionOutput);

 //	Surface - print out coordinates of points on surface for plotting plus check of pressure on surface

  SurfaceOutput = fopen(Surface_file,"w");
  fprintf(SurfaceOutput,  "# %s\n", Title);
  fprintf(SurfaceOutput,  "\n# Solution by %s\n", Method);
  fprintf(SurfaceOutput,  "\n# Surface of wave - trough-crest-trough,");
  fprintf(SurfaceOutput,  " note quadratic point spacing clustered around crest");
  fprintf(SurfaceOutput,  "\n# Non-dimensionalised with respect to depth");
  fprintf(SurfaceOutput,  "\n# X/d, eta/d, & check of surface pressure\n");

  for ( i=-Surface_points/2 ; i <= Surface_points/2; i++)
  {
   X = 2 * LL * (i * fabs(i)/Surface_points/Surface_points);	//NB Quadratic point spacing, clustered near crest
   eta = SolutionObject.Surface(X);

   SolutionObject.Point(X,eta,phi,u,v,dphidt,dudt,dvdt,ut,vt,ux,vx,uy,vy,Pressure,Bernoulli_check);
   fprintf(SurfaceOutput,  "\n%8.4lf\t%7.4f\t%7.0e", X, eta, Pressure);
  }
  fprintf(SurfaceOutput,  "\n\n");
  fclose(SurfaceOutput);

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
   X = 0.5 * LL * I/(Nprofiles);
   eta = SolutionObject.Surface(X);
   fprintf(FlowfieldOutput,  "\n\n# X/d = %8.4f, Phase = %6.1f°\n", X, X/LL*360);

   for(i=0 ; i <= Points; ++i)
   {
    y = (i)*eta/(Points);
    SolutionObject.Point(X,y,phi,u,v,dphidt,dudt,dvdt,ut,vt,ux,vx,uy,vy,Pressure,Bernoulli_check);
    fprintf(FlowfieldOutput,  "\n%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f",
        y, u, v, dphidt, ut, vt, ux, uy, Bernoulli_check);
   }
  }
  fprintf(FlowfieldOutput,  "\n\n");
  fclose(FlowfieldOutput);

  double Velo[Points+1], sum1, sum2, ucm;
  X = 0.;
  eta = SolutionObject.Surface(X);

  for(i=0 ; i <= Points; ++i)
  {
   y = (i)*eta/(Points);
   SolutionObject.Point(X,y,phi,u,v,dphidt,dudt,dvdt,ut,vt,ux,vx,uy,vy,Pressure,Bernoulli_check);
   Velo[i] = u;
  }

  for(i=1, sum1=0; i <= Points-1; i+=2) sum1 += Velo[i];
  for(i=2, sum2=0; i <= Points-2; i+=2) sum2 += Velo[i];
  ucm = (Velo[0]+4*sum1+2*sum2+Velo[Points])/3./Points;

  HHighest = Get_Highest(LL);

  AccuracyOutput = fopen(Accuracy_file,"a");
  fprintf(AccuracyOutput,"\n%s %2d\t%7.4f\t%8.3f\t%7.3f\t%3.0f\t%8.5f",
			Method, n, HH, LL, 0.5*SolutionObject.Get_z(2)/pow(SolutionObject.Get_z(1),3), SolutionObject.Get_z(2)/SolutionObject.Get_z(1)/HHighest*100., ucm);
  fclose(AccuracyOutput);
 }
 else if ( ModelType==CnoidalSolution )
 {
  int 		i,ii;
  double 	x,y,eta,eta_d,uu,vv,HHighest;
  double SSU=SolutionObject.Get_SU(),cc=SolutionObject.Get_c(),HH=SolutionObject.Get_H(),hh=SolutionObject.Get_h(),RR_d=SolutionObject.Get_R_d(),LL=SolutionObject.Get_L();

  // Output solution summary
  SolutionOutput = fopen(Solution_file,"w");
  fprintf(SolutionOutput,"//  %s\n", Title);
  SolutionObject.Output(monitor);
  SolutionObject.Output(SolutionOutput);
  fclose(SolutionOutput);

  SurfaceOutput = fopen(Surface_file,"w");
  fprintf(SurfaceOutput,"//  %s\n",Title);

  // Output free surface
  fprintf(SurfaceOutput, "\n# Cnoidal Theory - %s",Method);
  fprintf(SurfaceOutput, "\n");
  fprintf(SurfaceOutput, "\n# Surface of wave - trough-crest-trough,");
  fprintf(SurfaceOutput, " note quadratic point spacing clustered around crest");
  fprintf(SurfaceOutput, "\n# Non-dimensionalised with respect to depth");
  fprintf(SurfaceOutput, "\n# X/d, eta/d, & check of surface pressure\n");

  for (i=-Surface_points/2; i<=Surface_points/2; ++i)
  {
   x = (double)i/Surface_points;
   x = LL*2*x*fabs(x);
   eta = SolutionObject.eta_h(x/hh);
   eta_d = eta*hh;
   uu = SolutionObject.u_h(x/hh, eta);
   vv = SolutionObject.v_h(x/hh, eta);
   fprintf(SurfaceOutput,"\n%8.4f\t%8.4f\t%8.0e", x, eta_d, 0.5*(uu*uu+vv*vv)+eta-RR_d/hh);
  }

  fclose(SurfaceOutput);

  // Output flowfield summary
  FlowfieldOutput = fopen(Flowfield_file,"w");
  fprintf(FlowfieldOutput,"//  %s", Title);

  fprintf(FlowfieldOutput,"\n\n# Cnoidal Theory - %s",Method);
  fprintf(FlowfieldOutput,"\n# Velocity profiles\n");
  fprintf(FlowfieldOutput,"\n# All quantities are dimensionless with respect to g and/or d\n");
  fprintf(FlowfieldOutput,"\n#*********************");
  fprintf(FlowfieldOutput,"\n# y        u       v");
  fprintf(FlowfieldOutput,"\n# -     -------------");
  fprintf(FlowfieldOutput,"\n# d        sqrt(gd)");
  fprintf(FlowfieldOutput,"\n#*********************");

  for(ii=0;ii<=Nprofiles;ii++)
  {
   x = LL*0.5*(double)ii/Nprofiles;
   eta = SolutionObject.eta_h(x/hh);
   eta_d = eta*hh;

   fprintf(FlowfieldOutput,"\n\n# x/d = %8.4f, Phase = %6.1f°\n", x, x/LL*360);

   for(i=0;i<=Points;i++)
   {
    y = i*eta/Points;
    uu = SolutionObject.u_h(x/hh, y)*sqrt(hh)+cc;
    vv = SolutionObject.v_h(x/hh, y)*sqrt(hh);
    fprintf(FlowfieldOutput,"\n%7.4f\t%7.4f\t%7.4f",y*hh,uu,vv);
   }
  }

  fclose(FlowfieldOutput);

  AccuracyOutput = fopen(Accuracy_file,"a");

  // Highest wave - eqn (32) of Fenton (1990)
  HHighest = Get_Highest(LL);

  double Velo[Points+1], sum1, sum2, ucm;

  x = 0.;
  eta = SolutionObject.eta_h(x/hh);
  eta_d = eta*hh;

  for(i=0 ; i <= Points; ++i)
  {
   y = (i)*eta/(Points);
   uu = SolutionObject.u_h(x/hh, y)*sqrt(hh)+cc;
   Velo[i] = uu;
   }

   for(i=1, sum1=0; i <= Points-1; i+=2) sum1 += Velo[i];
   for(i=2, sum2=0; i <= Points-2; i+=2) sum2 += Velo[i];

   ucm = (Velo[0]+4*sum1+2*sum2+Velo[Points])/3./Points;
   fprintf(AccuracyOutput,"\n%s %2d\t%7.4f\t%8.3f\t%7.3f\t%3.0f\t%8.5f", Method, n, HH, LL, SSU, HH/HHighest*100, ucm);
   fclose(AccuracyOutput);
 }
 else
 {
  assert(false);
 }

}

#endif

