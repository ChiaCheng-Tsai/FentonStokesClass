
#ifndef Stokes_H
#define Stokes_H
// Stokes theory calculations

#include <math.h>
#include <stdio.h>
#include <process.h>
#include <string.h>
#include <conio.h>
#include <stdlib.h>
#include <assert.h>
#define	ANSI

#include "FentonInout.h"

class StokesClass
{
 private:
  double C[8], D[8], E[8], e[8];
 public:
  StokesClass(FentonInoutClass<StokesClass> & );
  double F(double kd,int Inoutn,double InoutT,double InoutH,double InoutCurrent,int InoutCurrent_criterion);
  void CDE(double kd);
  void zBY(double kd,int Inoutn,double* Inoutz,double* InoutB,double* InoutY);
};

// read data (n, Case, T, L, MaxH, Current_criterion, Current) from Inout
// and set Method, z[], B[], Y[], H, kd, kH, and SU to Inout

StokesClass::StokesClass(FentonInoutClass<StokesClass> & Inout)
{

 double rhs, rhs1, rhs2, kd1, kd2, kFM, omega, delta;
 int i;

 if(Inout.n>7)Inout.n=7;

 sprintf(Inout.Method, "%d-order Stokes theory", Inout.n);
 printf( "\n# Solution by %d-order Stokes theory", Inout.n);

 Inout.z = new double [2*Inout.n+10+1]; // +1 by Tsai
 Inout.Y = new double [Inout.n+1]; // +1 by Tsai
 Inout.B = new double [Inout.n+1]; // +1 by Tsai
// see page 359 in Fenton (1988) for the definition of z and B
// see Eq. (3.5) in Fenton (1999) for the definition of B
// see page 363 in Fenton (1988) for the definition of Y

 Inout.H = Inout.MaxH;

 iff(Inout.Case,Wavelength)
	{
 	 Inout.kd = 2. * pi / Inout.L;
	 Inout.kH = Inout.kd * Inout.H ;
	 CDE(Inout.kd);
	}

// Solve dispersion relation using bisection

 iff(Inout.Case,Period)
	{
	 omega = 2*pi/Inout.T;
	// Fenton & McKee for initial estimate
	 kFM = omega*omega*pow(1/tanh(pow(omega,1.5)),(2./3.));
	 kd1 = 0.5*kFM;
	 kd2 = 2*kFM;
	 delta = 1.;
	 Inout.kd = kFM; // added by Tsai
	 CDE(kd1);
	 rhs1 = F(kd1,Inout.n,Inout.T,Inout.H,Inout.Current,Inout.Current_criterion);
	 CDE(kd2);
	 rhs2 = F(kd2,Inout.n,Inout.T,Inout.H,Inout.Current,Inout.Current_criterion);
	 if (rhs1*rhs2 >= 0. ) {printf("\nInterval wrong\n");assert(false);}
	 while ( fabs(delta) > 1.e-7 * Inout.kd )
		{
		 Inout.kd = 0.5 * (kd1+kd2);
		 delta = kd2-kd1;
		 CDE(Inout.kd);
		 rhs = F(Inout.kd,Inout.n,Inout.T,Inout.H,Inout.Current,Inout.Current_criterion);
		 if ( rhs*rhs1 < 0.)
			{kd2 = Inout.kd ;}
		 else
			{kd1 = Inout.kd ; rhs1 = rhs;}
		}
	 Inout.kH = Inout.kd * Inout.H ;
	}

//	Put in array ready for standard output routine

// Diagnostic

 Inout.z[1] = Inout.kd;
 Inout.z[2] = Inout.kH;

 Inout.SU = 0.5*Inout.kH/pow(Inout.kd,3);
 printf("\n# Stokes-Ursell no.: %7.3f",Inout.SU);
 if( Inout.SU > 0.5)
	{
	printf(" > 1/2. Results are unreliable\n");
	}

 e[1] = 0.5 * Inout.kH;
 for ( i=2 ; i<=Inout.n ; i++ ) e[i] = e[i-1] * e[1];

// Calculate coefficients

 CDE(Inout.kd);
 zBY(Inout.kd,Inout.n,Inout.z,Inout.B,Inout.Y);

 Inout.z[7] = C[0] + e[2]*C[2] + e[4] * C[4] + e[6] * C[6]; // ubar
 Inout.z[8] = - e[2]*D[2] - e[4]*D[4]- e[6]*D[6];
 Inout.z[9] = 0.5 * C[0]*C[0] + e[2]*E[2] + e[4] * E[4]+ e[6] * E[6];

 if(Inout.Current_criterion==1)
	{
	 Inout.z[5] = Inout.Current*sqrt(Inout.kd);
	 Inout.z[4] = Inout.z[7] + Inout.z[5];
	 Inout.z[6] = Inout.z[4] + Inout.z[8]/Inout.kd - Inout.z[7];
	}

 if(Inout.Current_criterion==2)
	{
	Inout.z[6] = Inout.Current*sqrt(Inout.kd);
	Inout.z[4] = Inout.z[6] - Inout.z[8]/Inout.kd + Inout.z[7];
	Inout.z[5] = Inout.z[4] - Inout.z[7];
	}

 iff(Inout.Case,Wavelength) Inout.z[3] = 2*pi/Inout.z[4];
 iff(Inout.Case,Period) Inout.z[3] = Inout.T * sqrt(Inout.kd);

}

// Evaluates dispersion relation - used in iterative solution for wavelength

double StokesClass::F(double kd,int Inoutn,double InoutT,double InoutH,double InoutCurrent,int InoutCurrent_criterion)
{

 double e2, rhs, kH, CC[4];
 kH = kd * InoutH ;
 e2 = 0.5*kH*0.5*kH;
 rhs = InoutCurrent * sqrt(kd) - 2.*pi/InoutT/sqrt(kd);
 CC[0] = C[0];
 CC[1] = CC[2] = CC[3] = 0.;

 // Eq (8) in Fention(1990)
 if ( Inoutn >= 3 )
	CC[1] = e2*C[2];
 if ( Inoutn >= 5 )
	CC[2] = e2*e2*C[4];
 if ( Inoutn >= 7 )
	CC[3] = e2*e2*e2*C[6];

 // Eq (11) in Fention(1990)
 if(InoutCurrent_criterion==2)
	{
	if ( Inoutn >= 3 )
		CC[1] += e2*D[2]/kd;
	if ( Inoutn >= 5 )
		CC[2] += e2*e2*D[4]/kd;
	if ( Inoutn >= 7 )
		CC[3] += e2*e2*e2*D[6]/kd;
	}

 CC[1] = CC[0]+CC[1];
 CC[2] = CC[1]+CC[2];
 CC[3] = CC[2]+CC[3];

 return rhs+CC[3];
}

// Calculate coefficient arrays C[], D[] and E[] re-derived by Tsai as in Table 3-1 of Fenton (1990)
// give kd and solve C[], D[], and E[]

void StokesClass::CDE(double kd)
{
 C[0]=pow(tanh(kd),0.5);
 C[2]=pow(tanh(kd),0.5)*(((8 + cosh(4*kd))*pow((1/sinh(kd)),4))/16.);
 C[4]=pow(tanh(kd),0.5)*(((-168 - 622*cosh(2*kd) - 736*cosh(4*kd) - 111*cosh(6*kd) + 16*cosh(8*kd) + cosh(10*kd))*pow((1/sinh(kd)),10))/4096.);
 C[6]=pow(tanh(kd),0.5)*(((7698678 + 14301028*cosh(2*kd) + 11521240*cosh(4*kd) + 8040318*cosh(6*kd) + 4762392*cosh(8*kd) +
        2175570*cosh(10*kd) + 500934*cosh(12*kd) - 2071*cosh(14*kd) - 9310*cosh(16*kd) + 3*cosh(18*kd) + 18*cosh(20*kd))
       *pow((1/sinh(kd)),16))/(3.145728e6*(8 + 11*cosh(2*kd) + 6*cosh(4*kd))));

 D[2]=(1/sqrt(tanh(kd)))*(-1./2.);
 D[4]=(1/sqrt(tanh(kd)))*((8 + 5*cosh(2*kd) + 4*cosh(4*kd) + cosh(6*kd))*pow(1/sinh(kd),6))/128.;
 D[6]=(1/sqrt(tanh(kd)))*(-(pow(1/sinh(kd),13)*(2375*sinh(kd) + 12611*sinh(3*kd) + 19507*sinh(5*kd) + 23003*sinh(7*kd) + 15467*sinh(9*kd) + 5747*sinh(11*kd) + 1507*sinh(13*kd) - 79*sinh(15*kd) - 70*sinh(17*kd)))/
    (131072.*(2 + 3*cosh(2*kd))*(1 + 4*cosh(2*kd))));

 E[2]=((6 + 2*cosh(2*kd) + cosh(4*kd))*pow((1/sinh(kd)),3)*(1/cosh(kd)))/16.;
 E[4]=((-145 - 302*cosh(2*kd) - 296*cosh(4*kd) - 71*cosh(6*kd) + 3*cosh(8*kd) + cosh(10*kd))*pow((1/sinh(kd)),9)*(1/cosh(kd)))/
    2048.;
 E[6]=((1943286 + 3620356*cosh(2*kd) + 2912560*cosh(4*kd) + 2002512*cosh(6*kd) + 1139076*cosh(8*kd) + 496104*cosh(10*kd) +
       127743*cosh(12*kd) + 7640*cosh(14*kd) - 1954*cosh(16*kd) - 132*cosh(18*kd) + 9*cosh(20*kd))*pow((1/sinh(kd)),15))/
   (393216.*(27*cosh(kd) + 17*cosh(3*kd) + 6*cosh(5*kd)));

 return;
}

// Calculate coefficient arrays A[][] and B[][] re-derived by Tsai as in Table 3-1 of Fenton (1990)
// see page 359 in Fenton (1988) for the definition of z and B
// see Eq. (3.5) in Fenton (1999) for the definition of B
// see page 363 in Fenton (1988) for the definition of Y
// give kd & n and solve  z[], B[] and Y[]

void StokesClass::zBY(double kd,int Inoutn,double* Inoutz,double* InoutB,double* InoutY) // the same as Fenton's AB()
{

 int i, j;
 double A[8][8], BB[8][8];

 for(i=0; i<=7 ; i++) for(j=0; j<=7 ; j++) {A[i][j] = BB[i][j] = 0;}

 A[1][1]=(1/sinh(kd));
 A[2][2]=(3*pow((1/sinh(kd)),4))/8.;
 A[3][1]=-((23 - 7*cosh(2*kd) + 10*cosh(4*kd) + cosh(6*kd))*pow((1/sinh(kd)),7))/64.;
 A[3][3]=-((-11 + 2*cosh(2*kd))*pow((1/sinh(kd)),7))/64.;
 A[4][2]=-((281 + 111*cosh(2*kd) + 252*cosh(4*kd) + 7*cosh(6*kd) - 3*cosh(8*kd))*pow((1/sinh(kd)),10))/1536.;
 A[4][4]=((382 + 597*cosh(2*kd) - 174*cosh(4*kd) + 5*cosh(6*kd))*pow((1/sinh(kd)),10))/(3072.*(2 + 3*cosh(2*kd)));
 A[5][1]=((42975 + 78648*cosh(2*kd) + 63618*cosh(4*kd) + 39736*cosh(6*kd) + 19358*cosh(8*kd) + 5442*cosh(10*kd) + 1358*cosh(12*kd) + 2*cosh(14*kd) - 37*cosh(16*kd))*pow((1/sinh(kd)),13))/
   (16384.*(2 + 3*cosh(2*kd))*(1 + 4*cosh(2*kd)));
 A[5][3]=((-8280 - 17334*cosh(2*kd) - 8802*cosh(4*kd) - 4979*cosh(6*kd) + 408*cosh(8*kd) + 105*cosh(10*kd) + 2*cosh(12*kd))*pow((1/sinh(kd)),13))/(32768.*(2 + 3*cosh(2*kd)));
 A[5][5]=((7664 + 6890*cosh(2*kd) + 4496*cosh(4*kd) - 3119*cosh(6*kd) + 272*cosh(8*kd) - 3*cosh(10*kd))*pow((1/sinh(kd)),13))/(32768.*(8 + 11*cosh(2*kd) + 6*cosh(4*kd)));
 A[6][6]=((1244125 + 1928725*cosh(2*kd) + 1141396*cosh(4*kd) + 186684*cosh(6*kd) - 107900*cosh(8*kd) - 176150*cosh(10*kd) +
          37100*cosh(12*kd) - 1487*cosh(14*kd) + 7*cosh(16*kd))*pow((1/sinh(kd)),16))/
          (655360.*(2 + 3*cosh(2*kd))*(1 + 4*cosh(2*kd))*(9 + 16*cosh(2*kd) + 10*cosh(4*kd)));
 A[6][4]=-((11557605 + 21344025*cosh(2*kd) + 14643352*cosh(4*kd) + 8170653*cosh(6*kd) + 2538500*cosh(8*kd) +
         385795*cosh(10*kd) - 300960*cosh(12*kd) - 20314*cosh(14*kd) + 1359*cosh(16*kd) - 15*cosh(18*kd))*
       pow((1/sinh(kd)),16))/(3.93216e6*pow(2 + 3*cosh(2*kd),2)*(1 + 4*cosh(2*kd)));
 A[6][2]=((2497110 + 4577013*cosh(2*kd) + 3657026*cosh(4*kd) + 2329161*cosh(6*kd) + 1190596*cosh(8*kd) + 398411*cosh(10*kd) +
        81414*cosh(12*kd) - 3161*cosh(14*kd) - 1794*cosh(16*kd) + 24*cosh(18*kd))*pow((1/sinh(kd)),16))/
    (786432.*(8 + 11*cosh(2*kd) + 6*cosh(4*kd)));
 A[7][7]=((735085189 + 1227947077*cosh(2*kd) + 910795886*cosh(4*kd) + 382064899*cosh(6*kd) + 117632588*cosh(8*kd) -
        87264287*cosh(10*kd) - 49722049*cosh(12*kd) - 36440942*cosh(14*kd) + 16185007*cosh(16*kd) -
        1422031*cosh(18*kd) + 28723*cosh(20*kd) - 60*cosh(22*kd))*pow((1/sinh(kd)),19))/
    (9.437184e7*(2086 + 3823*cosh(2*kd) + 2984*cosh(4*kd) + 1906*cosh(6*kd) + 992*cosh(8*kd) + 369*cosh(10*kd) +
        90*cosh(12*kd)));
 A[7][5]=-((15290099368 + 28211179102*cosh(2*kd) + 22243034846*cosh(4*kd) + 14611400830*cosh(6*kd) + 7841536046*cosh(8*kd) +
         3118965289LL*cosh(10*kd) + 764085674*cosh(12*kd) - 63117932*cosh(14*kd) - 117813140*cosh(16*kd) -
         48456514*cosh(18*kd) + 2544568*cosh(20*kd) + 556767*cosh(22*kd) - 14850*cosh(24*kd) - 54*cosh(26*kd))*
       pow((1/sinh(kd)),19))/
    (7.5497472e7*pow(8 + 11*cosh(2*kd) + 6*cosh(4*kd),2)*(9 + 16*cosh(2*kd) + 10*cosh(4*kd)));
 A[7][3]=((408506090 + 762792180*cosh(2*kd) + 608494934*cosh(4*kd) + 412075566*cosh(6*kd) + 224487020*cosh(8*kd) +
        94046200*cosh(10*kd) + 25072830*cosh(12*kd) + 2588037*cosh(14*kd) - 995062*cosh(16*kd) - 158585*cosh(18*kd) +
        10060*cosh(20*kd) + 730*cosh(22*kd))*pow((1/sinh(kd)),19))/
    (4.194304e7*pow(2 + 3*cosh(2*kd),2)*(1 + 4*cosh(2*kd)));
 A[7][1]=-((1589018985462 + 3054193939488*cosh(2*kd) + 2709219165939*cosh(4*kd) + 2213997416553*cosh(6*kd) +
         1661336214690*cosh(8*kd) + 1138970490972*cosh(10*kd) + 708350113572*cosh(12*kd) + 395666005708*cosh(14*kd) +
         195724460131*cosh(16*kd) + 84036007536*cosh(18*kd) + 30388586828*cosh(20*kd) + 8831541434*cosh(22*kd) +
         1895617494*cosh(24*kd) + 249006441*cosh(26*kd) + 4179693*cosh(28*kd) - 3919437*cosh(30*kd) +
         172287*cosh(32*kd) + 115209*cosh(34*kd))*pow((1/sinh(kd)),19))/
    (7.5497472e7*pow(8 + 11*cosh(2*kd) + 6*cosh(4*kd),2)*
      (49 + 84*cosh(2*kd) + 63*cosh(4*kd) + 34*cosh(6*kd) + 15*cosh(8*kd)));

 for ( i = 1 ; i <= Inoutn ; i++ )
 {
  for ( Inoutz[Inoutn+10+i] = 0., j = ((i+1) % 2) + 1 ; j <= Inoutn ; j += 2 )
	 	Inoutz[Inoutn+10+i] += A[j][i] * e[j];

  Inoutz[Inoutn+10+i] *= C[0]	* cosh(i*kd);
  InoutB[i] = Inoutz[Inoutn+10+i];
 }

 BB[1][1]= 1.;
 BB[2][2]=((2 + cosh(2*kd))*(cosh(kd)/sinh(kd))*pow((1/sinh(kd)),2))/4.;
 BB[3][1]= (-3*(14 + 15*cosh(2*kd) + 6*cosh(4*kd) + cosh(6*kd))*pow((1/sinh(kd)),6))/256.;
 BB[3][3]= (3*(14 + 15*cosh(2*kd) + 6*cosh(4*kd) + cosh(6*kd))*pow((1/sinh(kd)),6))/256.;
 BB[4][2]=((-264*cosh(kd) - 262*cosh(3*kd) - 114*cosh(5*kd) - 9*cosh(7*kd) + cosh(9*kd))*pow((1/sinh(kd)),9))/768.;
 BB[4][4]=((753*cosh(kd) + 505*cosh(3*kd) + 234*cosh(5*kd) + 99*cosh(7*kd) + 26*cosh(9*kd) + 3*cosh(11*kd))*pow((1/sinh(kd)),9))/(1536.*(2 + 3*cosh(2*kd)));
 BB[5][1]=((2081002 + 3867653*cosh(2*kd) + 3021600*cosh(4*kd) + 1910415*cosh(6*kd) + 884040*cosh(8*kd) + 259709*cosh(10*kd) + 33824*cosh(12*kd) - 4177*cosh(14*kd) - 1266*cosh(16*kd))*
      pow((1/sinh(kd)),12))/(786432.*(2 + 3*cosh(2*kd))*(1 + 4*cosh(2*kd)));
 BB[5][3]=(-9*(9942 + 18596*cosh(2*kd) + 14465*cosh(4*kd) + 6930*cosh(6*kd) + 1866*cosh(8*kd) + 74*cosh(10*kd) - 33*cosh(12*kd))*pow((cosh(kd)/sinh(kd)),2)*pow((1/sinh(kd)),10))/
    (131072.*(2 + 3*cosh(2*kd)));
 BB[5][5]=(5*(148034 + 258553*cosh(2*kd) + 184128*cosh(4*kd) + 102939*cosh(6*kd) + 51816*cosh(8*kd) + 22849*cosh(10*kd) + 7552*cosh(12*kd) + 1579*cosh(14*kd) + 150*cosh(16*kd))*
      pow((1/sinh(kd)),12))/(1.572864e6*(2 + 3*cosh(2*kd))*(1 + 4*cosh(2*kd)));
 BB[6][6]=(3*(9890490*cosh(kd) + 8553720*cosh(3*kd) + 6427500*cosh(5*kd) + 4239300*cosh(7*kd) + 2504046*cosh(9*kd) +
          1346769*cosh(11*kd) + 663135*cosh(13*kd) + 275675*cosh(15*kd) + 92475*cosh(17*kd) + 22977*cosh(19*kd) +
          3643*cosh(21*kd) + 270*cosh(23*kd))*pow((1/sinh(kd)),15))/
          (2.62144e6*(8 + 11*cosh(2*kd) + 6*cosh(4*kd))*(9 + 16*cosh(2*kd) + 10*cosh(4*kd)));
 BB[6][4]=-((11834620*cosh(kd) + 10298830*cosh(3*kd) + 7728820*cosh(5*kd) + 4906890*cosh(7*kd) + 2578326*cosh(9*kd) +
         1088709*cosh(11*kd) + 355850*cosh(13*kd) + 80955*cosh(15*kd) + 8105*cosh(17*kd) - 888*cosh(19*kd) -
         217*cosh(21*kd))*pow((1/sinh(kd)),15))/(655360.*pow(2 + 3*cosh(2*kd),2)*(1 + 4*cosh(2*kd)));
 BB[6][2]=((6248194*cosh(kd) + 5393974*cosh(3*kd) + 3972042*cosh(5*kd) + 2419792*cosh(7*kd) + 1162944*cosh(9*kd) +
        382334*cosh(11*kd) + 57804*cosh(13*kd) - 1669*cosh(15*kd) - 989*cosh(17*kd) - 26*cosh(19*kd))*
      pow((1/sinh(kd)),15))/(524288.*(8 + 11*cosh(2*kd) + 6*cosh(4*kd)));
 BB[7][7]=(7*(50371994790 + 95863961142*cosh(2*kd) + 82789481811*cosh(4*kd) + 64953901428*cosh(6*kd) + 46626715989*cosh(8*kd) +
        30791907220*cosh(10*kd) + 18932238035*cosh(12*kd) + 10866296014*cosh(14*kd) + 5795909850*cosh(16*kd) +
        2790769026LL*cosh(18*kd) + 1165429125*cosh(20*kd) + 411851612*cosh(22*kd) + 117186859*cosh(24*kd) +
        24705788*cosh(26*kd) + 3355221*cosh(28*kd) + 216090*cosh(30*kd))*pow((1/sinh(kd)),18))/
    (2.415919104e10*(2 + 2*cosh(2*kd) + 3*cosh(4*kd))*
      (190 + 330*cosh(2*kd) + 222*cosh(4*kd) + 103*cosh(6*kd) + 30*cosh(8*kd)));
 BB[7][5]=(-5*(741113326128 + 1412241577590*cosh(2*kd) + 1220673956160*cosh(4*kd) + 955292747595*cosh(6*kd) +
        674672429454*cosh(8*kd) + 428097257894*cosh(10*kd) + 242640146062*cosh(12*kd) + 121906505321*cosh(14*kd) +
        53674003296*cosh(16*kd) + 20295321429*cosh(18*kd) + 6334240062*cosh(20*kd) + 1497759376*cosh(22*kd) +
        217827554*cosh(24*kd) + 4033999*cosh(26*kd) - 4564236*cosh(28*kd) - 567684*cosh(30*kd))*pow((1/sinh(kd)),18))/
    (4.831838208e9*pow(2 + 3*cosh(2*kd),2)*pow(1 + 4*cosh(2*kd),2)*(9 + 16*cosh(2*kd) + 10*cosh(4*kd)));
 BB[7][3]=(9*(7865383120 + 14797192490*cosh(2*kd) + 12289301020*cosh(4*kd) + 8940194302*cosh(6*kd) + 5605859458*cosh(8*kd) +
        2942913215LL*cosh(10*kd) + 1230777302*cosh(12*kd) + 376502213*cosh(14*kd) + 70568224*cosh(16*kd) +
        3501761*cosh(18*kd) - 1140210*cosh(20*kd) - 100461*cosh(22*kd) + 7566*cosh(24*kd))*pow((1/sinh(kd)),18))/
    (2.68435456e9*pow(2 + 3*cosh(2*kd),2)*(1 + 4*cosh(2*kd)));
 BB[7][1]=-((204682036568508 + 393228562989111*cosh(2*kd) + 348317247555996*cosh(4*kd) + 283950853772844*cosh(6*kd) +
         212300376997076*cosh(8*kd) + 144817641948188*cosh(10*kd) + 89451964321344*cosh(12*kd) +
         49502099109029*cosh(14*kd) + 24171043702496*cosh(16*kd) + 10181723004279*cosh(18*kd) +
         3572345158408*cosh(20*kd) + 983079076060*cosh(22*kd) + 187360159500*cosh(24*kd) + 16258146268*cosh(26*kd) -
         1871309540*cosh(28*kd) - 405408279*cosh(30*kd) + 67001220*cosh(32*kd) + 14007492*cosh(34*kd))*
       pow((1/sinh(kd)),18))/
    (9.663676416e9*pow(8 + 11*cosh(2*kd) + 6*cosh(4*kd),2)*
      (49 + 84*cosh(2*kd) + 63*cosh(4*kd) + 34*cosh(6*kd) + 15*cosh(8*kd)));

 for ( i = 1 ; i <= Inoutn ; i++ )
  for ( InoutY[i] = 0., j = ((i+1) % 2) + 1 ; j <= Inoutn ; j += 2 )
   InoutY[i] += BB[j][i] * e[j];

 return;
}

#endif

