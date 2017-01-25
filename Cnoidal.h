
#ifndef Cnoidal_H
#define Cnoidal_H
// Cnoidal theory calculations

#include <math.h>
#include <stdio.h>
#include <process.h>
#include <string.h>
#include <conio.h>
#include <stdlib.h>
#include <assert.h>
#define	ANSI

#include "FentonInout.h"
#include "FentonElliptic.h"

char 		*Level[7]={(char *)"",(char *)"Cnoidal 1st Order theory",(char *)"Cnoidal 2nd Order theory",
                       (char *)"Cnoidal 3rd Order theory",(char *)"Cnoidal 4th Order theory",
                       (char *)"Cnoidal 5th Order theory",(char *)"Cnoidal 5th order theory with Aitken convergence enhancement"};

class CnoidalClass
{
 private:
  // FentonInout input variables
  double Current,H,T,L;
  int Current_criterion,n;
  char Title[100],Case[20],Method[100];

  // Class internal variables
  double c,ce,cs,m,e,epsilon,k,h,alpha,delta,Ubar_d,Q_d,R_d,SU;

  double ee[7], mm[7];
  void Set_mmee(double Inoutm,double Inoute);
  double Ubar_h(int Inoutn,double Inoutm,double Inoute,double Inoutepsilon);
  double Q_h(int Inoutn,double Inoutm,double Inoute,double Inoutepsilon);
  double lambda_d(int Inoutn,double Inoutm,double Inoute,double InoutH);
  double lambda_series(int Inoutn,double Inoutm,double Inoute,double InoutH);
  double hoverd(int Inoutn,double Inoutm,double Inoute,double InoutH);
  double Alpha(int Inoutn,double Inoutm,double Inoute,double Inoutepsilon);
  double R_h(int Inoutn,double Inoutm,double Inoute,double Inoutepsilon);
  double Aitken(double *S, int j);
 public:
  CnoidalClass(){};
  CnoidalClass(FentonInout<CnoidalClass> & );
  double Get_h(){return h;};
  double Get_R_d(){return R_d;};
  double Get_L(){return L;};
  double Get_H(){return H;};
  double Get_SU(){return SU;};
  double Get_c(){return c;};
  int Get_n(){return n;};
  double eta_h(double x);
  double u_h(double x, double Y);
  double v_h(double x, double Y);
  void Output(FILE *Solution);

  // dummy member function
  double Get_kd(){assert(false);return 0;};
  double Get_kh(){assert(false);return 0;};
  double Get_kH(){assert(false);return 0;};
  double Get_z(int i){assert(false);return i*0;};
  double Surface(double x){assert(false);return x*0;};
  void Point(double X,double y,
             double &phi,double &u,double &v,double &dphidt,double &dudt,double &dvdt,double &ut,double &vt,
             double &ux,double &vx,double &uy,double &vy,double &Pressure,double &Bernoulli_check)
             {assert(false);return;};
};

// read data (n, Case, T, L, H, Current_criterion, Current) from Inout
// and set alpha, delta, m, KK, e, epsilon, Method, and others to Inout

CnoidalClass::CnoidalClass(FentonInout<CnoidalClass> & Inout)
{
 int iter;
 double KK,q1,UU1;

 Current=Inout.Current;
 Current_criterion=Inout.Current_criterion;
 H=Inout.H;
 T=Inout.T;
 L=Inout.L;
 n=Inout.n;

 strcpy(Case,Inout.Case);
 strcpy(Title,Inout.Title);

 if (Current_criterion==1) ce = Current;
 if (Current_criterion==2) cs = Current;
//Use direct iteration to solve for m
 printf("\n\nIteration for m (1-m shown)");

 iff(Case, Period)
   m = 1.-16.*exp(-T*sqrt(0.75*H));
 iff(Case, Wavelength)
   m = 1.-16.*exp(-L*sqrt(0.75*H));

 for(iter=1;iter<=5;iter++)
	{
	 e=E(m)/K(m);
  	 epsilon=H/hoverd(n,m,e,H);

     Set_mmee(m,e);

     KK=0;
	 iff(Case, Period)
		{
		 if (Current_criterion == 1)
			KK=T*(ce+Ubar_h(n,m,e,epsilon)*sqrt(hoverd(n,m,e,H)))*sqrt(3*H/m)/4/lambda_series(n,m,e,H);
		 if (Current_criterion == 2)
			KK=T*(cs+Q_h(n,m,e,epsilon)*pow(hoverd(n,m,e,H),1.5))*sqrt(3*H/m)/4/lambda_series(n,m,e,H);
		}
	 iff(Case, Wavelength)
		{
		 KK=L*sqrt(3*H/m)/4/lambda_series(n,m,e,H);
		}
	 q1=exp(-Fenton_pi*KK/K(1.-m));
	 m=pow((1-2*q1+2*pow(q1,4))/(1+2*q1+2*pow(q1,4)),4);
	 printf("\n%14.6e", 1.-m);
	}

 e = E(m)/K(m);

 Set_mmee(m,e);

 k = sqrt(m);
 h = hoverd(n,m,e,H);
 epsilon = H/h;
 alpha = Alpha(n,m,e,epsilon);
 delta = 4./3*alpha*alpha;
 Ubar_d = Ubar_h(n,m,e,epsilon)*sqrt(h);
 Q_d = Q_h(n,m,e,epsilon)*pow(h,1.5);

 if (Current_criterion==1)
   {
	c = ce+Ubar_d;
    cs = c-Q_d;
   }
 if (Current_criterion==2)
   {
    c = cs+Q_d;
    ce = c-Ubar_d;
   }

 iff(Case,Wavelength) T = L/c;
 iff(Case,Period) L = lambda_d(n,m,e,H);
 R_d = R_h(n,m,e,epsilon)*h;
 UU1 = H*L*L;
 SU = UU1/8/Fenton_pi/Fenton_pi;
 if(SU < 0.5) printf("\nThe program is being applied to a wave with Stokes-Ursell number less than 1/2. Results are unreliable");

 strcpy(Inout.Method,Level[n]);
 strcpy(Method,Inout.Method);

}

void CnoidalClass::Set_mmee(double Inoutm,double Inoute)
{
 int i;

 ee[1] = Inoute;
 mm[1] = Inoutm;

 for(i=2; i<=5 ; ++i)
 {
	ee[i] = ee[i-1]*Inoute;
	mm[i] = mm[i-1]*Inoutm;
 }
 return;
}

//######################################################################################
// Ubar_h ( H/h, m)
// U_bar_Height is U/sqrt(g.h) as a function of H/h, m, and e(m), Fenton (1999, eqn A.6)
//######################################################################################

double CnoidalClass::Ubar_h(int Inoutn,double Inoutm,double Inoute,double Inoutepsilon)
{
 double U_bar_Height[7];

 Set_mmee(Inoutm,Inoute);

 U_bar_Height[1] = 1+Inoutepsilon/Inoutm*(1./2-Inoute);
 U_bar_Height[2] = U_bar_Height[1]+pow(Inoutepsilon,2)/mm[2]*(-13./120-1./60*Inoutm-1./40*mm[2]+(1./3+1./12*Inoutm)*Inoute);
 U_bar_Height[3] = U_bar_Height[2]+pow(Inoutepsilon,3)/mm[3]*(-361./2100+1899./5600*Inoutm-2689./16800*mm[2]+13./280*mm[3]+(7./75-103./300*Inoutm+131./600*mm[2])*Inoute);
 U_bar_Height[4] = U_bar_Height[3]+pow(Inoutepsilon,4)/mm[4]*(2349./112000+29053./168000*Inoutm-1181./2100*mm[2]+11161./28000*mm[3]-273./3200*mm[4]+(29369./28000*mm[2]-15867./28000*mm[3]-5729./8400*Inoutm+1583./4200)*Inoute);
 U_bar_Height[5] = U_bar_Height[4]+pow(Inoutepsilon,5)/mm[5]*(1786123./16170000-32376301./103488000*Inoutm-87413873./776160000*mm[2]+474001783./517440000*mm[3]-71678951./97020000*mm[4]+97103./616000*mm[5]+(-61854593./35280000*mm[3]+35498549./35280000*mm[4]+444959./1260000*mm[2]+1196747./980000*Inoutm-691177./735000)*Inoute);

 U_bar_Height[6] = Aitken(U_bar_Height,5);

 return U_bar_Height[Inoutn];
}

//##############################################################################
// Q_h ( H/h, m)
// Q_Height is Q/sqrt(g.h^3) as a function of H/h and m, Fenton (1999, eqn A.4)
//##############################################################################

double CnoidalClass::Q_h(int Inoutn,double Inoutm,double Inoute,double Inoutepsilon)
{
 double Q_Height[7];

 Set_mmee(Inoutm,Inoute);

 Q_Height[1] = 1+Inoutepsilon/Inoutm*(-1./2+Inoutm);
 Q_Height[2] = Q_Height[1]+pow(Inoutepsilon,2)/mm[2]*(9./40-7./20*Inoutm-1./40*mm[2]);
 Q_Height[3] = Q_Height[2]+pow(Inoutepsilon,3)/mm[3]*(-11./140+69./1120*Inoutm+11./224*mm[2]+3./140*mm[3]);
 Q_Height[4] = Q_Height[3]+pow(Inoutepsilon,4)/mm[4]*(-16109./42000*mm[3]+59321./84000*mm[2]-123787./168000*Inoutm-871./22400*mm[4]+133687./336000);
 Q_Height[5] = Q_Height[4]+pow(Inoutepsilon,5)/mm[5]*(89101./1232000*mm[5]+7482007./16170000*mm[4]-4473257./5390000-347331631./517440000*mm[3]-21859819./36960000*mm[2]+163246841./103488000*Inoutm);

 Q_Height[6] = Aitken(Q_Height,5);

 return Q_Height[Inoutn];
}

//#################################################################################
// lambda_d ( H/d, m)
// Wavelength is lamda/d as a function of H/d, m, and e(m), Fenton (1999, eqn A.7)
//#################################################################################

double CnoidalClass::lambda_d(int Inoutn,double Inoutm,double Inoute,double InoutH)
{
 return K(Inoutm)*4/sqrt(3*InoutH/Inoutm)*lambda_series(Inoutn, Inoutm, Inoute,InoutH);
}

double CnoidalClass::lambda_series(int Inoutn,double Inoutm,double Inoute,double InoutH)
{
 double Wavelength[7];

 Set_mmee(Inoutm,Inoute);

 Wavelength[1] = 1;
 Wavelength[2] = Wavelength[1]+InoutH/Inoutm*(-3./2*Inoute+5./4-5./8*Inoutm);
 Wavelength[3] = Wavelength[2]+pow(InoutH,2)/mm[2]*(-15./32+15./32*Inoutm-21./128*mm[2]+(-1./16*Inoutm+1./8)*Inoute+3./8*ee[2]);
 Wavelength[4] = Wavelength[3]+pow(InoutH,3)/mm[3]*(341227./336000-341227./224000*Inoutm+984359./1344000*mm[2]-20127./179200*mm[3]+(-1471./1600-409./6400*mm[2]+1471./1600*Inoutm)*Inoute+(-7./64*Inoutm+7./32)*ee[2]+1./16*ee[3]);
 Wavelength[5] = Wavelength[4]+pow(InoutH,4)/mm[4]*(-105363683./37632000+105363683./18816000*Inoutm-306621467./75264000*mm[2]+95894101./75264000*mm[3]-1575087./28672000*mm[4]+(-2462811./448000*Inoutm+820937./224000-1086367./1792000*mm[3]+2728241./896000*mm[2])*Inoute+(-9501./6400-2679./25600*mm[2]+9501./6400*Inoutm)*ee[2]+(13./64-13./128*Inoutm)*ee[3]+3./128*ee[4]);

 Wavelength[6] =  Aitken(Wavelength, 5);

 return Wavelength[Inoutn];
}

//################################################
// h/d ( H, m)
// Fenton (1999, eqn A.8)
//################################################

double CnoidalClass::hoverd(int Inoutn,double Inoutm,double Inoute,double InoutH)
{
 double hd[7];

 Set_mmee(Inoutm,Inoute);

 hd[1] = 1+InoutH/Inoutm*(1-Inoute-Inoutm);
 hd[2] = hd[1]+pow(InoutH/Inoutm,2)*(-1./2+1./2*Inoutm+(1./2-1./4*Inoutm)*Inoute);
 hd[3] = hd[2]+pow(InoutH/Inoutm,3)*(133./200-399./400*Inoutm+133./400*mm[2]+(233./200*Inoutm-1./25*mm[2]-233./200)*Inoute+(1./2-1./4*Inoutm)*ee[2]);
 hd[4] = hd[3]+pow(InoutH/Inoutm,4)*(-122./75+244./75*Inoutm-1227./500*mm[2]+1241./1500*mm[3]+(481./150+6529./3000*mm[2]-573./2000*mm[3]-481./100*Inoutm)*Inoute+(52./25*Inoutm-57./400*mm[2]-52./25)*ee[2]+(1./2-1./4*Inoutm)*ee[3]);
 hd[5] = hd[4]+pow(InoutH/Inoutm,5)*(57231077./11760000-57231077./4704000*Inoutm+69379843./5880000*mm[2]-130123673./23520000*mm[3]+123967./120000*mm[4]+(126350477./5880000*Inoutm+2579201./490000*mm[3]-302159./1470000*mm[4]-26893043./1680000*mm[2]-126350477./11760000)*Inoute+(24361./4000*mm[2]-1779./2000*mm[3]-10347./800*Inoutm+3449./400)*ee[2]+(-123./400*mm[2]+649./200*Inoutm-649./200)*ee[3]+(1./2-1./4*Inoutm)*ee[4]);

 hd[6] = Aitken(hd, 5);

 return hd[Inoutn];
}

//####################################################################
// alpha ( H/h)
// Alpha is alpha as a function of H/h and m, Fenton (1999, eqn A.2)
//####################################################################

double CnoidalClass::Alpha(int Inoutn,double Inoutm,double Inoute,double Inoutepsilon)
{
 int i;
 double alpha[7];

 Set_mmee(Inoutm,Inoute);

 alpha[1] = 1;
 alpha[2] = alpha[1]+Inoutepsilon/Inoutm*(1./4-7./8*Inoutm);
 alpha[3] = alpha[2]+pow(Inoutepsilon,2)/mm[2]*(1./32-11./32*Inoutm+111./128*mm[2]);
 alpha[4] = alpha[3]+pow(Inoutepsilon,3)/mm[3]*(184711./1344000*mm[2]+114567./224000*Inoutm-126817./336000-149273./179200*mm[3]);
 alpha[5] = alpha[4]+pow(Inoutepsilon,4)/mm[4]*(13618217./25088000*mm[3]+22012297./28672000*mm[4]-34858533./25088000*mm[2]+509843./2508800+2777099./6272000*Inoutm);

 alpha[6] = Aitken(alpha,5);

 for(i=1; i <=6 ; ++i)
	alpha[i] =  sqrt(3*Inoutepsilon/4/Inoutm) * alpha[i];

 return alpha[Inoutn];
}

//########################################################################
// R_h ( H/h)
// R_Height is R/(g.h) as a function of H/h and m, Fenton (1999, eqn A.5)
//########################################################################

double CnoidalClass::R_h(int Inoutn,double Inoutm,double Inoute,double Inoutepsilon)
{
 double R_Height[7];

 Set_mmee(Inoutm,Inoute);

 R_Height[1] = 3./2+Inoutepsilon/Inoutm*(-1./2+Inoutm);
 R_Height[2] = R_Height[1]+pow(Inoutepsilon,2)/mm[2]*(-7./20*Inoutm+7./20-1./40*mm[2]);
 R_Height[3] = R_Height[2]+pow(Inoutepsilon,3)/mm[3]*(25./224*Inoutm-107./560+13./1120*mm[2]+13./280*mm[3]);
 R_Height[4] = R_Height[3]+pow(Inoutepsilon,4)/mm[4]*(-30823./42000*Inoutm+55331./84000*mm[2]+1214./2625-26833./84000*mm[3]-17./200*mm[4]);
 R_Height[5] = R_Height[4]+pow(Inoutepsilon,5)/mm[5]*(-270759631./258720000+24097./154000*mm[5]+21098053./64680000*mm[4]-202951241./517440000*mm[3]-864417./880000*mm[2]+198968527./103488000*Inoutm);

 R_Height[6] = Aitken(R_Height,5);

 return R_Height[Inoutn];
}

//########################
// Convergence enhancement
//########################

// Function to calculate Aitken transform

double CnoidalClass::Aitken(double *S, int j)
{
 double den, R;
 den  = S[j]+S[j-2]-2*S[j-1];
 if (fabs(den) < 1e-6) R = S[j];
 else R = S[j]-pow(S[j]-S[j-1],2)/den;
 return(R);
}

//####################
// eta_h (x/h)
//####################

double CnoidalClass::eta_h(double x)
{

//print(`"Eta" is eta/h as a function of H/h, m, and cn^2, Fenton (1999, eqn A.1)`);

 double Inoutalpha=alpha;
 double Inoutk=k;
 int Inoutn=n;
 double Inoutm=m;
 double Inoute=e;
 double Inoutepsilon=epsilon;

 int		i;
 double	Eta[7], C[11];

 C[1] = cn(Inoutalpha*x, Inoutk);

 for(i=2;i<=10;i++)
  C[i] = C[i-1]*C[1];

 Set_mmee(Inoutm,Inoute);

 Eta[1] = 1+C[2]*Inoutepsilon;
 Eta[2] = Eta[1]+pow(Inoutepsilon,2)/mm[2]*(-3./4*mm[2]*C[2]+3./4*mm[2]*C[4]);
 Eta[3] = Eta[2]+pow(Inoutepsilon,3)/mm[3]*((111./80*mm[3]-61./80*mm[2])*C[2]+(-53./20*mm[3]+61./80*mm[2])*C[4]+101./80*mm[3]*C[6]);
 Eta[4] = Eta[3]+pow(Inoutepsilon,4)/mm[4]*((59737./24000*mm[3]-4883./1600*mm[4]-302./375*mm[2])*C[2]+(-20791./4800*mm[3]+302./375*mm[2]+35551./4800*mm[4])*C[4]+(22109./12000*mm[3]-156611./24000*mm[4])*C[6]+17367./8000*mm[4]*C[8]);
 Eta[5] = Eta[4]+pow(Inoutepsilon,5)/mm[5]*((209511./32000*mm[5]-2209587./313600*mm[4]+3014947./1568000*mm[3]+684317./1568000*mm[2])*C[2]+(-2218593./112000*mm[5]-684317./1568000*mm[2]+3910057./224000*mm[4]-114211./24500*mm[3])*C[4]+(-490143./32000*mm[4]+40547./1600*mm[5]+4294557./1568000*mm[3])*C[6]+(-12800441./784000*mm[5]+7694543./1568000*mm[4])*C[8]+1331817./313600*mm[5]*C[10]);

 Eta[6] = Aitken(Eta,5);

 return Eta[Inoutn];
}


//####################
// u_h (x/h, y/h)
//####################

double CnoidalClass::u_h(double x, double Y)
{

 double Inoutalpha=alpha;
 double Inoutk=k;
 int Inoutn=n;
 double Inoutm=m;
 double Inoute=e;
 double Inoutdelta=delta;

 int i;
 double	C[11], uu[7], y[11];

//print(`"uu" is U/sqrt(g.h) as a function of delta, m, y/h and cn^2, Fenton (1999, eqn A.3.1)`);
 C[1] = cn(Inoutalpha*x, Inoutk);
 y[1] = Y;

 for(i=2;i<=10;i++)
 {
  C[i] = C[i-1]*C[1];
  y[i] = y[i-1]*y[1];
 }

 Set_mmee(Inoutm,Inoute);

 uu[1] = -1+(1./2-mm[1]+mm[1]*C[2])*Inoutdelta;
 uu[2] = uu[1]+pow(Inoutdelta,2)*(-79./40*mm[2]-19./40+79./40*mm[1]+C[2]*(-3./2*mm[1]+3*mm[2])-mm[2]*C[4]+(-3./4*mm[1]+3./4*mm[2]+C[2]*(-3*mm[2]+3./2*mm[1])+9./4*mm[2]*C[4])*y[2]);
 uu[3] = uu[2]+pow(Inoutdelta,3)*(55./112+7113./1120*mm[2]-2371./560*mm[3]-3471./1120*mm[1]+C[2]*(71./40*mm[1]-339./40*mm[2]+339./40*mm[3])+C[4]*(27./10*mm[2]-27./5*mm[3])+6./5*mm[3]*C[6]+(9./8*mm[1]-27./8*mm[2]+9./4*mm[3]+C[2]*(-27./2*mm[3]-9./4*mm[1]+27./2*mm[2])+C[4]*(-75./8*mm[2]+75./4*mm[3])-15./2*mm[3]*C[6])*y[2]+(-3./8*mm[3]-3./16*mm[1]+9./16*mm[2]+C[2]*(3./8*mm[1]+51./16*mm[3]-51./16*mm[2])+C[4]*(-45./8*mm[3]+45./16*mm[2])+45./16*mm[3]*C[6])*y[4]);
 uu[4] = uu[3]+pow(Inoutdelta,4)*(-11813./22400-382841./28000*mm[2]+108923./5600*mm[3]-108923./11200*mm[4]+31581./8000*mm[1]+C[2]*(-53327./42000*mm[1]+1192733./84000*mm[2]-39177./1120*mm[3]+13059./560*mm[4])+C[4]*(-13109./3000*mm[2]+12793./600*mm[3]-12793./600*mm[4])+C[6]*(-1763./375*mm[3]+3526./375*mm[4])-197./125*mm[4]*C[8]+(1017./160*mm[4]-213./160*mm[1]+123./16*mm[2]-1017./80*mm[3]+C[2]*(5967./80*mm[3]+213./80*mm[1]-483./16*mm[2]-1989./40*mm[4])+C[4]*(3231./160*mm[2]-15579./160*mm[3]+15579./160*mm[4])+C[6]*(729./20*mm[3]-729./10*mm[4])+189./10*mm[4]*C[8])*y[2]+(-27./16*mm[4]+27./8*mm[3]+9./32*mm[1]-63./32*mm[2]+C[2]*(-9./16*mm[1]-999./32*mm[3]+369./32*mm[2]+333./16*mm[4])+C[4]*(453./8*mm[3]-327./32*mm[2]-453./8*mm[4])+C[6]*(-915./32*mm[3]+915./16*mm[4])-315./16*mm[4]*C[8])*y[4]+(-3./160*mm[1]+57./320*mm[2]+51./320*mm[4]-51./160*mm[3]+C[2]*(279./80*mm[3]-93./40*mm[4]+3./80*mm[1]-99./80*mm[2])+C[4]*(567./80*mm[4]-567./80*mm[3]+189./160*mm[2])+C[6]*(-63./8*mm[4]+63./16*mm[3])+189./64*mm[4]*C[8])*y[6]);
 uu[5] = uu[4]+pow(Inoutdelta,5)*(57159./98560+327236467./17248000*mm[2]-884845613./17248000*mm[3]-57144683./2464000*mm[5]+57144683./985600*mm[4]-124831351./34496000*mm[1]+C[2]*(-144821./156800*mm[1]-34543./3136*mm[2]+14639941./196000*mm[3]-3566001./28000*mm[4]+3566001./56000*mm[5])+C[4]*(1131733./294000*mm[2]-26486863./588000*mm[3]+3137133./28000*mm[4]-1045711./14000*mm[5])+C[6]*(757991./73500*mm[3]-72731./1500*mm[4]+72731./1500*mm[5])+C[8]*(298481./36750*mm[4]-298481./18375*mm[5])+13438./6125*mm[5]*C[10]+(-39177./896*mm[4]+39177./2240*mm[5]+53327./56000*mm[1]-1299387./112000*mm[2]+9221./250*mm[3]+C[2]*(-11797957./56000*mm[3]-53327./28000*mm[1]+358171./8000*mm[2]-232269./1400*mm[5]+232269./700*mm[4])+C[4]*(-1628189./56000*mm[2]+29702871./112000*mm[3]+4638023./11200*mm[5]-13914069./22400*mm[4])+C[6]*(-192481./2000*mm[3]-893761./2000*mm[5]+893761./2000*mm[4])+C[8]*(11187./50*mm[5]-11187./100*mm[4])-5319./125*mm[5]*C[10])*y[2]+(1989./128*mm[4]-1989./320*mm[5]-4191./320*mm[3]-213./640*mm[1]+657./160*mm[2]+C[2]*(213./320*mm[1]+9753./80*mm[3]-3075./128*mm[2]+62649./640*mm[5]-62649./320*mm[4])+C[4]*(-139149./640*mm[3]+13563./640*mm[2]-112023./320*mm[5]+336069./640*mm[4])+C[6]*(68643./640*mm[3]+330183./640*mm[5]-330183./640*mm[4])+C[8]*(-5481./16*mm[5]+5481./32*mm[4])+1701./20*mm[5]*C[10])*y[4]+(9./320*mm[1]-387./640*mm[2]-333./128*mm[4]+171./80*mm[3]+333./320*mm[5]+C[2]*(-423./20*mm[5]-4077./160*mm[3]+423./10*mm[4]-9./160*mm[1]+693./160*mm[2])+C[4]*(1461./16*mm[5]-4383./32*mm[4]+54*mm[3]-267./64*mm[2])+C[6]*(-2541./16*mm[5]+2541./16*mm[4]-987./32*mm[3])+C[8]*(7875./64*mm[5]-7875./128*mm[4])-567./16*mm[5]*C[10])*y[6]+(-9./8960*mm[1]+153./4480*mm[2]-279./4480*mm[5]-81./640*mm[3]+279./1792*mm[4]+C[2]*(14769./8960*mm[3]-333./1280*mm[2]+6219./4480*mm[5]-6219./2240*mm[4]+9./4480*mm[1])+C[4]*(4293./448*mm[4]-3321./896*mm[3]-1431./224*mm[5]+459./1792*mm[2])+C[6]*(567./256*mm[3]+2997./256*mm[5]-2997./256*mm[4])+C[8]*(-1215./128*mm[5]+1215./256*mm[4])+729./256*mm[5]*C[10])*y[8]);

 uu[6] =  uu[5];

// NB - the Aitken velocities were a bit irregular, so I did not apply them

 return uu[Inoutn];
}

//####################
// v_h (delta, x, y)
//####################

double CnoidalClass::v_h(double x, double Y)
{

 double Inoutalpha=alpha;
 double Inoutk=k;
 int Inoutn=n;
 double Inoutm=m;
 double Inoute=e;
 double Inoutdelta=delta;

 int		i;
 double	C[11], vv[7], y[11], Lead, S, D;

 C[1] = cn(Inoutalpha*x, Inoutk);
 y[1] = Y;

 for(i=2;i<=10;i++)
 {
  C[i] = C[i-1]*C[1];
  y[i] = y[i-1]*y[1];
 }

 Set_mmee(Inoutm,Inoute);

 S = sn(Inoutalpha*x, Inoutk);
 D = dn(Inoutalpha*x, Inoutk);

 Lead = y[1]*mm[1]*C[1]*S*D*sqrt(3)*pow(Inoutdelta,3./2);

 vv[1] =  1;
 vv[2] = vv[1]+((1./2-mm[1]+(3./2)*mm[1]*C[2])*y[2]-2*mm[1]*C[2]+3*mm[1]-3./2)*Inoutdelta;
 vv[3] = vv[2]+(((27./16)*mm[2]*C[4]+((9./8)*mm[1]-(9./4)*mm[2])*C[2]-(51./80)*mm[1]+(51./80)*mm[2]+3./40)*y[4]+(-(15./2)*mm[2]*C[4]+(-(25./4)*mm[1]+(25./2)*mm[2])*C[2]-3./4-(9./2)*mm[2]+(9./2)*mm[1])*y[2]+(18./5)*mm[2]*C[4]+((27./5)*mm[1]-(54./5)*mm[2])*C[2]+(339./40)*mm[2]+71./40-(339./40)*mm[1])*pow(Inoutdelta,2);
 vv[4] = vv[3]+(((27./16)*mm[3]*C[6]+(-(27./8)*mm[3]+(27./16)*mm[2])*C[4]+(-(81./40)*mm[2]+(81./40)*mm[3]+(27./80)*mm[1])*C[2]-(99./560)*mm[1]+(279./560)*mm[2]+3./560-(93./280)*mm[3])*y[6]+(-(63./4)*mm[3]*C[6]+(-(549./32)*mm[2]+(549./16)*mm[3])*C[4]+(-(327./80)*mm[1]-(453./20)*mm[3]+(453./20)*mm[2])*C[2]-(999./160)*mm[2]+(333./80)*mm[3]-9./80+(369./160)*mm[1])*y[4]+((126./5)*mm[3]*C[6]+((729./20)*mm[2]-(729./10)*mm[3])*C[4]+((5193./80)*mm[3]+(1077./80)*mm[1]-(5193./80)*mm[2])*C[2]+(1989./80)*mm[2]-(161./16)*mm[1]-(663./40)*mm[3]+71./80)*y[2]-(788./125)*mm[3]*C[6]+(-(1763./125)*mm[2]+(3526./125)*mm[3])*C[4]+(-(12793./300)*mm[3]-(13109./1500)*mm[1]+(12793./300)*mm[2])*C[2]+(1192733./84000)*mm[1]-(39177./1120)*mm[2]+(13059./560)*mm[3]-53327./42000)*pow(Inoutdelta,3);
 vv[5] = vv[4]+(((405./256)*mm[4]*C[8]+((135./64)*mm[3]-(135./32)*mm[4])*C[6]+((999./256)*mm[4]+(189./256)*mm[2]-(999./256)*mm[3])*C[4]+(-(369./448)*mm[2]+(477./224)*mm[3]+(51./896)*mm[1]-(159./112)*mm[4])*C[2]+(691./4480)*mm[4]+1./4480-(37./1280)*mm[1]-(691./2240)*mm[3]+(1641./8960)*mm[2])*y[8]+(-(405./16)*mm[4]*C[8]+(-(1125./32)*mm[3]+(1125./16)*mm[4])*C[6]+(-(423./32)*mm[2]+(1089./16)*mm[3]-(1089./16)*mm[4])*C[4]+((108./7)*mm[2]-(4383./112)*mm[3]-(267./224)*mm[1]+(1461./56)*mm[4])*C[2]-9./1120+(99./160)*mm[1]-(423./140)*mm[4]-(4077./1120)*mm[2]+(423./70)*mm[3])*y[6]+((1701./20)*mm[4]*C[8]+((5481./40)*mm[3]-(5481./20)*mm[4])*C[6]+((990549./3200)*mm[4]+(205929./3200)*mm[2]-(990549./3200)*mm[3])*C[4]+((13563./1600)*mm[1]-(112023./800)*mm[4]+(336069./1600)*mm[3]-(139149./1600)*mm[2])*C[2]+(9753./400)*mm[2]-(62649./1600)*mm[3]-(615./128)*mm[1]+(62649./3200)*mm[4]+213./1600)*y[4]+(-(1773./25)*mm[4]*C[8]+(-(3729./25)*mm[3]+(7458./25)*mm[4])*C[6]+((893761./2000)*mm[3]-(893761./2000)*mm[4]-(192481./2000)*mm[2])*C[4]+(-(4638023./11200)*mm[3]+(4638023./16800)*mm[4]-(1628189./84000)*mm[1]+(9900957./56000)*mm[2])*C[2]+(358171./24000)*mm[1]-(11797957./168000)*mm[2]-(77423./1400)*mm[4]-53327./84000+(77423./700)*mm[3])*y[2]+(13438./1225)*mm[4]*C[8]+(-(1193924./18375)*mm[4]+(596962./18375)*mm[3])*C[6]+((757991./24500)*mm[2]-(72731./500)*mm[3]+(72731./500)*mm[4])*C[4]+((3137133./14000)*mm[3]-(26486863./294000)*mm[2]+(1131733./147000)*mm[1]-(1045711./7000)*mm[4])*C[2]-(34543./3136)*mm[1]+(14639941./196000)*mm[2]-(3566001./28000)*mm[3]+(3566001./56000)*mm[4]-144821./156800)*pow(Inoutdelta,4);

 vv[6] =  vv[5];

// NB - the Aitken velocities were a bit irregular, so I did not apply them

 return vv[Inoutn]*Lead;
}

//#############################
// Output for CnoidalSolution
//#############################

void CnoidalClass::Output(FILE *Solution)
{
 fprintf(Solution,"\nHeight/Depth;%6.3f", H);
 iff (Case,Wavelength) fprintf(Solution,"  Length/Depth;%8.3f", L);
 iff (Case,Period) fprintf(Solution,"  Dimensionless Period;%8.3f", T);
 fprintf(Solution,"\nCurrent criterion; %d  Magnitude;%6.3lf\n", Current_criterion,Current);
 fprintf(Solution,"\n# Cnoidal Theory - %s",Method);
 fprintf(Solution,"\n\n# Stokes-Ursell No.;   %7.3f", SU);
 fprintf(Solution,"\n# Elliptic parameter m; %10.7f\n", m);
 fprintf(Solution,"\n# Integral quantities\n");
 fprintf(Solution,"\n# Solution non-dimensionalised by g & mean depth");
 fprintf(Solution,"\n 1 Wave length                       %8.4f", L);
 fprintf(Solution,"\n 2 Wave height                       %8.4f", H);
 fprintf(Solution,"\n 3 Period                            %8.4f", T);
 fprintf(Solution,"\n 4 Wave speed                        %8.4f", c);
 fprintf(Solution,"\n 5 Eulerian current                  %8.4f", ce);
 fprintf(Solution,"\n 6 Stokes current                    %8.4f", cs);
 fprintf(Solution,"\n 7 Mean fluid speed in frame of wave %8.4f", Ubar_d);
 fprintf(Solution,"\n 8 Discharge                         %8.4f", Q_d);
 fprintf(Solution,"\n 9 Bernoulli constant                %8.4f\n", R_d);
 fprintf(Solution,"\n");
}

#endif

