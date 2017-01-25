
#ifndef FentonElliptic_H
#define FentonElliptic_H

//	Evaluation of elliptic integrals and functions, via theta functions.

//	NOTE:  The program for theta_3 uses a different approach.  The others are
//	working satisfactorily too, however

//	Used as: K(m), E(m), Z(u,m),sn(x,m), cn(x,m), dn(x,m), and
//				theta_1(x,q) etc. if desired.

#include	<math.h>

#include "FentonHeaders.h"

double theta_1(double, double);
double theta_2(double, double);
double theta_3(double, double);
double theta_4(double, double);

double K(double m)
{
 int i,j=1;
 double cont,e,Kd,m4,q,t;

 if ( m < 0.5 )
	{
	 m = 1.-m ;
	 j = -1;
	}
 m4 = pow(m,0.25);
 e = 0.5 * (1.-m4)/(1.+m4) ;
 q = e + 2.*pow(e,5.) + 15.*pow(e,9.) + 150.*pow(e,13.) ;
 cont = t = 1.;
 for ( i=1 ; fabs(cont) > 1.e-18 ; i++ )
	{
	 cont = 2. * pow(q, (float)(i*i));
	 t = t + cont ;
	}
 Kd = 0.5*Fenton_pi*t*t ;
 e=(Kd/Fenton_pi)*log(1./q) ;
 if( j > 0 ) return(e);
 else return(Kd);
}

double E(double m)
{
 int i,j=1;
 double cont,e,k,m4,q,t;

 if ( m > 0.5 )
	{
	 m = 1.-m ;
	 j = -1;
	}
 m4 = pow(1.-m,0.25);
 e = 0.5 * (1.-m4)/(1.+m4) ;
 q = e + 2.*pow(e,5.) + 15.*pow(e,9.) + 150.*pow(e,13.) ;
 cont = t = 1./24.;
 for ( i=1 ; fabs(cont) > 1.e-18 ; i++ )
	{
	 e = pow(q,2.*i);
	 cont = e/pow(1.-e,2.);
	 t = t - cont ;
	}
 k = K(m);
 e = k*((2.-m)/3.+2.*pow(Fenton_pi/k,2.)*t);
 if( j > 0 )
	return(e);
 else
	{
	 q = K(1.-m);
	 e = 0.5*Fenton_pi/k + q - e*q/k;
	 return(e);
	}
}

//**********************************************
// Elliptic functions
//**********************************************

// Bring back to range [-2K,+2K]

double Shift(double u, double K)
{
 int N;
 double uu;
 N = trunc(u/(4*K));
 uu = u-N*4*K;
 if (uu < -2*K) uu+=4*K;
 if (uu > 2*K) uu -= 4*K;
 return uu;
}

double sn(double u, double m)
{
 double k,kd,term,z,q;

 k = K(m);
 kd = K(1.-m);
 q = exp(-Fenton_pi*kd/k);
 z = 0.5*Fenton_pi*Shift(u,k)/k;
 term = pow(m,-0.25)*theta_1(z,q)/theta_4(z,q);
 return(term);
}

double cn(double u, double m)
{
 double k,kd,q,term,z;

 k = K(m);
 kd = K(1.-m);
 q = exp(-Fenton_pi*kd/k);
 z = 0.5*Fenton_pi*Shift(u,k)/k;
 term = pow((1.-m)/m,0.25)*theta_2(z,q)/theta_4(z,q);
 return(term);
}

double dn(double u, double m)
{
 double k,kd,q,term,z;

 k = K(m);
 kd = K(1.-m);
 q = exp(-Fenton_pi*kd/k);
 z = 0.5*Fenton_pi*Shift(u,k)/k;
 term = pow((1.-m),0.25)*theta_3(z,q)/theta_4(z,q);
 return(term);
}

double theta_1(double z,double q)
{
 int i;
 double K_on_KD,qd,sign,sum,term,factor;

 factor = -1.;
 i = floor(z/Fenton_pi);
 z = z - Fenton_pi*i;
 if ( 2*(i/2) == i ) factor = 1.;
 K_on_KD = -Fenton_pi/log(q);
 qd = exp(-Fenton_pi*K_on_KD);

 if ( q < 0.04321 )
	{
 	 sum = 2.*pow(q,0.25)*sin(z);
	 term = 1.;
	 sign = -1.;
	 for ( i=1 ; fabs(term) > 1.e-16*fabs(sum) ; i++ )
		{
		 term = 2.*sign*pow(q,0.25+i+i*i)*sin(z*(i+i+1));
		 sum = sum + term ;
		 sign = -sign ;
		}
	}
 else
	{
	 sum = 0.;
	 term = 1.;
	 sign = 1.;
	 for ( i=0 ; fabs(term) > 1.e-16*fabs(sum) ; i++ )
		{
		 term = sign*pow(qd,0.25+i+i*i)*sinh(z*(i+i+1)*K_on_KD);
		 sum = sum + term ;
		 sign = -sign ;
		}
	 sum = sum * 2. * sqrt(K_on_KD) * exp(-z*z*K_on_KD/Fenton_pi);
	}
 return(sum*factor);
}

double theta_2(double z,double q)
{
 int i;
 double K_on_KD,qd,sign,sum,term,factor;

 factor = -1.;
 i = floor(z/Fenton_pi);
 z = z - Fenton_pi*i;
 if ( 2*(i/2) == i ) factor = 1.;
 K_on_KD = -Fenton_pi/log(q);
 qd = exp(-Fenton_pi*K_on_KD);

 if ( q < 0.04321 )
	{
	 sum = 2.*pow(q,0.25)*cos(z);
	 term = 1.;
	 for ( i=1 ; fabs(term) > 1.e-16*fabs(sum) ; i++ )
		{
		 term = 2.*pow(q,0.25+i+i*i)*cos(z*(i+i+1));
		 sum = sum + term ;
		}
	}
 else
	{
	 sum = 1.;
	 term = 1.;
	 sign = -1.;
	 for ( i=1 ; fabs(term) > 1.e-16*fabs(sum) ; i++ )
		{
		 term = 2.*sign*pow(qd,(float)i*i)*cosh(z*(i+i)*K_on_KD);
		 sum = sum + term ;
		 sign = -sign ;
		}
 	 sum = sum * sqrt(K_on_KD) * exp(-z*z*K_on_KD/Fenton_pi);
	}
 return(sum*factor);
}

double theta_3(double z,double q)
{
 int i, i1, i2, ia, ib;
 double K_on_KD, sum,t1,term;

 K_on_KD = -Fenton_pi/log(q);

 if ( q < 0.04321 )
	{
 	 sum = 1.;
	 for ( i=1 ; i<=4 ; i++ )
		{
		 term = 2.*pow(q,(float)i*i)*cos(z*(i+i));
		 sum = sum + term ;
		}
	}
 else
	{
	 t1 = sqrt(log(1.e-16)/(-Fenton_pi*K_on_KD));
	 i1 = t1 + z / Fenton_pi;
	 i2 = -t1 - z / Fenton_pi;
	 if ( i1 >= i2 ) {ia = i2-1; ib = i1+1;}
	 else { ia = i1-1 ; ib = i2 + 1;}
	 sum = 0.;
	 for ( i=ia; i <= ib ; i++ )
		{
		 t1 = i - z / Fenton_pi ;
		 term = exp(-Fenton_pi*K_on_KD*t1*t1);
		 sum = sum + term ;
		}
 	 sum = sum * sqrt(K_on_KD);
	}

 return(sum);
}

double theta_4(double z,double q)
{
 int i;
 double K_on_KD,qd,sign,sum,term;

 i = floor(z/Fenton_pi);
 z = z - Fenton_pi*i;
 K_on_KD = -Fenton_pi/log(q);
 qd = exp(-Fenton_pi*K_on_KD);

 sum = 1.;
 if ( q < 0.04321 )
	{
	 sign = -1.;
	 for ( i=1 ; i<=4 ; i++ )
		{
		 term = 2.*sign*pow(q,(float)i*i)*cos(z*(i+i));
		 sum = sum + term ;
		 sign = -sign ;
		}
	}
 else
	{
	 sum = 0.;
	 term = 1.;
	 for ( i=0 ; fabs(term) > 1.e-16*fabs(sum) ; i++ )
		{
	 	 term = pow(qd,0.25+i+i*i)*cosh(z*(i+i+1)*K_on_KD);
		 sum = sum + term ;
		}
	 sum = sum * 2. * sqrt(K_on_KD) * exp(-z*z*K_on_KD/Fenton_pi);
	}
 return(sum);
}

double Z(double u, double m)
{
 int i;
 double k,kd,q,qd,sign,sum,term,t,z;

 k = K(m);
 kd = K(1.-m);
 q = exp(-Fenton_pi*kd/k);
 qd = exp(-Fenton_pi*k/kd);
 z = .5 * Fenton_pi * u / k;
 i = floor(z/Fenton_pi);
 z = z - Fenton_pi*i;
 t = theta_4(z,q);
 sum = 0.;
 if ( m <= 0.5 )
	{
	 sign = 1.;
	 for ( i=1 ; i<=4 ; i++ )
		{
		 term = 4.*i*sign*pow(q,(float)i*i)*sin(z*(i+i));
		 sum = sum + term ;
		 sign = -sign ;
		}
	}
 else
	{
	 sum = 0.;
	 term = 1.;
	 for ( i=0 ; fabs(term) > 1.e-16*fabs(sum) ; i++ )
		{
		 term = (i+i+1)*pow(qd,0.25+i+i*i)*sinh(z*(i+i+1)*k/kd);
		 sum = sum + term ;
		}
	 sum = sum * 2. * sqrt(k/kd) * exp(-z*z*k/(kd*Fenton_pi));
	 sum = (sum - 2. * z * t/Fenton_pi)*k/kd;
	}
 sum = .5 * Fenton_pi * sum /( k * t ) ;
 return(sum);
}

#endif
