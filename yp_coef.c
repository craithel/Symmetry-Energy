/*---------------------------------------------------------------------------
 * Written by Carolyn Raithel.
 * May 17, 2019.
 * For citations, please reference Raithel & Ozel (2019): arXiv...
 * These routines provide the coefficients a,b,and c for the series
 * expansion of the neutron excess term, defined in Raithel & Ozel (2019):
 *
 *		   (1-2Y_p)^2 = a + b u + c u^2 + ...
 *
 * The terms are derived in the accompanying Mathematica notebook Yp_expansion_public.nb
 *
 * To compile: gcc yp_coefs.c -lm
---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define n0 0.16		// nuclear saturation density, in fm^-3
#define hc 197.327	// hbar*c, in MeV fm

double aa(double S0)
{
	return  pow(pow(2,0.6666666666666666)*pow(hc,2)*n0*pow(M_PI,1.3333333333333333) - 
     		pow(2,0.3333333333333333)*hc*pow(M_PI,0.6666666666666666)*pow(24*n0 + n0*sqrt(576 + 
		(2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)),0.6666666666666666)*S0,2)/
   		(256.*pow(24*n0 + n0*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)),0.6666666666666666)*pow(S0,4));
}

double bb(double S0, double L0)
{
	return  (pow(hc,2)*pow(M_PI,1.3333333333333333)*(-(L0*n0) + S0)*
		(pow(2,0.3333333333333333)*hc*n0*pow(M_PI,0.6666666666666666) - 
       		pow(n0*(24 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))),0.6666666666666666)*S0)*
     		(pow(2,0.3333333333333333)*pow(hc,4)*pow(n0,1.3333333333333333)*pow(M_PI,2.6666666666666665)*
		(32 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))) - 
       		pow(hc,3)*n0*pow(M_PI,2)*(16 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)))*
		pow(24 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)),0.6666666666666666)*S0 + 
       		384*pow(2,0.3333333333333333)*hc*pow(n0,0.3333333333333333)*pow(M_PI,0.6666666666666666)*
		(24 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)))*pow(S0,3) - 
       		192*pow(24 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)),1.6666666666666667)*pow(S0,4)))/
   		(128.*pow(2,0.3333333333333333)*pow(24 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/
		pow(S0,3)),1.6666666666666667)*pow(S0,5)*(pow(hc,3)*n0*pow(M_PI,2) + 288*pow(S0,3)));
}

double cc(double S0, double L0, double Ksym)
{
	return  (pow(hc,2)*pow(M_PI,1.3333333333333333)*(3*pow(S0,3)*pow(-(L0*n0) + S0,2)*
        	pow(pow(2,0.3333333333333333)*pow(hc,4)*pow(n0,2)*pow(M_PI,2.6666666666666665)*
		(32 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))) - 
          	pow(hc,3)*n0*pow(M_PI,2)*(16 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/
		pow(S0,3)))*pow(24*n0 + n0*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)),0.6666666666666666)*
           	S0 + 384*pow(2,0.3333333333333333)*hc*n0*pow(M_PI,0.6666666666666666)*
		(24 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)))*pow(S0,3) - 
          	192*pow(n0,0.6666666666666666)*pow(24 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/
		pow(S0,3)),1.6666666666666667)*pow(S0,4),2) - 
       		2*(pow(2,0.3333333333333333)*hc*n0*pow(M_PI,0.6666666666666666) - 
		pow(n0*(24 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))),0.6666666666666666)*S0)*
        	(5308416*pow(2,0.3333333333333333)*hc*n0*pow(M_PI,0.6666666666666666)*
		(24 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)))*pow(S0,9)*
           	(-3*pow(L0,2)*pow(n0,2) + n0*(4*L0 + Ksym*n0)*S0 + pow(S0,2)) - 
          	2654208*pow(n0,0.6666666666666666)*pow(24 + sqrt(576 + (2*pow(hc,3)*
		n0*pow(M_PI,2))/pow(S0,3)),1.6666666666666667)*pow(S0,10)*
           	(-2*pow(L0,2)*pow(n0,2) + n0*(2*L0 + Ksym*n0)*S0 + 2*pow(S0,2)) + 
          	pow(2,0.3333333333333333)*pow(hc,10)*pow(n0,4)*pow(M_PI,6.666666666666667)*
		(-5*pow(L0,2)*pow(n0,2) + 2*n0*(3*L0 + Ksym*n0)*S0 + 3*pow(S0,2)) - 
          	pow(hc,9)*pow(n0,3)*pow(M_PI,6)*pow(24*n0 + n0*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/
		pow(S0,3)),0.6666666666666666)*S0*
           	(-5*pow(L0,2)*pow(n0,2) + 2*n0*(3*L0 + Ksym*n0)*S0 + 3*pow(S0,2)) - 
          	2304*pow(hc,3)*n0*pow(M_PI,2)*pow(24*n0 + n0*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/
		pow(S0,3)),0.6666666666666666)*pow(S0,7)*
           	(-2*pow(L0,2)*pow(n0,2)*(312 + 11*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))) + 720*L0*n0*S0 + 
             	3*Ksym*pow(n0,2)*(88 + 3*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)))*S0 + 
		26*L0*n0*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))*S0 + 432*pow(S0,2) + 
             	14*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))*pow(S0,2)) - 
          	4*pow(hc,6)*pow(n0,2)*pow(M_PI,4)*pow(24*n0 + n0*sqrt(576 + 
		(2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)),0.6666666666666666)*pow(S0,4)*
           	(-(pow(L0,2)*pow(n0,2)*(1224 + 25*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)))) +
		1488*L0*n0*S0 + 10*Ksym*pow(n0,2)*(48 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)))*S0 +
		30*L0*n0*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))*S0 + 696*pow(S0,2) + 
             	15*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))*pow(S0,2)) + 
          	2304*pow(2,0.3333333333333333)*pow(hc,4)*pow(n0,2)*pow(M_PI,2.6666666666666665)*pow(S0,6)*
           	(-8*pow(L0,2)*pow(n0,2)*(156 + 5*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))) +
		1584*L0*n0*S0 + 3*Ksym*pow(n0,2)*(152 + 5*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)))*S0 + 
		50*L0*n0*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))*S0 + 576*pow(S0,2) + 
             	20*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))*pow(S0,2)) + 
          	4*pow(2,0.3333333333333333)*pow(hc,7)*pow(n0,3)*pow(M_PI,4.666666666666667)*pow(S0,3)*
           	(-(pow(L0,2)*pow(n0,2)*(1704 + 35*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)))) + 2064*L0*n0*S0 + 
             	14*Ksym*pow(n0,2)*(48 + sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)))*S0 +
		42*L0*n0*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))*S0 + 984*pow(S0,2) + 
             	21*sqrt(576 + (2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3))*pow(S0,2)))))/
   		(1536.*pow(2,0.3333333333333333)*pow(n0,0.6666666666666666)*pow(24 + sqrt(576 + 
		(2*pow(hc,3)*n0*pow(M_PI,2))/pow(S0,3)),2.6666666666666665)*pow(S0,9)*
     		pow(pow(hc,3)*n0*pow(M_PI,2) + 288*pow(S0,3),2));
}




int main( int argc, char *argv[])
{

	printf("a=%f \n", aa(32.));			// calculate a for S0=32 MeV
	printf("b=%f \n", bb(32., 40.));		// calculate b for S0=32 MeV, L0=40 MeV
	printf("c=%f \n", cc(32., 40., -110.));		// calculate b for S0=32 MeV, L0=40 MeV, Ksym=-110 MeV

	return 0;
}

