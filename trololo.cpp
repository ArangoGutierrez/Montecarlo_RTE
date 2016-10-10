#include <cstdio>
#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
using namespace std;
/*#include <thread>*/

//		Welcome to my
// 	Personal Library for Code Testing
//		Arango Gutierrez

const long double pi = 3.14159265358979323846;// Defines pi

/*
double random_rngint(int tr,double ri,double rs)// random numbers 
{
	time_t now;time(&now);
	long double X= now,M = 2147483648,n,o;
	int a = 1103515245,c = 12345,i;
	for(i=0; i<tr; i++) 
	{
	do 	{
		X = fmod ((a*X+c),M);// Linear congruence Method
		n = (X/M--);
		o =(n*rs);
		}while (o >rs || o <ri);
	}  
	return o;
}

double random_exp(int tr,int la)// random numbers with -Ln|1-y|/la
{
	int i;
	double e;
	for(i=0; i<tr; i++) 
	{
	e = -(log(1-random_int(i)))/(la);
	}  
	return e;
}

double random_radexp(int n)// Exponential radians for 0 to 2*Pi with
{
	int i;
	double te;
	for(i=1; i<(n+1); i++) 
	{
	te = 2*pi*random_exp(i,5);
	}  	
	return te;
}
double random_rad(int n)// Exponential radians for 0 to 2*Pi with
{
	int i;
	double teta;
	for(i=1; i<(n+1); i++) 
	{
	teta = 2.*pi*random_int(i);
	}  	
	return teta;
}
double random_deg(int n,double lbd, double rbd)// degrees from lbd to rbd
{
	long double M = 2147483648, X=1.0;
	int a = 1103515245,c = 12345,i;
	double min=100,max=0,e,rb,lb,d;
	lb = lbd*pi/180.;
	rb = rbd*pi/180.;
	for(i=0; i<n; i++) 
	{
	do{
		X = fmod ((a*X+c),M);
		e = 2.*pi*(X/(M-1.));
		d = ((e)*180/pi);
	} while (e >= rb || e <= lb); 
	}  	
	return e;
}*/
/* Newton Method function for Rayleigh Scattering*//*
long double NM(long double y) 
{	
	int i;
	long double teta,ti=2*pi;
	for(i=1;i<6;i++)
		{
		teta =(ti-((6*pi*(((3*ti+0.5*sin(2*ti))/(6*pi))-y))/(cos(2*ti)+3)));
		ti=teta;
		}
	return teta;//In Radians [0,2*pi)
}*//*
void vctr(int i,int a, vector<int>* v) //Vector passed by reference
{
  v->at(i)=a+i;
} */

/*random numbers Function*/
double random_int()
{	
	clock_t t;
	t=clock();
	double tc = (((float)t)/CLOCKS_PER_SEC);
	long double X = (pi*t)+(pi/tc),M = 2147483648,n;
	int a = 1103515245,c = 12345,i;
	for(i=1;i<5;i++)
	{
		X = fmod ((a*X+c),M);// Linear congruence
		n = (X/(M-1));
		}  
	return n;
}
/* Newton Method function for arcos()*//*
long double arcos(long double l, long double x1,long double x2,vector<long double>* b) 
{	
	double beta,bi=pi/2;
	for(int i=1;i<4;i++)
		{
		beta = bi-(((-(x1-x2)/l)+cos(bi))/(-sin(bi)));
		bi = beta;
		}
	b->at(0)=beta;//In Radians [0,2*pi)
}*/
/*Color RGB spectre*//*
void spectral_color(double &r,double &g,double &b,double l) // RGB <0,1> <- lambda l <400,700> [nm]
    {
    double t;  r=0.0; g=0.0; b=0.0;
         if ((l>=400.0)&&(l<410.0)) { t=(l-400.0)/(410.0-400.0); r=    +(0.33*t)-(0.20*t*t); }
    else if ((l>=410.0)&&(l<475.0)) { t=(l-410.0)/(475.0-410.0); r=0.14         -(0.13*t*t); }
    else if ((l>=545.0)&&(l<595.0)) { t=(l-545.0)/(595.0-545.0); r=    +(1.98*t)-(     t*t); }
    else if ((l>=595.0)&&(l<650.0)) { t=(l-595.0)/(650.0-595.0); r=0.98+(0.06*t)-(0.40*t*t); }
    else if ((l>=650.0)&&(l<700.0)) { t=(l-650.0)/(700.0-650.0); r=0.65-(0.84*t)+(0.20*t*t); }
         if ((l>=415.0)&&(l<475.0)) { t=(l-415.0)/(475.0-415.0); g=             +(0.80*t*t); }
    else if ((l>=475.0)&&(l<590.0)) { t=(l-475.0)/(590.0-475.0); g=0.8 +(0.76*t)-(0.80*t*t); }
    else if ((l>=585.0)&&(l<639.0)) { t=(l-585.0)/(639.0-585.0); g=0.84-(0.84*t)           ; }
         if ((l>=400.0)&&(l<475.0)) { t=(l-400.0)/(475.0-400.0); b=    +(2.20*t)-(1.50*t*t); }
    else if ((l>=475.0)&&(l<560.0)) { t=(l-475.0)/(560.0-475.0); b=0.7 -(     t)+(0.30*t*t); }
    }
*/
int main() 
{	
	FILE *dskw1;
	char FileName[50];
	int file,i;
	double lambda;
  	file=sprintf(FileName,"./OutputData/experimento.dat");
	file++;
  	dskw1=fopen(FileName,"w+");
	/*double nrolls=100000;*/
#pragma omp parallel for private(i)
	for(i=0;i<100000;i++)
	{
	lambda=random_int();
	printf("%f\n",lambda);
	fprintf(dskw1,"%f\n",lambda);
	}
	return 0;
}












