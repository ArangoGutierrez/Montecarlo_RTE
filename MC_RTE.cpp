/*
	Monte-Carlo solution for Radiative transfer equation

		Arango,C;Arguelles,A;Reina,H 
		CIBioFi-QuanTIC
		  
		Last updated: June 28, 2016

	This program is in development. HTC settup

	Centro de Investigación e Innovación en Bioinformática y Fotónica
	https://http://cibiofi.univalle.edu.co/

*/
#include <cstdio>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
using namespace std;

#define pi 3.14159265358979323846
#define gr 980.6160			//Gravity in cm/s
#define Av 6.022140874*pow(10,23)	//Avogadro's Number

/* random number Function Generator by using a linear congruential generator (LCG)
is an algorithm that yields a sequence of pseudo-randomized numbers calculated
with a discontinuous piecewise linear equation.	*/
double random_eng()
{	
	clock_t t;
	t=clock();
	double tc;
	tc = (((float)t)/CLOCKS_PER_SEC);
	long double X = fmod ((pi*t)+(pi/tc),(pi/tc));/*Seed*/
	double M = 2147483648,n;
	int a = 1103515245,c=12345,i;
	for(i=1;i<5;i++)
	{
		X = fmod ((a*X+c),M);// Linear congruence
		n = (X/(M-1));
		}  
	return n;
}
/* Newton Method function for Rayleigh Scattering angle*/
long double NM(long double y) 
{	
	int i;
	long double teta,ti=2*pi;
	for(i=1;i<5;i++)
		{
		teta =(ti-((6*pi*(((3*ti+0.5*sin(2*ti))/(6*pi))-y))/(cos(2*ti)+3)));
		ti=teta;
		}
	return teta;//In Radians [0,2*pi)
}
/*Scattering calculations*/
void SC(long double lambda, vector<double>* tao) //Vector passed by reference
{
	double P=101325.0;//1 Atm aprox in pa
	long double sigma;//Rayleigh Cross-Section
	long double nv;//Refractive index
	long double Fk;//King Factor for depolarization
	long double R=8.3144598,T=288.15;
  	long double Ns = (2.546899*pow(10,19));//Molecular density Molecules Bodhaine et al 1999
  	Fk=1.034+((3.17*pow(10.0,-4.0))*(1/pow(lambda*0.001,2)));//king factor for a 100% N2 atmosphere
	nv=((8060.51+2480990/(132.274-pow(lambda*0.001,-2))+17455.7/(39.32957-pow(lambda*0.001,-2)))/pow(10,8))+1;
	sigma=((24*pow(pi,3))/(pow(lambda*0.0000001,4)*pow(Ns,2)))*(pow((pow(nv,2)-1),2)/pow((pow(nv,2)+2),2))*Fk;
	tao->at(0)=(sigma*0.0001*P*Av)/(R*T);
	//printf("tao: %Lf \n",1/((sigma*0.0001*P*Av)/(R*T)));	//Optical Depth
}
/*Brownian motion for a single particle*/
void BM(int lambda,double SunPos, vector<double>* brmt) //Vector passed by reference
{		
	double xo,yo,xa,ya,xb,yb,l,vx,vy,dg,y,x,m,te,b,ltao,tsd;
	double teta;
	int f,t,i;
	double re=16367444.7;	//Earth center to Exosphere 16377.4447Km in m
	double rt=6377444.7;	//Earth center to Troposphere 6377.4447Km in m
	double rsea=6367444.7;	//Earth center to Sea lvl  6367.4447km in m
	double rcali=6368444.7;	//Earth center to Cali lvl  6368.4447km in m
	//Photons (Xo,Yo) in Troposphere
	dg=(1.531877998)+SunPos*0.003243194083333;	//Astronomical Horizon
	xo=re*cos(dg);
	yo=re*sin(dg);
	xa=xo;
	ya=yo;
	//Zenith Dependent with a 2 degree range
	tsd=0.130899694;/*Degree step*/
	te = ((SunPos*tsd)+pi)+(0.017453293*(random_eng())*pow(-1,int (random_eng()*10)));	
	vector<double> tao (1);
	SC(lambda, &tao);
	ltao=tao[0];
	l = (-1.0/ltao)*log(1.0-random_eng());
	vx = l*cos(te);
	vy = l*sin(te);
	//Browniam movement
	do{
	xb=xa;
	yb=ya;
	xa += vx;
	ya += vy;
	te = fmod (te + NM(random_eng()),(2*pi));
	l = (-1.0/ltao)*log(1.0-random_eng());
	vx = l*cos(te);
	vy = l*sin(te);
	}while(((xa*xa)+(ya*ya)<(pow(re,2))) && (ya>=rcali) && ((pow(xa,2)+pow(ya,2))>=(pow(rt,2))));	
		if(((pow(xa,2)+pow((ya),2))<=(pow(rt,2))))
		{	
		m=(-ya+yb)/(-xa+xb);
		b=ya-(m*xa);
		y=((2*b/pow(m,2))+sqrt(pow((2*b)/pow(m,2),2)-(4*(1+(1/pow(m,2)))*(-pow(rt,2)+(pow(b,2)/pow(m,2))))))/(2*(1+(1/pow(m,2))));
		x=(y-b)/m;
		teta = atan2(y,x);
		brmt->at(0)=teta;
		}
		else brmt->at(0)=-1;
}
/*MonteCarlo*/
int main(int argc, char *argv[]) 
{	
	int nthreads, tid,i;	

	int lambda=atoi(argv[1]);
	int SunPos=atoi(argv[2]);
	FILE *dskw1;
	char FileName[50];
	int file;
  	file=sprintf(FileName,"./OutputData/MC_RTE_v1.2_%d_%d.dat",SunPos,lambda);
	file++;
  	dskw1=fopen(FileName,"w+");
	double teta;
	/*double nrolls=1000000;*/
	for(i=0;i<1000;i++)
	{
	vector<double> brmt (1);
	BM(lambda,SunPos, &brmt);
	teta=brmt[0];
	if(teta>=0)	
	fprintf(dskw1,"%d %lf\n",lambda,teta*180/pi);
	}
	return 0;
}
