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
	double re=6970948.7;	//Earth center to Thermosphere 6970.9487 km in m
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
double AM0v_to_i(int lambda) // Energy Correction function
{
	double iters; //Number of Iterations 
	int j	=	lambda-380;
float vlambda[322][2] = {	
{380,1202},	{381,1082},	{382,791.3},	{383,684.1},
{384,959.7},	{385,1008},	{386,1007},	{387,1004},
{388,984.3},	{389,1174},	{390,1247},	{391,1342},
{392,1019},	{393,582.3},	{394,1026},	{395,1314},
{396,854.5},	{397,928.8},	{398,1522},	{399,1663},
{400,1682},	{401,1746},	{402,1759},	{403,1684},
{404,1674},	{405,1667},	{406,1589},	{407,1628},
{408,1735},	{409,1715},	{410,1532},	{411,1817},
{412,1789},	{413,1756},	{414,1737},	{415,1734},
{416,1842},	{417,1665},	{418,1684},	{419,1701},
{420,1757},	{421,1797},	{422,1582},	{423,1711},
{424,1767},	{425,1695},	{426,1698},	{427,1569},
{428,1587},	{429,1475},	{430,1135},	{431,1686},
{432,1646},	{433,1731},	{434,1670},	{435,1723},
{436,1929},	{437,1806},	{438,1567},	{439,1825},
{440,1713},	{441,1931},	{442,1980},	{443,1909},
{444,1973},	{445,1821},	{446,1891},	{447,2077},
{448,1973},	{449,2027},	{450,2144},	{451,2109},
{452,1941},	{453,1970},	{454,1979},	{455,2034},
{456,2077},	{457,2100},	{458,1971},	{459,2009},
{460,2040},	{461,2055},	{462,2104},	{463,2040},
{464,1976},	{465,2042},	{466,1921},	{467,2015},
{468,1994},	{469,1990},	{470,1877},	{471,2018},
{472,2041},	{473,1991},	{474,2051},	{475,2016},
{476,1956},	{477,2075},	{478,2009},	{479,2076},
{480,2035},	{481,2090},	{482,2023},	{483,2019},
{484,1969},	{485,1830},	{486,1625},	{487,1830},
{488,1914},	{489,1960},	{490,2007},	{491,1896},
{492,1896},	{493,1888},	{494,2058},	{495,1926},
{496,2017},	{497,2018},	{498,1866},	{499,1970},
{500,1857},	{501,1812},	{502,1894},	{503,1934},
{504,1869},	{505,1993},	{506,1961},	{507,1906},
{508,1919},	{509,1916},	{510,1947},	{511,1997},
{512,1867},	{513,1861},	{514,1874},	{515,1900},
{516,1669},	{517,1726},	{518,1654},	{519,1828},
{520,1831},	{521,1906},	{522,1823},	{523,1894},
{524,1958},	{525,1930},	{526,1674},	{527,1828},
{528,1897},	{529,1918},	{530,1952},	{531,1963},
{532,1770},	{533,1923},	{534,1858},	{535,1990},
{536,1871},	{537,1882},	{538,1904},	{539,1832},
{540,1769},	{541,1881},	{542,1825},	{543,1879},
{544,1879},	{545,1901},	{546,1879},	{547,1833},
{548,1863},	{549,1895},	{550,1862},	{551,1871},
{552,1846},	{553,1882},	{554,1898},	{555,1897},
{556,1821},	{557,1846},	{558,1787},	{559,1808},
{560,1843},	{561,1824},	{562,1850},	{563,1861},
{564,1854},	{565,1798},	{566,1829},	{567,1887},
{568,1810},	{569,1860},	{570,1769},	{571,1823},
{572,1892},	{573,1876},	{574,1867},	{575,1830},
{576,1846},	{577,1857},	{578,1783},	{579,1828},
{580,1838},	{581,1853},	{582,1873},	{583,1857},
{584,1860},	{585,1783},	{586,1830},	{587,1848},
{588,1750},	{589,1612},	{590,1813},	{591,1787},
{592,1808},	{593,1796},	{594,1773},	{595,1782},
{596,1805},	{597,1780},	{598,1757},	{599,1774},
{600,1746},	{601,1751},	{602,1719},	{603,1787},
{604,1776},	{605,1763},	{606,1759},	{607,1757},
{608,1743},	{609,1744},	{610,1703},	{611,1746},
{612,1705},	{613,1683},	{614,1713},	{615,1713},
{616,1609},	{617,1707},	{618,1724},	{619,1707},
{620,1734},	{621,1690},	{622,1713},	{623,1666},
{624,1656},	{625,1632},	{626,1697},	{627,1697},
{628,1697},	{629,1677},	{630,1658},	{631,1639},
{632,1645},	{633,1651},	{634,1653.5},	{635,1656},
{636,1655},	{637,1654},	{638,1652.5},	{639,1651},
{640,1632.5},	{641,1614},	{642,1617.5},	{643,1621},
{644,1624},	{645,1627},	{646,1615},	{647,1603},
{648,1580.5},	{649,1558},	{650,1582},	{651,1606},
{652,1602.5},	{653,1599},	{654,1565.5},	{655,1532},
{656,1458},	{657,1384},	{658,1466.5},	{659,1549},
{660,1560},	{661,1571},	{662,1563},	{663,1555},
{664,1557.5},	{665,1560},	{666,1547.5},	{667,1535},
{668,1540.5},	{669,1546},	{670,1531},	{671,1516},
{672,1518.5},	{673,1521},	{674,1515.5},	{675,1510},
{676,1509},	{677,1508},	{678,1503},	{679,1498},
{680,1495},	{681,1492},	{682,1485.5},	{683,1479},
{684,1467},	{685,1455},	{686,1461},	{687,1467},
{688,1464},	{689,1461},	{690,1454.5},	{691,1448},
{692,1448},	{693,1448},	{694,1442},	{695,1436},
{696,1426},	{697,1416},	{698,1420.5},	{699,1425},
{700,1405.5},	{701,1386}
};
iters=vlambda[j][1];
return iters;
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
	iters=AM0v_to_iters(lambda);
	double teta,iters;
	/*double nrolls=1000000;*/
	for(i=0;i<iters;i++)
	{
	vector<double> brmt (1);
	BM(lambda,SunPos, &brmt);
	teta=brmt[0];
	if(teta>=0)	
	fprintf(dskw1,"%d %lf\n",lambda,teta*180/pi);
	}
	return 0;
}
