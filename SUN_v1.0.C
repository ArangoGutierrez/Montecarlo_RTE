#include<stdio.h>
#include<math.h>

#define wl_m 200//188.153 
#define wl_M 2500//1100.726
#define wl_s 0.45
double hc_KT_nm = 1e3*6.62606957*2.99792458/(1.3806488*5.525);
double hc2_nm   = 6.62606957*2.99792458*2.99792458;

double SUN_Radiance(double wl){
  return 2.e17*hc2_nm/(pow(wl,5)*(exp(hc_KT_nm/wl)-1.));
  
}

int main(){
  FILE *dskw;
  dskw=fopen("SUN_Radiance.dat","w+");
  for(double wl=wl_m;wl<=wl_M;wl+=wl_s)
    fprintf(dskw,"%g %g\n",wl,SUN_Radiance(wl));
  fclose(dskw);
  return 0;
}
