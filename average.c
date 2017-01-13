#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nlamax (16*8*8*8)

int    iener[4*nlamax+1];
double g[4*nlamax+1];

int main(void){
  int nx, nla;
  double x,beta;
  int ix,itmin,itmax;
  double tmin,tmax;

  int ie;
  double gmax;
  double sum,a;

  int count,iflag;
  double a0,ae,ae2,aecp,cv,as,as2;  //edited
  double weight;
  double fdummy;
  int idummy;

  char * fname1;
  char cdummy[80];
  FILE *fpw;

  scanf("%lf %lf",&tmin,&tmax);

  fname1 = "measure.dat";

  fpw=fopen(fname1,"r");

  fscanf(fpw,"# linear size of cube     L= %d,  number of sites  N= %d\n",
         &idummy,&nla);
  printf("# number of sites  %d\n",nla);

  fscanf(fpw,"# concentration of holes %lf,  number of holes %d\n",
         &fdummy,&idummy);
  nla = nla - idummy;
  printf("# number of spins  %d\n",nla);

  for (ie=1; ie<=6; ie++){ //edited, changed 5 to 6 value
    fscanf(fpw,"%[^\n]%*c",cdummy);
  }

  ie=0;
  while(fscanf(fpw,"%d %lf %lf %lf",&iener[ie],&g[ie],&fdummy,&fdummy) != EOF)
    {
if(ie<1){
      printf("# %d %e\n",iener[ie],g[ie]);
}
iener[ie]*=-1;
      ie++;
    }

  fclose(fpw);

  count=ie;

/*
  gmax=-1000;
  for(ie=0; ie<count; ie++){
    if(g[ie]>gmax){
      gmax=g[ie];
    }
  }
  sum=0;
  for(ie=0; ie<count; ie++){
    sum += exp(g[ie]-gmax);
  }
  a = nla*log(2) - gmax - log(sum);
  for(ie=0; ie<count; ie++){
    g[ie]+=a;
  }

// printf("# %f %f\n",gmax,a);
*/

  itmin=tmin*20;
  itmax=tmax*20;

  for(ix=itmin; ix<=itmax; ix++){
    x=0.05*ix;
    beta=1/x;

    gmax=-1000;

    for(ie=0; ie<count; ie++){
      if(g[ie]-beta*(iener[ie]-iener[0])>gmax){
        gmax=g[ie]-beta*(iener[ie]-iener[0]);
      }
    }

    a0=0;
    ae=0;
    ae2=0;
    aecp=0; //edited
    as=0;
    as2=0;  //edited

    for(ie=0; ie<count; ie++){
        weight=exp(g[ie]-beta*(iener[ie]-iener[0])-gmax);
        a0  += weight;
        ae  += weight*iener[ie];
        ae2 += weight*iener[ie]*iener[ie];
        as  += weight*g[ie];
    }

    aecp=ae;  //edited

    ae  /= nla;
    ae2 /= (double)nla*(double)nla;
    as  /= nla;

    ae  /= a0;
    ae2 /= a0;
    aecp /= a0;  //edited
    as  /= a0;

    as2 = (log(a0) + gmax - beta*iener[0] + aecp*beta)/nla; //edited

    cv = beta*beta*(ae2-ae*ae)*nla;

    printf("%f %e %e %e %e\n",x,ae,cv,as,as2); //edited

  }
}
