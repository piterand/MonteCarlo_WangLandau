/*
        wl_diluted_ising_triangular_modified.c

        Wang-Landau method for diluted Ising model (triangular lattice)

           programmed by YO          2008/04/07  (original for square lattice)
           programmed by YO          2016/07/13

        nla=nx*ny     : number of lattice sites
        p               : concentration of holes [0,1]

        g[E+3*nla]        : log of DOS
               -3*nla <= E <= nla  (asymmetric because of frustration)
        g[ie]  0 <= ie <= nemax (= 4 * nla)

        n      : number of iteration (0-nfinal)
        nfinal : final number of n   (default 24)
        f      : modification factor for g[ie] (1/2**n)
        factor : criterion of flat   (default 0.80) 

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void period();
void spinset();
void energyset();
void mc();
void single();
void gupdate();
void normalize();

#define nx 24
#define ny nx
#define nla (nx*ny)
#define nemax (4*nla)

  int isp[nla];
  int isp0[nla];                                 // modified
  int nn[6*nla];

  double g[nemax+1];
  double visit[nemax+1];
  double hist[nemax+1];
  int nonzero[nemax+1];

  double f,factor;
  int nfinal;

  int energy;

  int iflag;
  double gmin;

  double p;
  int nhole;

  int pstyle = 0;  // 0 for fixed p 
                   // 1 for average p

int main(void)
{
    int iri1, iri2, ie;

    factor = 0.8;
    nfinal = 24;

    printf("# Wang-Landau calculation of Ising model ");
    printf("on the triangular lattice\n");

    printf("# Input data are random number seeds (iri1,iri2) ");
    printf("and concentration of holes (p) \n");

    scanf("%d %d",&iri1,&iri2);
    printf("# random number seeds %10d %10d\n",iri1,iri2);

    scanf("%lf",&p);

//    srand48(iri1);
    srand(iri1);

    period();
    spinset();
    printf("# linear size of cube     L= %2d,  number of sites  nla= %d\n",
            nx,nla);
    printf("# concentration of holes  %5.3f,  number of holes  %d,  number of spins  %d\n",p,nhole,nla-nhole);
    printf("# actual concentration    %5.3f",(double)nhole/nla);
    if(pstyle==0){
      printf("   p fixed style\n");
    } else {
      printf("   average p style\n");
    }

    energyset();
    printf("# initial energy = %d\n",energy);    // modified

    for(ie=0; ie<=nemax; ie++){
      g[ie]=0;
      hist[ie]=0;
    }

//  srand48(iri2);
    srand(iri2);
    mc();
    normalize();

    for(ie=0; ie<=nemax; ie++){   
      if (nonzero[ie] == 1) {
        printf("%d  %e  %e  %e\n",ie-3*nla,g[ie],g[ie]/(nla-nhole),hist[ie]); 
      }
    }

}

void period()
/*
       periodic boundary conditions for triangular lattice
                   system size = nx*ny*nz
*/
{
    int la, ix;

    for (la=0; la <= nla-1; la++){
      ix=((int)(la/nx))*nx;
      nn[la]       = (la+1)%nx   +ix;
      nn[la+nla]   = (la-1+nx)%nx+ix;
      nn[la+2*nla] = (la+nx) % nla;
      nn[la+3*nla] = (la-nx+nla)% nla;
      nn[la+4*nla] = ((la+1)%nx+ix+nx)% nla;
      nn[la+5*nla] = ((la-1+nx)%nx+ix-nx+nla)% nla;
    }

    for (la=nla-nx; la <= nla-1; la++){
      ix=((int)(la/nx))*nx;
      nn[la+2*nla] = ((la+nx) % nla - ny/2 + nx)%nx;
      nn[la+4*nla] = (((la+1)%nx+ix+nx)% nla - ny/2 + nx)%nx;
    }
    for (la=0; la <= nx-1; la++){
      ix=((int)(la/nx))*nx;
      nn[la+3*nla] = ((la-nx+nla)% nla + ny/2)%nx + nla-nx;
      nn[la+5*nla] = (((la-1+nx)%nx+ix-nx+nla)% nla + ny/2)%nx + nla-nx;
    }
}

void spinset()
/*
        set initial spins
*/
{
    int la, la1;
    int ihole;

    for (la=0; la <= nla-1; la++){
      isp[la]=1;
    }

    if(pstyle==0){
//
//  p fixed
//
      nhole = (nla*p+0.5);
      ihole = 0;
      for (la1=0; la1 <= 50*nla-1; la1++){
        if(ihole>=nhole){break;}
        //la=lrand48()%nla;
        la=rand()%nla;
        if(isp[la]!=0)
        {
          isp[la]=0;
          ihole++;
        }
      }

//      for (la=0; la <= nla-1; la++){
//        printf("%d %d\n",la,isp[la]);
//      }

//      printf("%d %d\n",nhole,ihole);

      if(ihole != nhole){printf("stop\n"); exit(1);}

    } else {
//
//  average p
//
      nhole = 0;
      for (la=0; la <= nla-1; la++){
        if(
                //drand48()<p
                (rand()/RAND_MAX) < p
                ){
          isp[la]=0;
          nhole++;
        }
//      printf("%d %d\n",la,isp[la]);
      }
    }

    for (la=0; la <= nla-1; la++){               // modified
      isp0[la]=isp[la];                          // modified
    }                                            // modified

}

void energyset()
/*
        set initial energy
*/
{
    int la, isp1, i;

    energy=0;
    for (la=0; la <= nla-1; la++){
      isp1=isp[la];
      energy += -isp1*(isp[nn[la]]+isp[nn[la+2*nla]]
                      +isp[nn[la+4*nla]]);
    }
//    printf("# initial energy = %d\n",energy);  // modified
}

void mc()
/*
        monte carlo update
*/
{
  int ie,n;
  int la;                                        // modified
  int check,flag;
  int step, totalstep;
  int count;
  double sum;

  char filename[20];
  FILE *fpw;

/*   initialization  */
  totalstep=0;
  f=1;

  for(ie=0; ie<=nemax; ie++){
    nonzero[ie]=0;
  }

  for( n = 0; n <= nfinal; n++){

    flag=0;
    step=0;

    for(ie=0; ie<=nemax; ie++){ 
      visit[ie]=0;
    }

    while(flag == 0){

      if(n<16 && step%200000==0){                // modified
        for (la=0; la <= nla-1; la++){           // modified
          isp[la]=isp0[la];                      // modified
        }                                        // modified
        energyset();                             // modified
      }                                          // modified

      single();

      step++;

      if(step%1000==0){

        for(ie=0; ie<=nemax; ie++){
          if(visit[ie] > 0) {nonzero[ie]=1;}
        }

        count=0;
        sum=0;
        for(ie=0; ie<=nemax; ie++){
          if(nonzero[ie]==1) {
            count++;
            sum+=visit[ie];
          }
        }

        check=1; 
        for(ie=0; ie<=nemax; ie++){
          if(nonzero[ie]==1) {
            if(visit[ie] < factor*(sum/count)){check=0;}
          }
        }

        if(check==1){flag++;}
      }
    }

    gupdate();

    totalstep += step;

    printf("# n=%2d    MCS=%9d\n",n,totalstep);

    sprintf(filename,"f%d%d.dat",n/10,n%10);
    fpw = fopen(filename,"w+");
    fprintf(fpw,"# n=%2d  step=%9d  MCS=%9d\n",n,step,totalstep);
    for (ie=0; ie<=nemax; ie++){
      if (nonzero[ie] == 1) {
        fprintf(fpw,"%d  %e  %e  %e\n",ie-3*nla,g[ie],visit[ie],hist[ie]);
      }
    }
    fclose(fpw); 

    f = f/2;
  }
  printf("# final   MCS=%9d\n",totalstep);

}

void single()
/*   single spin flip */
{
  int la,la1,isp1,i;
  int energyn;
  double ga,gb;

  for(la1=0; la1 <= nla-1; la1++){
    //la=lrand48()%nla;
      la=rand()%nla;
    isp1 = isp[la];

    if(isp[la]!=0){
      energyn = energy;
      energyn += 2*isp1*(isp[nn[la]]+isp[nn[la+nla]]
                      +isp[nn[la+2*nla]]+isp[nn[la+3*nla]]
                      +isp[nn[la+4*nla]]+isp[nn[la+5*nla]]);

//      if((energyn+3*nla) <= nemax){
        ga = g[(energy+3*nla)];
        gb = g[(energyn+3*nla)];

        if(exp(ga-gb) >
                //drand48()
                rand()/RAND_MAX
                ){
          isp[la] = -isp1;
          energy  = energyn;
        }

        g[(energy+3*nla)]     += f;
        visit[(energy+3*nla)] += 1;
        hist[(energy+3*nla)]  += 1;
//      }
    }
  }
}

void gupdate()
{
  int ie;
  double gmin;

/* set min of g[ie] as 1 */
  gmin=10000000;

  for (ie=0; ie<=nemax; ie++){
    if (nonzero[ie] == 1) {
      if(g[ie] < gmin) {
        gmin = g[ie];
      }
    }
  }

  for (ie=0; ie<=nemax; ie++){
    if (nonzero[ie] == 1) {
      g[ie] += -gmin;
    }
  }
}

void normalize()
{
  int ie;
  double gmax, sum, a;

  gmax = -1000;
  for(ie=0; ie<nemax; ie++){
    if (nonzero[ie] == 1) {
      if(g[ie]>gmax){
        gmax = g[ie];
      }
    }
  }

  sum=0;
  for(ie=0; ie<nemax; ie++){
    if (nonzero[ie] == 1) {
      sum += exp(g[ie]-gmax);
    }
  }

  a = (nla-nhole)*log(2) - gmax - log(sum);

  for(ie=0; ie<nemax; ie++){
    if (nonzero[ie] == 1) {
      g[ie] += a;
    }
  }
}

