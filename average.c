#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nlamax (16*8*8*8*8*8) // !не знаю какой нужно задать размер?

double    iener[4*nlamax+1];
double g[4*nlamax+1];

int main(void){
    int nla;
    double x,beta;
    int ix,itmin,itmax;
    double tmin,tmax;

    int ie;
    double gmax;
    double sum,a;

    int count;
    double a0,ae,ae2,aecp,cv,as,as2;  //edited
    double weight;
    double fdummy;

    char cdummy[10000];
    FILE *fpw;

    char filename[100];
    printf("# Please, input target filename: ");
    if( scanf("%s",filename)==1){}
    else{printf("Error! Can not find the file!");return 0;}


    fpw=fopen(filename,"r");
    if(fpw == 0){
        printf("Error with open file");
        return 0;
    }


    printf("# Enter min and max temperature \n");
    if(scanf("%lf %lf",&tmin,&tmax)==2){}
    else{printf("Error! Failed to read min and max temperature!"); return 0;}






    printf("# Enter number of spins  \n");
    if(scanf("%d",&nla)==1){}
    else{printf("Error! Failed to read number of spins!"); return 0;}


    char c;

    for(c=fgetc(fpw);c=='\n'||c=='\r'||c=='#';c=fgetc(fpw)) //пропуск комментариев
        fgets(cdummy,10000,fpw);

    fseek(fpw,-1,SEEK_CUR);       // сдвиг курсора на один символ назад

    ie=0;
// !обязательно! первые две колонки должны быть энергия и g[E], остальное неважно
    while(fscanf(fpw,"%lf %lf",&iener[ie],&g[ie]) != EOF)
    {
        c=fgetc(fpw);
        if(c=='\n'||c=='\r'){}
            else
            fgets(cdummy,10000,fpw);

        if(ie<1){
            printf("# %lf %e\n",iener[ie],g[ie]);
        }
        ie++;
    }

    fclose(fpw);

    count=ie;


    itmin=tmin*1000;
    itmax=tmax*1000;

    for(ix=itmin; ix<=itmax; ix++){
        x=0.001*ix;
        beta=1./x;

        gmax=-1000.;

        for(ie=0; ie<count; ie++){
            if(g[ie]-beta*(iener[ie]-iener[0])>gmax){
                gmax=g[ie]-beta*(iener[ie]-iener[0]);
            }
        }

        a0=0;
        ae=0;
        ae2=0;
        aecp=0; //edited

        as2=0;  //edited

        for(ie=0; ie<count; ie++){
            weight=exp(g[ie]-beta*(iener[ie]-iener[0])-gmax);
            a0  += weight;
            ae  += weight*iener[ie];
            ae2 += weight*iener[ie]*iener[ie];

        }

        aecp=ae;  //edited

        ae  /= (double)nla;
        ae2 /= (double)nla*(double)nla;


        ae  /= a0;
        ae2 /= a0;
        aecp /= a0;  //edited


        as2 = (log(a0) + gmax - beta*iener[0] + aecp*beta)/(double)nla; //edited

        cv = beta*beta*(ae2-ae*ae)*(double)nla;

        printf("%f %e %e %e \n",x,ae,cv,as2); //edited

    }
}
