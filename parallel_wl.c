/**
        parallel_wl.c

        Wang-Landau method for different magnetic systems

           programmed by:
            Makarov Aleksandr
            Andriushchenko Petr
            Shevchenko Yuriy


*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "mpi.h"


unsigned n;                     // количество спинов
signed char *spins;             //массив направления спинов. По умолчанию +1. n - число считанных спинов (число незакомментированных строк в csv-файле).
unsigned short *a_neighbours;   //число соседей каждого спина. Считается как число энергий в соответствующей строке в csv-файле.
unsigned short *neighbours;     // соседи каждого спина
unsigned int *sequencies;       //для каждого спина описывает, с какого ключа в массиве energies[] начинают описываться парные энергии
double *energies;               //сами энергии из файла. Описывается как одномерный массив. Длина массива - число парных энергий в csv-файле.
double *intervals;              //массив интервалов
double *intervalsE;             //массив интервалов значений
int intervalsNum=0;             //число значений интервалов из файла
double emin, emax;              //минимумы и максимумы энергии
double e;                       //текущая энергия системы
unsigned eCount=0;              //число пар энергий
unsigned histSize=0;            //число элементов в гистограммах

double *g;
unsigned *visit;
unsigned *hist;
int *nonzero;

double f;                       // Модификационный фактор (уменьшается с каждым WL-шагом)
double factor;                  // Критерий плоскости гистограммы H
unsigned nfinal;                // число WL-циклов

#define PRECISION 1e2             // Точность 1eX, где X - Сколько знаков учитывать в энергии после запятой
// (1e0 - 0 знаков после запятой (для модели Изинга), 1e100 - 100 знаков после запятой)

#define DEBUG true

void rotate(int spin);          // Считает энергию системы
void complete();

void mc();
void single();
void normalize();


#include "common.c"

int main(int argc, char **argv)
{
    MPI_Init(&argc,&argv);      //инициализация mpi
    int rank, size;
    
    MPI_Comm_size(MPI_COMM_WORLD, &size); //получение числа процессов
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //текущий id процесса
	printf("size = %d, rank = %d\n", size, rank);
    
    
    int seed=0;                 // Random seed
    printf("# Please, input random number seed:  ");
    scanf("%u",&seed);
    seed += rank;
    char filename[100];
    printf("# Please, input target filename: ");
    scanf("%s",filename);

    if (!readCSV(filename)){
        printf("# Error! File '%s' is unavaliable!\n",filename);
        return 0;
    }

    
    char intervalsFile[50] = "csv_examples/intervals.csv";          ///new
    
    if (!readCSVintervals(intervalsFile)){                          ///new
        printf("# Error! File '%s' is unavaliable!\n", intervalsFile);
        return 0;
    }
    
    for(int i=0; i < intervalsNum; ++i){
        printf("%lf \n", intervals[i]);
    }
    
    for(int i=0; i < intervalsNum; ++i){
        printf("%lf \n", intervalsE[i]);
    }
    
    printf("\n");

#ifdef DEBUG
    printf("# spins:");
    for (unsigned i=0;i<n;i++){
        printf("%d,",spins[i]);
    }
    printf("\n");

    printf("# a_neighbours:");
    for (unsigned i=0;i<n;i++){
        printf("%d,",a_neighbours[i]);
    }
    printf("\n");

    printf("# sequencies:");
    for (unsigned i=0;i<n;i++){
        printf("%d,",sequencies[i]);
    }
    printf("\n");

    printf("# neighbours:");
    for (unsigned i=0;i<eCount;i++){
        printf("%d,",neighbours[i]);
    }
    printf("\n");

    printf("# energies:");
    for (unsigned i=0;i<eCount;i++){
        printf("%f,",energies[i]);
    }
    printf("\n");
#endif


    printf("\n# e = %lf, emin = %lf, emax = %lf\n",e,emin,emax);

    printf("\n# e = %lf\n",e);
    
    printf("# initial energy = %lf\n",e);
    
    factor = 0.8; // Критерий плоскости гистограммы H
    nfinal = 24; // число WL-циклов
    
    unsigned ie;
    for(ie=0; ie<=eCount; ie++){
      g[ie]=0;
      hist[ie]=0;
    }
    
    srand(seed);
    mc();
    normalize();

    for(ie=0; ie<=histSize; ie++){
      if (nonzero[ie] == 1) {
        printf("%e  %e  %e  %d\n",(double)(ie+emin)/PRECISION,g[ie],g[ie]/n,hist[ie]);
      }
    }
    
    complete();
    
    
    MPI_Finalize();
    return 0;
}

void rotate(int spin){
    double dE=0;
    spins[spin] *= -1;
    for(unsigned i = sequencies[spin]; i<sequencies[spin]+a_neighbours[spin]; ++i){
        dE += energies[i]*spins[neighbours[i]]*spins[spin]*2;
    }
    e += dE;
}

// clean arrays
void complete(){
    free(spins);
    free(a_neighbours);
    free(neighbours);
    free(sequencies);
    free(energies);

    free(g);
    free(visit);
    free(hist);
    free(nonzero);
    free(intervals);
    free(intervalsE);
}

/*
        monte carlo update
*/
void mc()
{
  unsigned ie,tt;
  int check,flag;
  int step, totalstep;
  int count;
  double sum;

/*   initialization  */
  totalstep=0;
  f=1;

  for(ie=0; ie<=histSize; ie++){
    nonzero[ie]=0;
  }

  for( tt = 0; tt <= nfinal; tt++){

    flag=0;
    step=0;

    for(ie=0; ie<=histSize; ie++){
      visit[ie]=0;
    }

    while(flag == 0){

      single();

      step++;

      if(step%1000==0){

        for(ie=0; ie<=histSize; ie++){
          if(visit[ie] > 0) {nonzero[ie]=1;}
        }

        count=0;
        sum=0;
        for(ie=0; ie<=histSize; ie++){
          if(nonzero[ie]==1) {
            count++;
            sum+=visit[ie];
          }
        }

        check=1; 
        for(ie=0; ie<=histSize; ie++){
          if(nonzero[ie]==1) {
            if(visit[ie] < factor*(sum/count)){check=0;}
          }
        }

        if (false && step%100000) //написать true для дебаг-вывода в файл каждые 100000 шагов
            dumpArrays();

        if(check==1){flag++;}
      }
    }

    gupdate();

    totalstep += step;

    printf("# n=%2d    MCS=%9d\n",tt,totalstep);

    f = f/2;
  }
  printf("# final   MCS=%9d\n",totalstep);

}


/*   single spin flip */
void single()
{
    unsigned la,la1;
    double energyOld;
    double ga,gb;

    int eoKey, enKey;
    
    for(la1=0; la1 <= n-1; la1++){
        la=rand()%n;
        energyOld = e;
        rotate(la);

        eoKey = (int)((energyOld-emin)*PRECISION);
        enKey = (int)((e-emin)*PRECISION);
        
        ga = g[eoKey];
        gb = g[enKey];
        
        if(exp(ga-gb) <= (double)rand()/RAND_MAX){
            spins[la] *= -1;
            e = energyOld;
            enKey = eoKey;
        }
        
        g[enKey]     += f;
        visit[enKey] += 1;
        hist[enKey]  += 1;
    }
}

void normalize()
{
    unsigned ie;
    double gmax, sum, a;
    
    gmax = -1000;
    for(ie=0; ie<histSize; ie++){
        if (nonzero[ie] == 1) {
            if(g[ie]>gmax){
                gmax = g[ie];
            }
        }
    }
    
    sum=0;
    for(ie=0; ie<histSize; ie++){
        if (nonzero[ie] == 1) {
            sum += exp(g[ie]-gmax);
        }
    }
    
    a = n*log(2) - gmax - log(sum);
    
    for(ie=0; ie<histSize; ie++){
        if (nonzero[ie] == 1) {
            g[ie] += a;
        }
    }
}
