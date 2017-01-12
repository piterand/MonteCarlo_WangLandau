#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

int n = 16; // количество спинов
signed char *spins; //массив направления спинов. По умолчанию +1. n - число считанных спинов (число незакомментированных строк в csv-файле).
unsigned short *a_neighbours; //число соседей каждого спина. Считается как число энергий в соответствующей строке в csv-файле.
unsigned short *neighbours; // соседи каждого спина
unsigned int *sequencies; //для каждого спина описывает, с какого ключа в массиве energies[] начинают описываться парные энергии
float *energies; //сами энергии из файла. Описывается как одномерный массив. Длина массива - число парных энергий в csv-файле.
float emin, emax; //минимумы и максимумы энергии
float e; //текущая энергия системы

#define PRECISION 3 //Сколько знаков учитывать в энергии после запятой


void rotate(int spin); // Считает энергию системы


int main(void)
{
    spins=(signed char *) malloc(n*sizeof(signed char));
    /*
    spins[n]={0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,0,
           1,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,
           0,1,0,1,0,0,1,0,0,0,0,0,0,0,1,0,
           1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,
           1,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,
           0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,
           0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,
           0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,0,
           0,0,0,0,1,0,0,0,0,1,0,1,1,0,0,0,
           0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,0,
           0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,
           0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,
           1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,
           0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,
           0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
           0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0};
*/
    
    a_neighbours=(unsigned short *) malloc(n*sizeof(unsigned short));
    neighbours=(unsigned short *) malloc(n*n*sizeof(unsigned short));     //поменять размер
    sequencies=(unsigned int *) malloc(n*sizeof(unsigned int));
    energies=(float *) malloc(n*n*sizeof(float ));                        //поменять размер
    
    int neighbours_temp[64] = {1,3,4,12,//0
                               0,2,5,13,//1
                               1,3,6,14,//2
                               0,2,7,15,//3
                               0,5,7,8,//4
                               1,4,6,9,//5
                               2,5,7,10,//6
                               3,4,6,11,//7
                               4,9,11,12,//8
                               5,8,10,13,//9
                               6,9,11,14,//10
                               7,8,10,15,//11
                               0,8,13,15,//12
                               1,9,12,14,//13
                               2,10,13,15,//14
                               3,11,12,14//15
                              };
    
    for(int i = 0; i<n; ++i){
        spins[i] = 1;
        a_neighbours[i] = 4;
        sequencies[i] = 4*i;
    }
    
    emax = 0;
    e = 0;
    
    for(int i = 0; i<n*4; ++i){
        neighbours[i] = neighbours_temp[i];
        energies[i] = 1;
        emax += abs(energies[i]);
        e += energies[i];
    }
    
    emax/=2;
    e/=2;
    emin = -emax;
    
    printf("\ne = %lf, emin = %lf, emax = %lf\n",e,emin,emax);
    
    rotate(5);
//    rotate(5);
    rotate(1);
    rotate(4);
    rotate(8);
//    rotate(5);
//    rotate(1);
//    rotate(4);
//    rotate(8);
    
    printf("\ne = %lf\n",e);
    
    printf("\nGoodbye World!\n");
    
    free(spins);
    free(a_neighbours);
    free(neighbours);
    free(sequencies);
    free(energies);
    return 0;
}

void rotate(int spin){
    float dE=0;
    spins[spin] *= -1;
    for(int i = sequencies[spin]; i<sequencies[spin]+a_neighbours[spin]; ++i){
        dE += energies[i]*spins[neighbours[i]]*spins[spin]*2;
//        printf("\ndE = %lf\n",dE);
    }
    e += dE;
}
