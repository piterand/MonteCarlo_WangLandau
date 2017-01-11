/*
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


int n; // количество спинов
signed char *spins; //массив направления спинов. По умолчанию +1. n - число считанных спинов (число незакомментированных строк в csv-файле).
unsigned short *a_neighbours; //число соседей каждого спина. Считается как число энергий в соответствующей строке в csv-файле.
unsigned short *neighbours; // соседи каждого спина
unsigned int *sequencies; //для каждого спина описывает, с какого ключа в массиве energies[] начинают описываться парные энергии
float *energies; //сами энергии из файла. Описывается как одномерный массив. Длина массива - число парных энергий в csv-файле.
float emin, emax; //минимумы и максимумы энергии
float e; //текущая энергия системы

#define PRECISION 3 //Сколько знаков учитывать в энергии после запятой


void readCSV(char* filename);
void rotate(int spin); // Считает энергию системы
void complete();


int main(void)
{


    fclose(f);

    printf("Testmass %s",testmas);
}


void readCSV(char *filename){
    // reserve memory for arrays
    spins=(signed char *) malloc(n*sizeof(signed char));
    a_neighbours=(unsigned short *) malloc(n*sizeof(unsigned short));
    neighbours=(unsigned short *) malloc(n*sizeof(unsigned short));     //поменять размер
    sequencies=(unsigned int *) malloc(n*sizeof(unsigned int));
    energies=(float *) malloc(n*sizeof(float ));                        //поменять размер


    // read data


    float *testmas;

    FILE *f = fopen("csv_examples/simplest_exmple.csv","r");
    fscanf(f,"%[,]%s",testmas);
    testmas=(float *) malloc(100*sizeof(float));

}

void rotate(int spin){

}

void complete(){
    // clean arrays
    free(spins);
    free(a_neighbours);
    free(neighbours);
    free(sequencies);
    free(energies);
}
