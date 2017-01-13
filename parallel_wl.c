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
#include <stdbool.h>


int n; // количество спинов
signed char *spins; //массив направления спинов. По умолчанию +1. n - число считанных спинов (число незакомментированных строк в csv-файле).
unsigned short *a_neighbours; //число соседей каждого спина. Считается как число энергий в соответствующей строке в csv-файле.
unsigned short *neighbours; // соседи каждого спина
unsigned int *sequencies; //для каждого спина описывает, с какого ключа в массиве energies[] начинают описываться парные энергии
double *energies; //сами энергии из файла. Описывается как одномерный массив. Длина массива - число парных энергий в csv-файле.
float emin, emax; //минимумы и максимумы энергии
float e; //текущая энергия системы
unsigned eCount=0; //число пар энергий

#define PRECISION 3 //Сколько знаков учитывать в энергии после запятой


void readCSV(char* filename);
void rotate(int spin); // Считает энергию системы
void complete();


int main(void)
{
    readCSV("csv_examples/HC_30_lr.csv");

    printf("\n");
    printf("spins:");
    for (unsigned i=0;i<n;i++){
        printf("%d,",spins[i]);
    }
    printf("\n");

    printf("a_neighbours:");
    for (unsigned i=0;i<n;i++){
        printf("%d,",a_neighbours[i]);
    }
    printf("\n");

    printf("sequencies:");
    for (unsigned i=0;i<n;i++){
        printf("%d,",sequencies[i]);
    }
    printf("\n");

    printf("neighbours:");
    for (unsigned i=0;i<eCount;i++){
        printf("%d,",neighbours[i]);
    }
    printf("\n");

    printf("energies:");
    for (unsigned i=0;i<eCount;i++){
        printf("%f,",energies[i]);
    }
    printf("\n");

    
    
    printf("\ne = %lf, emin = %lf, emax = %lf\n",e,emin,emax);
    
    rotate(5);
    rotate(1);
    rotate(4);
    rotate(8);
    
    printf("\ne = %lf\n",e);
    
    complete();
}


void readCSV(char *filename){

    char c; //считанный из файла символ
    char symb[100]; //символ энергии в текстовом файле

    //get system sizes
    bool isFirstLine=true;
    n=0;
    FILE *file2 = fopen(filename, "r");
    int fpos = 1, lastFpos=0;
    do{
        c = fgetc(file2);
        while(c=='#'){
            do c = fgetc(file2); while (c != '\n');
            c = fgetc(file2);
        }
        if (isFirstLine && c==';')
            ++n;
        if (c=='\n')
            isFirstLine=false;

        if (c==';' || c=='\n') {
            if (fpos-1 != lastFpos)
                ++eCount;
            lastFpos = fpos;
        }

        fpos++;
    } while (c != EOF);
    ++n;
    fclose(file2);


    // read data







    FILE *file = fopen(filename, "r");

    bool firstSymbolInLine=true, skipFlag=false;
    double parsedNumber;
    int numInSymb=0;
    symb[0]='\0';
    int row=0; //line number in file (not account the commented lines)
    int col=0; //column number in line (taking to accound the ';' symbols)
    int neighCount=0; //
    int energyNum=0; //holds actual count of previously parsed energies
    e = 0;
    emax = 0;
    
    do {
        c = fgetc(file);

        if (firstSymbolInLine && c=='#'){ //if it is comment, skip the line
            skipFlag=true; //skip to end of line
        }
        firstSymbolInLine=false;

        if (!skipFlag){
            if (c==';' || c=='\n' || c == EOF){ //if we found a number, process it
                if (numInSymb!=0){
                    sscanf( symb, "%lf", &parsedNumber );
                    neighbours[energyNum] = col;
                    energies[energyNum] = parsedNumber;
                    e += parsedNumber;
                    emax += fabs(parsedNumber);


                    printf("%f\t",parsedNumber);
                    numInSymb=0;
                    ++neighCount;
                    ++energyNum;
                }
                ++col;
            } else {
                symb[numInSymb] = c;
                symb[numInSymb+1] = '\0';
                ++numInSymb;
            }

            if (c=='\n' || c == EOF){
                a_neighbours[row] = neighCount;
                sequencies[row] = energyNum-neighCount;
                col=0;
                neighCount=0;
                spins[row]=1;
                ++row;
                printf("\n");
            }
        }

        if (c=='\n'){ //if it is newline, mark the flag
            firstSymbolInLine=true;
            skipFlag=false;
        }
    } while (c != EOF);
    
    emax/=2;
    e/=2;
    emin = -emax;
    
    fclose(file);

}

void rotate(int spin){
    float dE=0;
    spins[spin] *= -1;
    for(int i = sequencies[spin]; i<sequencies[spin]+a_neighbours[spin]; ++i){
        dE += energies[i]*spins[neighbours[i]]*spins[spin]*2;
    }
    e += dE;
}

void complete(){
    // clean arrays
    free(spins);
    free(a_neighbours);
    free(neighbours);
    free(sequencies);
    free(energies);
}
