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

    // reserve memory for arrays
    spins=(signed char *) malloc(n*sizeof(signed char));
    a_neighbours=(unsigned short *) malloc(n*sizeof(unsigned short));
    neighbours=(unsigned short *) malloc(n*sizeof(unsigned short));     //поменять размер
    sequencies=(unsigned int *) malloc(n*sizeof(unsigned int));
    energies=(float *) malloc(n*sizeof(float ));                        //поменять размер

       FILE *file = fopen("csv_examples/simplest_exmple.csv", "r");

       const char sym = 59;
       char cdum[1000];
       char cn;
       int count=0;


       while(cn = fgetc(file)=='#')
       {
            fscanf(file,"%[^\n]%*c",cdum);
       }
       fseek(file,-2,SEEK_CUR);
       cn = fgetc(file);
       printf("%c",cn);

       int coursor=ftell(file);




       while((cn = fgetc(file)) != '\n')
       {
           if (cn == sym)
                  count++;
       }

       printf("%u",count);

      // fclose(file);
/////////////////////////////////////////////////////////

       n=(count+1)^2;
       for(int i=0;i<n;i++)
           *(spins+i)=1;
       //file = fopen("csv_examples/simplest_exmple.csv", "r");
       fseek(file,coursor,SEEK_SET);
// записываем энергии

       for(int i=0;i<20;i++)
       {
       fscanf(file,"%d",energies);
         printf("\n%d",energies[i]);
       }

//       while(! feof(file))
//       {
//            fscanf(file,"%[^\n]%*c%[^\n]%*c",cdum);
//       }
//       for(int i=0;i<10;i++)
//           printf("\n%d\n",cdum[i]);
//        cn = fgetc(file);
//       printf("\n%c\n",cn);



       fclose(file);


       //       while(! feof(file))
       //       {
       //                fscanf(file, "%*[^\n]%*c");
       //                cnt++;
       //          }

//       while( fscanf(file, "%79[^\n]\n", array) == 1 )
//       {
//           printf("%s\n", array);
//       }

//    float *testmas;
//    testmas=(float *) malloc(100*sizeof(float));

//    FILE *f = fopen("simplest_exmple.csv","r");

//    fscanf(f,"%s%d",testmas);


//    testmas[0]=2.78;
//    testmas[1]=2;
//    testmas[2]=3;

//    printf ("\n Массив %c [%i] \n",testmas,3);
//    for (int i=0;i<4;i++)
//         printf ("%4.2f\n",*(testmas+i));

    //for(int i=0;i<3;i++)

        //printf("%d\n",testmas[i]);
     //fclose(f);
     return 0;



}


void readCSV(char *filename)
{


    // read data






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
