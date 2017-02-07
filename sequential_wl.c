/**
        sequential_wl.c


        Sequentional Wang-Landau method for different magnetic systems

           programmed by:
            Makarov Aleksandr
            Andriushchenko Petr
            Shevchenko Yuriy


*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>



unsigned n;                     //количество спинов
signed char *spins;             //массив направления спинов. По умолчанию +1. n - число считанных спинов (число незакомментированных строк в csv-файле).
unsigned short *a_neighbours;   //число соседей каждого спина. Считается как число энергий в соответствующей строке в csv-файле.
unsigned short *neighbours;     //соседи каждого спина
unsigned int *sequencies;       //для каждого спина описывает, с какого ключа в массиве energies[] начинают описываться парные энергии
double *energies;               //сами энергии из файла. Описывается как одномерный массив. Длина массива - число парных энергий в csv-файле.
double emin, emax;              //минимумы и максимумы энергии
double e;                       //текущая энергия системы
unsigned eCount=0;              //число пар энергий
unsigned histSize=0;            //число элементов в гистограммах

double *g;                      // гистограмма g[E], только хранятся !логарифмы!.
unsigned *visit;                //
unsigned *hist;                 // Гистограмма H
int *nonzero;                   // Массив, была ли хоть раз посещена данная энергия 1 или 0

double f;                       // Модификационный фактор (уменьшается с каждым WL-шагом)
double factor = 0.8;            // Критерий плоскости гистограммы H
unsigned nfinal = 24;           // число WL-циклов


#define PRECISION 1e2           // Точность 1eX, где X - Сколько знаков учитывать в энергии после запятой
                                // (1e0 - 0 знаков после запятой (для модели Изинга), 1e100 - 100 знаков после запятой)

//#define DEBUG true             // Режим дебага для отлавливания утечек памяти

void rotate(int spin);          // Считает энергию системы
void complete();

#include "common.c"             //подключить файл с общими функциями для параллельного и последовательного варианта WL

void mc();
void single();
void normalize();



int main(void)
{
    unsigned i;
    unsigned long seed=0;       // Random seed
    printf("# Please, input random number seed from 1 to 4 294 967 295:  ");
    if (scanf("%lu",&seed) == 1){}
    else{
        printf("# Error! Failed to read integer seed!\n");
        return 0;
    }


    char filename[100];
    printf("# Please, input target filename: ");
    if( scanf("%s",filename)==1){}
    else{
        printf("# Error! Failed to read file name!\n");
        return 0;
    }

    if (!readCSV(filename)){
        printf("# Error! File '%s' is unavaliable!\n",filename);
        return 0;
    }

#ifdef DEBUG
    printf("# spins:");
    for (i=0;i<n;i++){
        if (i>=n || i<0) printf("Error with memory working1");
        printf("%d,",spins[i]);
    }
    printf("\n");

    printf("# a_neighbours:");
    for (i=0;i<n;i++){
        if (i>=n || i<0) printf("Error with memory working2");
        printf("%d,",a_neighbours[i]);
    }
    printf("\n");

    printf("# sequencies:");
    for (i=0;i<n;i++){
        if (i>=n || i<0) printf("Error with memory working3");
        printf("%d,",sequencies[i]);
    }
    printf("\n");

    printf("# neighbours:");
    for (i=0;i<eCount;i++){
        if (i>=eCount || i<0) printf("Error with memory working4");
        printf("%d,",neighbours[i]);
    }
    printf("\n");

    printf("# energies:");
    for (i=0;i<eCount;i++){
        if (i>=eCount || i<0) printf("Error with memory working5");
        printf("%f,",energies[i]);
    }
    printf("\n");
#endif

    printf("\n#  initial energy = %lf, emin = %lf, emax = %lf\n",e,emin,emax);
    
    unsigned ie;
    for(ie=0; ie<histSize; ie++){
#ifdef DEBUG
        if (ie>=histSize || ie<0) printf("Error with memory working6");
#endif
        g[ie]=0;
        hist[ie]=0;
    }
    
    srand(seed);

    mc();
    normalize();
    
    printf("# e  g[ie]  g[ie]/n  hist[ie]\n");
    for(ie=0; ie<histSize; ie++){
        if (nonzero[ie] == 1) {
#ifdef DEBUG
            if (ie>=histSize || ie<0)
                printf("Error with memory working7");
#endif
            printf("%e  %e  %e  %d\n",(double)ie/PRECISION+emin,g[ie],g[ie]/n,hist[ie]);
        }
    }
    
    complete(); //очистка памяти
}

/// Переворот спина, подсчет изменения энергии
void rotate(int spin)
{
    unsigned i;
    double dE=0;
    spins[spin] *= -1;
    for(i = sequencies[spin]; i<sequencies[spin]+a_neighbours[spin]; ++i){
#ifdef DEBUG
        if (spin>=n || spin<0) printf("Error with memory working10");
        if (i>=eCount || i<0) printf("Error with memory working11");
        if (neighbours[i]>=eCount || neighbours[i]<0) printf("Error with memory working111");
#endif
        dE += energies[i]*spins[neighbours[i]]*spins[spin]*2;
    }
    e += dE;
}


/// Очистка памяти
void complete()
{
    // clean arrays
    free(spins);
    free(a_neighbours);
    free(neighbours);
    free(sequencies);
    free(energies);

    free(g);
    free(visit);
    free(hist);
    free(nonzero);
}


/// Монтекарло шаг
void mc(){
/*
        monte carlo update
*/
    unsigned ie,tt; //итераторы
    int check,flag;
    long long step, totalstep;
    int count;
    double sum;

    /*   initialization  */
    totalstep=0;
    f=1;

    for(ie=0; ie<histSize; ie++){       //обнуляем массив nonzero
#ifdef DEBUG
        if (ie>=histSize || ie<0) printf("Error with memory working12");
#endif
        nonzero[ie]=0;
    }


    for( tt = 0; tt <= nfinal; tt++){    // WL цикл

        flag=0;
        step=0;

        for(ie=0; ie<histSize; ie++){
#ifdef DEBUG
            if (ie>=histSize || ie<0) printf("Error with memory working5");
#endif
            visit[ie]=0;                  // обнуляем массив visit
        }

        while(flag == 0){

            single();

            step++;

            if(step%1000==0){             // каждые 1000 шагов

                for(ie=0; ie<histSize; ie++){
#ifdef DEBUG
                    if (ie>=histSize || ie<0) printf("Error with memory working5");
#endif
                    if(visit[ie] > 0) {nonzero[ie]=1;}        // проверяем, появились ли новые энергии
                }

                count=0;
                sum=0;
                for(ie=0; ie<histSize; ie++){
#ifdef DEBUG
                    if (ie>=histSize || ie<0) printf("Error with memory working5");
#endif
                    if(nonzero[ie]==1) {
                        count++;                        // подсчитываем количество встреченных энергий
                        sum+=visit[ie];                 // сумма посещений всех энергий
                    }
                }

                check=1;
                for(ie=0; ie<histSize; ie++){           // проверка на плоскоту
#ifdef DEBUG
                    if (ie>=histSize || ie<0) printf("Error with memory working5");
#endif
                    if(nonzero[ie]==1) {
                        if(visit[ie] < factor*(sum/count)){check=0;}    // sum/count = среднее значение количества посещений энергий
                    }                                                   // если количество посещений хоть одной энергии меньше среднего значения * factor, то проверка провалилась
                }

                if (false && step%100000) //написать true для дебаг-вывода в файл каждые 100000 шагов
                    dumpArrays();

                if(check==1){flag++;}
            }
        }

        gupdate();

        totalstep += step;

        printf("# n=%2d    MCS=%9d\n",tt,totalstep);    // ! на самом деле тут totalstep*n MCS, так как в функции single цикл по n

        f = f/2;
    }
    printf("# final   MCS=%9d\n",totalstep);

}

void single(){
// n spins flips.

    unsigned la,la1;                // итераторы
    double energyOld;               // старая энергия
    double ga,gb;                   // g[старой энергии] и g[новой энергии]

    int eoKey, enKey;               // номер столбика гистограммы энергий старой и новой
    
    for(la1=0; la1 <= n-1; la1++){  //цикл выполняется n раз, не знаю почему
        la=rand()%n;                // выбираем случайный спин
        energyOld = e;              // записываем старую энергию
        rotate(la);                 // переворачиваем выбранный спин

        eoKey = (int)((energyOld-emin)*PRECISION); //вычисляем номер столбика гистограммы для старой энергии
        enKey = (int)((e-emin)*PRECISION);         //вычисляем номер столбика гистограммы для новой энергии
#ifdef DEBUG
        if (eoKey>=histSize || eoKey<0) printf("Error with memory working12");
        if (enKey>=histSize || enKey<0) printf("Error with memory working12");
        if (la>=n || la<0) printf("Error with memory working12");
#endif

        ga = g[eoKey];              // g[старой энергии]
        gb = g[enKey];              // g[новой энергии]
        
        if(exp(ga-gb) <= (double)rand()/RAND_MAX){      // условия переворота, если не принимаем, то заходим внутрь цикла
            spins[la] *= -1;        // не принимаем новую конфигурацию, обратно переворачиваем спин
            e = energyOld;          // обратно записываем старую энергию
            enKey = eoKey;          // берем старый столбик гистограммы
        }

        g[enKey]     += f;          // прибавляем f в текущий столбик гистограммы (так как тут хрянятся логарифмы)
        visit[enKey] += 1;
        hist[enKey]  += 1;
    }
}

// нормализация гистограммы
void normalize()
{
    unsigned ie;
    double gmax, sum, a;
    
    gmax = -1000;
    for(ie=0; ie<histSize; ie++){
#ifdef DEBUG
        if (ie>=histSize || ie<0) printf("Error with memory working5");
#endif
        if (nonzero[ie] == 1) {
            if(g[ie]>gmax){
                gmax = g[ie];
            }
        }
    }
    
    sum=0;
    for(ie=0; ie<histSize; ie++){
#ifdef DEBUG
        if (ie>=histSize || ie<0) printf("Error with memory working5");
#endif
        if (nonzero[ie] == 1) {
            sum += exp(g[ie]-gmax);
        }
    }
    
    a = n*log(2) - gmax - log(sum);
    
    for(ie=0; ie<histSize; ie++){
#ifdef DEBUG
        if (ie>=histSize || ie<0) printf("Error with memory working5");
#endif
        if (nonzero[ie] == 1) {
            g[ie] += a;
        }
    }
}
