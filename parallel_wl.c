/**
        parallel_wl.c

        Parallel Wang-Landau method for different magnetic systems

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

// для параллельного кода
double *intervals;              //массив интервалов
double *intervalsE;             //массив интервалов значений
int intervalsNum=0;             //число значений интервалов из файла
double emin_for_current_rank, emax_for_current_rank; //минимумы и максимумы энергии для конкретного процесса
unsigned int rank, size;

// для функции exchanage
char exchange_buffer[10000];    //буфер обмена !! Должен быть больше количества спинов в системе!!
double exchange_energy;         //обмениваемое значение энергии
double exchange_Ge_a;           //обмениваемое значение G(энергии)
double exchange_Ge_b;           //обмениваемое значение G(энергии)
signed char *exchange_spins;
bool exchange_status;


double *g;
unsigned *visit;
unsigned *hist;
int *nonzero;

double f;                       // Модификационный фактор (уменьшается с каждым WL-шагом)
double factor = 0.8;            // Критерий плоскости гистограммы H
unsigned nfinal = 24;           // число WL-циклов

#define PRECISION 1e1           // Точность 1eX, где X - Сколько знаков учитывать в энергии после запятой
                                // (1e0 - 0 знаков после запятой (для модели Изинга), 1e100 - 100 знаков после запятой)

//#define DEBUG true            // Что бы отключить(включить) режим дебага нужно закомментировать (раскомментировать) эту строку.

int readCSVintervals(char *filename); //считывает интервалы из файла
void rotate(int spin);          // Считает энергию системы
void complete();

void mc(double eFrom, double eTo);
void single(double eFrom, double eTo);
bool exchange(unsigned a, unsigned b);
void normalize();

#include "common.c"

int main(int argc, char **argv)
{
    MPI_Init(&argc,&argv);      //инициализация mpi

    MPI_Comm_size(MPI_COMM_WORLD, &size); //получение числа процессов
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //текущий id процесса
    MPI_Status status;

    printf("#size = %d, rank = %d\n", size, rank);

    int seed=0;                 // Random seed
    char filename[300];         // целевой файл с энергиями
    char filenameinterval[300]; // целевой файл с интервалами

    if(rank==0)
    {
        printf("# Please, input random number seed:  ");
        scanf("%u",&seed);

        //printf("# Please, input target energy filename: ");
        //scanf("%s",filename);

        //printf("# Please, input target intervals filename: ");
        //scanf("%s",filenameinterval);

        strcpy(filename,"/home/petr/scienceworks/Programs_with_Git/wanglandauparallel/csv_examples/square_ising_4x4.csv");
        strcpy(filenameinterval, "/home/petr/scienceworks/Programs_with_Git/wanglandauparallel/csv_examples/intervals.csv");


    }

    MPI_Bcast(&seed,1, MPI_INT, 0, MPI_COMM_WORLD);                 // рассылаем seed
    MPI_Bcast(filename,300, MPI_CHAR, 0, MPI_COMM_WORLD);           // рассылаем имя файла с энергиями
    MPI_Bcast(filenameinterval,300, MPI_CHAR, 0, MPI_COMM_WORLD);   // рассылаем имя файла с интервалами

    seed += rank;

    printf("#rank = %d, seed = %d\n", rank,seed);

    if (!readCSV(filename)){
        printf("# Error!! File '%s' is unavaliable!\n",filename);
        return 0;
    }


    if (!readCSVintervals(filenameinterval)){                          ///new
        printf("# Error! File '%s' is unavaliable!\n", filenameinterval);
        return 0;
    }

    // тут нужно написать распределение интервалов по процессам, пока делаем вручную
    //printf("\n!!!intervalsNum=%d\n",intervalsNum);

    if(size>intervalsNum/2){
        if(rank=0){
            printf("\n!!!Error, please enater corresponding number of process equal number of intervals = %d",intervalsNum);
        }
        return 0;
    }

    exchange_spins=(signed char *) malloc(n*sizeof(signed char));   //массив спинов для обмена

    if(rank<intervalsNum/2){
        emin_for_current_rank=intervalsE[2*rank];
        emax_for_current_rank=intervalsE[2*rank+1];
        printf("\n!!!my rank=%d, emin=%f    emax=%f\n",rank,emin_for_current_rank, emax_for_current_rank);
    }

#ifdef DEBUG

    for(int i=0; i < intervalsNum; ++i){
        printf("%lf \n", intervals[i]);
    }

    for(int i=0; i < intervalsNum; ++i){
        printf("%lf \n", intervalsE[i]);
    }

    printf("\n");


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


 //   printf("\n# e = %lf, emin = %lf, emax = %lf\n",e,emin,emax);

//    printf("# initial energy = %lf\n",e);


    unsigned ie;
    for(ie=0; ie<histSize; ++ie){
      g[ie]=0;
      hist[ie]=0;
    }

    srand(seed);
    mc(emin,emax);
    normalize();
    exchange(0,1);
    MPI_Barrier(MPI_COMM_WORLD);

    // вывод
//    printf("# e  g[ie]  g[ie]/n  hist[ie]\n");
//    for(ie=0; ie<histSize; ie++){
//      if (nonzero[ie] == 1) {
//        printf("%e  %e  %e  %d\n",(double)ie/PRECISION+emin,g[ie],g[ie]/n,hist[ie]);
//      }
//    }

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
void mc(double eFrom, double eTo)
{
  unsigned ie,tt;
  int check,flag;
  long long step, totalstep; step, totalstep;
  int count;
  double sum;

/*   initialization  */
  totalstep=0;
  f=1;

  for(ie=0; ie<histSize; ie++){
    nonzero[ie]=0;
  }

  for( tt = 0; tt <= nfinal; tt++){

    flag=0;
    step=0;

    for(ie=0; ie<histSize; ie++){
      visit[ie]=0;
    }

    while(flag == 0){

      single(eFrom,eTo);


      step++;

      if(step%1000==0){

        for(ie=0; ie<histSize; ie++){
          if(visit[ie] > 0) {nonzero[ie]=1;}
        }

        count=0;
        sum=0;
        for(ie=0; ie<histSize; ie++){
          if(nonzero[ie]==1) {
            count++;
            sum+=visit[ie];
          }
        }

        check=1;
        for(ie=0; ie<histSize; ie++){
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
void single(double eFrom, double eTo){
    unsigned la,la1;        // итераторы
    double energyOld;       // старая энергия
    double ga,gb;           // g[старой энергии] и g[новой энергии]

    int eoKey, enKey;       // номер столбика гистограммы энергий старой и новой


    for(la1=0; la1 <= n-1; la1++){  //цикл выполняется n раз, не знаю почему
        la=rand()%n;            // выбираем случайный спин
        energyOld = e;          // записываем старую энергию
        rotate(la);             // переворачиваем выбранный спин

        eoKey = (int)((energyOld-emin)*PRECISION); //вычисляем номер столбика гистограммы для старой энергии
        enKey = (int)((e-emin)*PRECISION);         //вычисляем номер столбика гистограммы для новой энергии


        ga = g[eoKey];          // g[старой энергии]
        gb = g[enKey];          // g[новой энергии]

        if(exp(ga-gb) <= (double)rand()/RAND_MAX){      // условия переворота, если не принимаем, то заходим внутрь цикла
            spins[la] *= -1;        // не принимаем новую конфигурацию, обратно переворачиваем спин
            e = energyOld;          // обратно записываем старую энергию
            enKey = eoKey;          // берем старый столбик гистограммы
        }

        if(e < eFrom || e > eTo){      // условия переворота, если не принимаем, то заходим внутрь цикла
            spins[la] *= -1;        // не принимаем новую конфигурацию, обратно переворачиваем спин
            e = energyOld;          // обратно записываем старую энергию
            enKey = eoKey;          // берем старый столбик гистограммы
        }

        g[enKey]     += f;          // прибавляем f в текущий столбик гистограммы (так как тут хрянятся логарифмы)
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

int readCSVintervals(char *filename){
    int numerOfStrings = 0;
    char c;                         //считаный из файла символ
    char symb[100000];                 //символ энергии в текстовом файле

    //get system sizes
    bool isFirstLine=true;
    FILE *file = fopen(filename, "r");

    if (!file)
        return 0;

    for(c=fgetc(file);c=='\n'||c=='\r'||c=='#';c=fgetc(file)){ //пропуск комментариев
        fgets(symb,100000,file);
    }

    fseek(file,-1,SEEK_CUR);       // сдвиг курсора на один символ назад
    int coursor=ftell(file);       // положение курсора начала данных

    do{
        c = fgetc(file);
        if (c=='\n' && isFirstLine){
            isFirstLine=false;
        }
        if ((c=='\n' && !isFirstLine) || c == EOF){
            isFirstLine=false;
            numerOfStrings++;
        }
    } while (c != EOF);

    intervals = (double *) calloc(numerOfStrings*2,sizeof(double));
    intervalsE = (double *) calloc(numerOfStrings*2,sizeof(double));
    unsigned i;

    // read data

    fseek(file,coursor,SEEK_SET);      //устанавливаем курсор в начало данных

    double parsedNumber;
    int numInSymb=0;
    symb[0]='\0';
    intervalsNum = 0;

    do {
        c = fgetc(file);

        if (c==';' || c=='\n' || c == EOF){ //if we found a number, process it
            if (numInSymb!=0){
                sscanf(symb, "%lf", &parsedNumber);
                intervals[intervalsNum] = parsedNumber;
                intervalsE[intervalsNum] = (emax-emin) * parsedNumber + emin;

                numInSymb=0;
                ++intervalsNum;
            }
        } else {
            symb[numInSymb] = c;
            symb[numInSymb+1] = '\0';
            ++numInSymb;
        }
    } while (c != EOF);

    fclose(file);
    int tr=0;
    //debug
//    for(tr=0;tr<numerOfStrings*2;++tr)
//    {
//        printf("\nintervals[%d]=%f",tr,intervals[tr]);
//    }
//    for(tr=0;tr<numerOfStrings*2;++tr)
//    {
//        printf("\nintervalsE[%d]=%f",tr,intervalsE[tr]);
//    }

    return 1;
}



bool exchange(unsigned a, unsigned b){

    int position = 0;
    int current_energy = 0;
    int current_energy2 = 0;
    MPI_Status status;
    exchange_status = 1;
    double exchange_probobility;
    double exchange_probobility_final;
    double exchange_rand;
    double exchange_emin;
    double exchange_emax;

    if (rank == a)
     {
         current_energy = (int)((e-emin)*PRECISION);
         printf("#1Smy rank=%d   I send my e_a=%f    g_a[e_a]=%f\n", rank,e, g[current_energy]);//debug

         MPI_Pack(&e, 1,MPI_DOUBLE, exchange_buffer, 100, &position, MPI_COMM_WORLD);
         MPI_Pack(&g[current_energy], 1,  MPI_DOUBLE, exchange_buffer, 100, &position, MPI_COMM_WORLD);
         MPI_Pack(&emin_for_current_rank, 1,MPI_DOUBLE, exchange_buffer, 100, &position, MPI_COMM_WORLD);
         MPI_Pack(&emax_for_current_rank, 1,MPI_DOUBLE, exchange_buffer, 100, &position, MPI_COMM_WORLD);
         MPI_Send(exchange_buffer, position, MPI_PACKED, b, 100, MPI_COMM_WORLD);   // упаковкаи отсылка e_a и g_a(e_a)

         //1st check to exit
         MPI_Bcast(&exchange_status,1, MPI_BYTE, b, MPI_COMM_WORLD);    // recive signal, if 0 -> exit.
         if(exchange_status==0)
             return 0;
         else{

             MPI_Recv(exchange_buffer, 10000, MPI_PACKED, b, 10000, MPI_COMM_WORLD, &status); // получаем от b e_b, g_b(e_b), g_b(e_a), и массив spins
             MPI_Unpack(exchange_buffer, 10000, &position, &exchange_energy, 1, MPI_DOUBLE, MPI_COMM_WORLD);
             MPI_Unpack(exchange_buffer, 10000, &position, &exchange_Ge_b, 1,  MPI_DOUBLE, MPI_COMM_WORLD);
             MPI_Unpack(exchange_buffer, 10000, &position, &exchange_Ge_a, 1,  MPI_DOUBLE, MPI_COMM_WORLD);
             MPI_Unpack(exchange_buffer, 10000, &position, exchange_spins, n,  MPI_SIGNED_CHAR, MPI_COMM_WORLD);
             printf("#2Rmy rank=%d   exch_E=%f   g_b[e_b]=%f g_b[e_a]=%f\n", rank,exchange_energy, exchange_Ge_b,exchange_Ge_a);//debug

             current_energy2=(int)((exchange_energy-emin)*PRECISION);
             exchange_probobility =(g[current_energy]*exchange_Ge_b)/(exchange_Ge_a*g[current_energy2]);
             if(exchange_probobility<1)
                 exchange_probobility_final=exchange_probobility;
             else
                 exchange_probobility_final=1;
             printf("#3My rank=%d   exchange_probobility =(g[current_energy] =  %f   * exchange_Ge_b=%f)  / exchange_Ge_a=%f * g[current_energy2]=%f = %f\n", rank,g[current_energy], exchange_Ge_b,exchange_Ge_a,g[current_energy2],exchange_probobility);//debug
             exchange_rand=(double)rand()/RAND_MAX;

             //2st check to exit
             if(exchange_probobility_final<exchange_rand){  // fail
                 exchange_status=0;
                 MPI_Bcast(&exchange_status,1, MPI_BYTE, a, MPI_COMM_WORLD);    // recive signal, 0 -> exit.
                 printf("2nd Cancel,exchange_probobility_final = %f    exchange_rand = %f",exchange_probobility_final,exchange_rand);//debug
                 return 0;
             }
             else{  //continue
                 printf("Succses,exchange_probobility_final = %f    exchange_rand = %f",exchange_probobility_final,exchange_rand);//debug
                 MPI_Bcast(&exchange_status,1, MPI_BYTE, a, MPI_COMM_WORLD);    // send signal 1

                 MPI_Send(spins, n, MPI_SIGNED_CHAR, b, 1010, MPI_COMM_WORLD);  // send spins_a

                 for (unsigned i=0;i<n;i++){
                     spins[i]=exchange_spins[i];    // change spins_b -> spins_a
                 }
                 e = exchange_energy;

                 g[current_energy2]     += f;
                 visit[current_energy2] += 1;
                 hist[current_energy2]  += 1;



             }



         }

     }

     if (rank == b)
     {
         current_energy = (int)((e-emin)*PRECISION);
         MPI_Recv(exchange_buffer, 100, MPI_PACKED, a, 100, MPI_COMM_WORLD, &status); // получаем e_a, g_a(e_a)
         MPI_Unpack(exchange_buffer, 100, &position, &exchange_energy, 1, MPI_DOUBLE, MPI_COMM_WORLD);
         MPI_Unpack(exchange_buffer, 100, &position, &exchange_Ge_a, 1,  MPI_DOUBLE, MPI_COMM_WORLD);
         MPI_Unpack(exchange_buffer, 100, &position, &exchange_emin, 1, MPI_DOUBLE, MPI_COMM_WORLD);
         MPI_Unpack(exchange_buffer, 100, &position, &exchange_emax, 1,  MPI_DOUBLE, MPI_COMM_WORLD);



         printf("#1Rmy rank=%d   My e_b=%f g_b[e_b]=%f, and I recive e_a=%f    g_a[e_a]=%f\n", rank,e,g[current_energy],exchange_energy,exchange_Ge_a);//debug

         //1st проверка на поподание энергий в окно
         if(exchange_energy>emax_for_current_rank||exchange_energy<emin_for_current_rank||e<exchange_emin||e>exchange_emax){ // если не попала в интервал энергии текущего процесса, то выйти из ф-ции
             exchange_status=0;
             MPI_Bcast(&exchange_status,1, MPI_BYTE, b, MPI_COMM_WORLD);    // рассылка сигнала exchange_status =0 выхода из фунуции.
             printf("\n1st Cancel\n");
             return 0;
         }

         else{
             MPI_Bcast(&exchange_status,1, MPI_BYTE, b, MPI_COMM_WORLD);    // рассылка exchange_status=1 все ок, продолжаем

             current_energy2 = (int)((exchange_energy-emin)*PRECISION); // g_b(e_a)

             printf("#2Smy rank=%d I send my exch_E=%f   g_b[e_b]=%f g_b[e_a]=%f\n", rank,e, g[current_energy],g[current_energy2]);//debug
             MPI_Pack(&e, 1,MPI_DOUBLE, exchange_buffer, 10000, &position, MPI_COMM_WORLD);
             MPI_Pack(&g[current_energy], 1,  MPI_DOUBLE, exchange_buffer, 10000, &position, MPI_COMM_WORLD);
             MPI_Pack(&g[current_energy2], 1,  MPI_DOUBLE, exchange_buffer, 10000, &position, MPI_COMM_WORLD);
             MPI_Pack(spins, n,  MPI_SIGNED_CHAR, exchange_buffer, 10000, &position, MPI_COMM_WORLD);
             MPI_Send(exchange_buffer, position, MPI_PACKED, a, 10000, MPI_COMM_WORLD); // отсылаем e_b, g_b(e_b), g_b(e_a), и массив spins

             //2nd check
             MPI_Bcast(&exchange_status,1, MPI_BYTE, a, MPI_COMM_WORLD);    // recive signal,
             if(exchange_status==0)
                 return 0;  // break
             else{
                 // continue
                 MPI_Recv(spins, n, MPI_SIGNED_CHAR, a, 1010, MPI_COMM_WORLD, &status);   // recive and change spins_a -> spins_b

                 e = exchange_energy;   // change energy

                 g[current_energy2]     += f;
                 visit[current_energy2] += 1;
                 hist[current_energy2]  += 1;

             }


         }
     }


     if (rank!=a && rank!=b){
         //1st check
         MPI_Bcast(&exchange_status,1, MPI_BYTE, b, MPI_COMM_WORLD);
         if(exchange_status==0)
             return 0;
         //2nd check
         MPI_Bcast(&exchange_status,1, MPI_BYTE, a, MPI_COMM_WORLD);    // recive signal,
         if(exchange_status==0)
             return 0;

     }



    return true;
}

