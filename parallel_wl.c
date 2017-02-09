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
#include <string.h>
#include <mpi.h>


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
double *intervalsE;             //массив интервалов значений энергий
unsigned intervalsNum=0;             //число значений интервалов из файла
double emin_for_current_rank, emax_for_current_rank; //минимумы и максимумы энергии для конкретного процесса
int rank, size;
int wl_step_count;              //текущий WL шаг
int wl_check_end;               //текущий WL шаг


// для функции exchanage
char exchange_buffer[10000];    //буфер обмена !! Должен быть больше количества спинов в системе!!
double exchange_energy;         //обмениваемое значение энергии
double exchange_Ge_a;           //обмениваемое значение G(энергии)
double exchange_Ge_b;           //обмениваемое значение G(энергии)
signed char *exchange_spins;
bool exchange_status;
int tdist;                      // по скольку процессов на каждый интервал


double *g;
unsigned *visit;
unsigned *hist;
int *nonzero;

double f;                       // Модификационный фактор (уменьшается с каждым WL-шагом)
double factor = 0.8;            // Критерий плоскости гистограммы H
unsigned nfinal = 24;           // число WL-циклов

int PRECISION; //!!!теперь задается пользователем при запуске программы
                                // Точность 1eX, где X - Сколько знаков учитывать в энергии после запятой
                                // (1e0 - 0 знаков после запятой (для модели Изинга), 1e100 - 100 знаков после запятой)

//#define DEBUG true            // Что бы отключить(включить) режим дебага нужно закомментировать (раскомментировать) эту строку.

int readCSVintervals(char *filename); //считывает интервалы из файла
void rotate(int spin);          // Считает энергию системы
void complete();
void showResult();              // Выводит гистограммы на экран

void mc(double eFrom, double eTo);      // фнукция запуска WL Монте-Карло в заданном интервале
void single(double eFrom, double eTo);  // функция переворота спина и попытка принятия новой системы в заданном интервале
bool exchange(int a, int b);  // функция обмена конфигурациями между потоками
void normalize();                       // нормализация гистограммы в конце рассчета

#include "common.c"

int main(int argc, char **argv)
{
    MPI_Init(&argc,&argv);      //инициализация mpi

    MPI_Comm_size(MPI_COMM_WORLD, &size); //получение числа процессов
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //текущий id процесса

    printf("#size = %d, rank = %d\n", size, rank);

    unsigned long seed=0;       // Random seed
    int prec=0;                 // Точность 1eX, где X - Сколько знаков учитывать в энергии после запятой
    char filename[300];         // целевой файл с энергиями
    char filenameinterval[300]; // целевой файл с интервалами

    if(rank==0)
    {
        printf("# Please, input random number seed from 1 to 4 294 967 295:  ");
        seed=1000;
        if (scanf("%lu",&seed) == 1){}
        else{
            printf("# Error! Failed to read integer seed!\n");
            return 0;
        }

        printf("# Please, chose precision X from 0 to 5(for example), where X - amount of numbers after dot. If you add 1, precision will iincrease 10 times: ");
        if (scanf("%u",&prec) == 1){}
        else{
            printf("# Error! Failed to read integer precision!\n");
            return 0;
        }
        prec = 0;
        PRECISION = pow(10,prec);   // !!Задание точности
        printf("# Precision = %d\n",PRECISION);


        printf("# Please, input target energy filename: ");
        scanf("%s",filename);

        printf("# Please, input target intervals filename: ");
        scanf("%s",filenameinterval);

    }
    MPI_Bcast(&seed,1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);                 // рассылаем seed
    MPI_Bcast(&PRECISION,1, MPI_INT, 0, MPI_COMM_WORLD);                 // рассылаем seed
    MPI_Bcast(filename,300, MPI_CHAR, 0, MPI_COMM_WORLD);           // рассылаем имя файла с энергиями
    MPI_Bcast(filenameinterval,300, MPI_CHAR, 0, MPI_COMM_WORLD);   // рассылаем имя файла с интервалами

    seed += rank;

    printf("#rank = %d, seed = %lu\n", rank,seed);

    if (!readCSV(filename)){
        printf("# Error!! File '%s' is unavaliable!\n",filename);
        return 0;
    }


    if (!readCSVintervals(filenameinterval)){                          ///new
        printf("# Error! File '%s' is unavaliable!\n", filenameinterval);
        return 0;
    }

    exchange_spins=(signed char *) malloc(n*sizeof(signed char));   //массив спинов для обмена

    //////////////// распределение интервалов по процессам, пока делаем вручную
    printf("\n!!!intervalsNum=%d\n",intervalsNum);

    if((size<(int)(intervalsNum/2)) || (size % (intervalsNum/2))!=0){

            printf("\n!!!Error, please enter number of process larger then number of intervals >= %d",intervalsNum);

        return 0;
    }

    //int tdrop = (int)size/(intervalsNum/2);
    //printf("tdrop=%d",tdrop);


//    if(rank>=intervalsNum/2){
//        emin_for_current_rank=intervalsE[2*rank];
//        emax_for_current_rank=intervalsE[2*rank+1];
//        printf("\n!!!my rank=%d, emin=%f    emax=%f\n",rank,emin_for_current_rank, emax_for_current_rank);
//    }



//    if(rank==0){
//        for(int i=0;i<intervalsNum;i++)
//            printf("#!intervalsE[%d]=%f\n",i,intervalsE[i]);

//    }

    tdist=(int)size/(intervalsNum/2);   // по скольку процессов на каждый интервал
    //printf("\n#!!!my rank=%d, intEmin=%d, intEmax=%d",rank,(int)((rank)/tdist)*2,(int)((rank)/tdist)*2+1);

    emin_for_current_rank=intervalsE[(int)((rank)/tdist)*2];
    emax_for_current_rank=intervalsE[(int)((rank)/tdist)*2+1];
    //printf("\n#!!!my rank=%d, emin=%f    emax=%f\n",rank,emin_for_current_rank, emax_for_current_rank);

    ////////////////

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

    //fflush(stdout);
    mc(emin,emax);
    MPI_Barrier(MPI_COMM_WORLD);

    normalize();

    //вывод
    showResult();
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

/// Полный пересчет суммарной энергии системы, для актуализации
void recalcE(){
    unsigned i,j,is,js;
    e=0;
    for (i=0; i<n; ++i){
        for (j=sequencies[i]; j<sequencies[i]+a_neighbours[i]; ++j){
            is=spins[i];
            js=spins[neighbours[j]];
            if (is!=js)
                e-=energies[j];
            else
                e+=energies[j];
        }
    }
    e/=2.;
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
  unsigned check,flag1,flag2,reduce_flag;
  long long step, totalstep;
  int count;
  double sum;
  unsigned iterator_for_exchange,rand_for_exchange,rand_for_exchange2;


/*   initialization  */
  totalstep=0;
  f=1;

  for(ie=0; ie<histSize; ie++){
    nonzero[ie]=0;
  }
  flag2=false;
  tt=0;
  do{

    flag1=0;
    step=0;

    for(ie=0; ie<histSize; ie++){
      visit[ie]=0;
    }

    while(flag1 == 0){

      single(eFrom,eTo);


      step++;



      if (step%1000==0){              // каждые 10000 пересчитываем суммарную энергию
          recalcE();
          if(flag2==false)
          {
              MPI_Barrier(MPI_COMM_WORLD);
              if(tdist==1){
                  for(iterator_for_exchange=0;iterator_for_exchange<(intervalsNum/2)-1;++iterator_for_exchange){
                      //printf("\n\n My rank = %d, Exchange(%d,%d)",rank,iterator_for_exchange,iterator_for_exchange+1);
                      exchange(iterator_for_exchange,iterator_for_exchange+1);
                      MPI_Barrier(MPI_COMM_WORLD);
                  }
              }
              else{
                  for(iterator_for_exchange=0;iterator_for_exchange<(intervalsNum/2-1);++iterator_for_exchange){
                      if(rank==0)
                      {
                          rand_for_exchange=rand()%tdist;
                          rand_for_exchange2=rand()%tdist;
                      }
                      MPI_Bcast(&rand_for_exchange,1, MPI_INT, 0, MPI_COMM_WORLD);                 // рассылаем номера обмен. ранков
                      MPI_Bcast(&rand_for_exchange2,1, MPI_INT, 0, MPI_COMM_WORLD);
                      //printf("\n\n My rank = %d, Exchange(%d,%d)",rank,tdist*iterator_for_exchange+rand_for_exchange,tdist*(iterator_for_exchange+1)+rand_for_exchange2);
                      exchange(tdist*iterator_for_exchange+rand_for_exchange,tdist*(iterator_for_exchange+1)+rand_for_exchange2);
                      MPI_Barrier(MPI_COMM_WORLD);
                  }
              }

              MPI_Allreduce(&tt,&reduce_flag,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
              if(reduce_flag>=nfinal)
                  flag2=true;
          }

      }
      //MPI_Barrier(MPI_COMM_WORLD);



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

        if(check==1){flag1++;}
      }
    }

    gupdate();

    totalstep += step;



    printf("# My rank = %d, n=%d    MCS=%lld\n",rank,tt,totalstep);
    fflush(stdout);

    f = f/2;
    tt++;
  } while(flag1*flag2==false);
  printf("\nMy rank = %d, I Am Here\n",rank);
  printf("# final   MCS=%lld\n",totalstep);
  fflush(stdout);


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

        if(fabs(e)<1E-11){
            e=0;
            enKey = (int)((e-emin)*PRECISION);
        }


        ga = g[eoKey];          // g[старой энергии]
        gb = g[enKey];          // g[новой энергии]

        if(e < eFrom || e > eTo){      // условия переворота, если не принимаем, то заходим внутрь цикла
            spins[la] *= -1;        // не принимаем новую конфигурацию, обратно переворачиваем спин
            e = energyOld;          // обратно записываем старую энергию
            enKey = eoKey;          // берем старый столбик гистограммы
        }
        else{
            if(exp(ga-gb) <= (double)rand()/RAND_MAX){      // условия переворота, если не принимаем, то заходим внутрь цикла
                spins[la] *= -1;        // не принимаем новую конфигурацию, обратно переворачиваем спин
                e = energyOld;          // обратно записываем старую энергию
                enKey = eoKey;          // берем старый столбик гистограммы
            }
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



bool exchange(int a, int b){

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

    if (rank == (int)a)
     {
         current_energy = (int)((e-emin)*PRECISION);
         //printf("#1Smy rank=%d   I send my e_a=%f    g_a[e_a]=%f\n", rank,e, g[current_energy]);//debug

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
             //printf("#2Rmy rank=%d   exch_E=%f   g_b[e_b]=%f g_b[e_a]=%f\n", rank,exchange_energy, exchange_Ge_b,exchange_Ge_a);//debug

             current_energy2=(int)((exchange_energy-emin)*PRECISION);
             exchange_probobility =(g[current_energy]*exchange_Ge_b)/(exchange_Ge_a*g[current_energy2]);
             if(exchange_probobility<1)
                 exchange_probobility_final=exchange_probobility;
             else
                 exchange_probobility_final=1;
             //printf("#3My rank=%d   exchange_probobility =(g[current_energy] =  %f   * exchange_Ge_b=%f)  / exchange_Ge_a=%f * g[current_energy2]=%f = %f\n", rank,g[current_energy], exchange_Ge_b,exchange_Ge_a,g[current_energy2],exchange_probobility);//debug
             exchange_rand=(double)rand()/RAND_MAX;

             //2st check to exit
             if(exchange_probobility_final<exchange_rand){  // fail
                 exchange_status=0;
                 MPI_Bcast(&exchange_status,1, MPI_BYTE, a, MPI_COMM_WORLD);    // recive signal, 0 -> exit.
                 //printf("2nd Cancel,exchange_probobility_final = %f    exchange_rand = %f",exchange_probobility_final,exchange_rand);//debug
                 return 0;
             }
             else{  //continue
                 //printf("Succses,exchange_probobility_final = %f    exchange_rand = %f",exchange_probobility_final,exchange_rand);//debug
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



         //printf("#1Rmy rank=%d   My e_b=%f g_b[e_b]=%f, and I recive e_a=%f    g_a[e_a]=%f\n", rank,e,g[current_energy],exchange_energy,exchange_Ge_a);//debug

         //1st проверка на поподание энергий в окно
         if(exchange_energy>emax_for_current_rank||exchange_energy<emin_for_current_rank||e<exchange_emin||e>exchange_emax){ // если не попала в интервал энергии текущего процесса, то выйти из ф-ции
             exchange_status=0;
             MPI_Bcast(&exchange_status,1, MPI_BYTE, b, MPI_COMM_WORLD);    // рассылка сигнала exchange_status =0 выхода из фунуции.
             //printf("\n1st Cancel\n");
             return 0;
         }

         else{
             MPI_Bcast(&exchange_status,1, MPI_BYTE, b, MPI_COMM_WORLD);    // рассылка exchange_status=1 все ок, продолжаем

             current_energy2 = (int)((exchange_energy-emin)*PRECISION); // g_b(e_a)

             //printf("#2Smy rank=%d I send my exch_E=%f   g_b[e_b]=%f g_b[e_a]=%f\n", rank,e, g[current_energy],g[current_energy2]);//debug
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

void showResult(){
    MPI_Status status;


    MPI_Barrier(MPI_COMM_WORLD);
    // Вывод сделан именно в таком формате для совместимости с программой склейки гистограмм
    printf("# intervals=%d\n",intervalsNum);
    printf("# gaps=%d\n",histSize);
    printf("# walkers=%d\n",size);
    printf("# nfrom=");
    for (unsigned i=0;i<intervalsNum;i+=2){
        printf("%e,",intervals[i]);
    }
    printf("\n");
    printf("# nto=");
    for (unsigned i=1;i<intervalsNum;i+=2){
        printf("%e,",intervals[i]);
    }
    printf("\n");

    if (rank==0){
        double *g2, *hist2;
        int *nonzero2;
        g2 = (double *) malloc(histSize*sizeof(double));
        hist2 = (double *) malloc(histSize*sizeof(double));
        nonzero2 = (int *) malloc(histSize*sizeof(int));
        for (int i=0; i < size; ++i){
            printf("-----\n");
            printf("%d\n",i);

            if (i==rank){
                for (unsigned j=0; j<histSize; ++j){
                    g2[j]=g[j];
                    hist2[j]=hist[j];
                    nonzero2[j]=nonzero[j];
                }
            } else {
                MPI_Recv(&g2,histSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
                MPI_Recv(&hist2,histSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
                MPI_Recv(&nonzero2,histSize,MPI_INT,i,0,MPI_COMM_WORLD,&status);
            }

            for (unsigned ie=0;ie<histSize;++ie){
                if (nonzero2[ie] == 1)
                    printf("%e  %e  %e  %d\n",(double)ie/PRECISION+emin,g[ie],g[ie]/n,hist[ie]);
            }
        }
        free(g2);
        free(hist2);
        free(nonzero2);
    } else {
        MPI_Send(&g,histSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        MPI_Send(&hist,histSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
