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


unsigned n;                          // количество спинов
signed char *spins;             //массив направления спинов. По умолчанию +1. n - число считанных спинов (число незакомментированных строк в csv-файле).
unsigned short *a_neighbours;   //число соседей каждого спина. Считается как число энергий в соответствующей строке в csv-файле.
unsigned short *neighbours;     // соседи каждого спина
unsigned int *sequencies;       //для каждого спина описывает, с какого ключа в массиве energies[] начинают описываться парные энергии
double *energies;               //сами энергии из файла. Описывается как одномерный массив. Длина массива - число парных энергий в csv-файле.
double *intervals;              //массив интервалов
double *intervalsE;             //массив интервалов значений
int intervalsNum=0;             //число значений интервалов из файла
double emin, emax;               //минимумы и максимумы энергии
double e;                        //текущая энергия системы
unsigned eCount=0;              //число пар энергий
unsigned histSize=0;            //число элементов в гистограммах

double *g;                      // гистограмма g[E], только хранятся !логарифмы!.
unsigned *visit;
unsigned *hist;
int *nonzero;

double f;           // Модификационный фактор (уменьшается с каждым WL-шагом)
double factor;      // Критерий плоскости гистограммы H
unsigned nfinal;    // число WL-циклов

#define PRECISION 1e0             // Точность 1eX, где X - Сколько знаков учитывать в энергии после запятой
                                  // (1e0 - 0 знаков после запятой (для модели Изинга), 1e100 - 100 знаков после запятой)


int readCSV(char* filename);
int readCSVintervals(char* filename);
void rotate(int spin);          // Считает энергию системы
void complete();

void mc();
void single();
void gupdate();
void normalize();
void dumpArrays();


int main(void)
{
    unsigned long seed=0;                 // Random seed
    printf("# Please, input random number seed from 1 to 4 294 967 295:  ");
    scanf("%u",&seed);
    char filename[100];
    printf("# Please, input target filename: ");
    scanf("%s",filename);

    if (!readCSV("csv_examples/square_ising_4x4.csv")){
        printf("# Error! File '%s' is unavaliable!\n",filename);
        return 0;
    }

    
    char intervalsFile[50] = "csv_examples/intervals.csv";          ///new
    
    if (!readCSVintervals(intervalsFile)){                          ///new
        printf("# Error! File '%s' is unavaliable!\n", intervalsFile);
        return 0;
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



    printf("\n# e = %lf, emin = %lf, emax = %lf\n",e,emin,emax);

    if (false){         // !DEBUG если true - загнать модель изинга в минимум
        rotate(1);
        rotate(3);
        rotate(4);
        rotate(6);
        rotate(9);
        rotate(11);
        rotate(12);
        rotate(14);
    }

    printf("\n# e = %lf\n",e);
    
    printf("# initial energy = %lf\n",e);
    
    //*
    
    factor = 0.8;       // Критерий плоскости гистограммы H
    nfinal = 24;        // число WL-циклов
    
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

    
    
     complete(); //  все пашет)
}

/// Функция чтения файла с энергиями
int readCSV(char *filename){

    char c;                         //считанный из файла символ
    char symb[100];                 //символ энергии в текстовом файле

    //get system sizes
    bool isFirstLine=true;
    n=0;
    FILE *file = fopen(filename, "r");

    if (!file)
        return 0;

    int fpos = 1, lastFpos=0;

    while(c = fgetc(file)=='#')     //пропуск комментариев
           {
                fscanf(file,"%[^\n]%*c",symb);
           }
     fseek(file,-1,SEEK_CUR);       // сдвиг курсора на один символ назад
     int coursor=ftell(file);       // положение курсора начала данных

    do{
        c = fgetc(file);
//        while(c=='#'){
//            do c = fgetc(file2); while (c != '\n');           // нет необходимости, только если у нас не будет комментариев прямо посреди данных, но можно оставить
//            c = fgetc(file2);
//        }
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

    // reserve memory for arrays
    spins=(signed char *) malloc(n*sizeof(signed char));
    a_neighbours=(unsigned short *) malloc(n*sizeof(unsigned short));
    neighbours=(unsigned short *) malloc(eCount*sizeof(unsigned short));
    sequencies=(unsigned int *) malloc(n*sizeof(unsigned int));
    energies = (double *) malloc(eCount*sizeof(double));


    // read data

    fseek(file,coursor,SEEK_SET);      //устанавливаем курсор в начало данных

    bool firstSymbolInLine=true, skipFlag=false;
    double parsedNumber;
    int numInSymb=0;
    symb[0]='\0';
    int row=0;                  //line number in file (not account the commented lines)
    int col=0;                  //column number in line (taking to accound the ';' symbols)
    int neighCount=0;           //
    int energyNum=0;            //holds actual count of previously parsed energies
    e = 0;
    emax = 0;                   // сумма всех взаимодействий с положительным знаком

    do {
        c = fgetc(file);

//        if (firstSymbolInLine && c=='#'){ //if it is comment, skip the line
//            skipFlag=true; //skip to end of line
//        }                                           // нет необходимости, только если у нас не будет комментариев прямо посреди данных, но можно оставить
        firstSymbolInLine=false;

//        if (!skipFlag){
            if (c==';' || c=='\n' || c == EOF){ //if we found a number, process it
                if (numInSymb!=0){
                    sscanf( symb, "%lf", &parsedNumber );
                    neighbours[energyNum] = col;
                    energies[energyNum] = parsedNumber;
                    e += parsedNumber;
                    emax += fabs(parsedNumber);

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
            }


        if (c=='\n'){ //if it is newline, mark the flag
            firstSymbolInLine=true;
//            skipFlag=false;
        }
    } while (c != EOF);

    emax/=2;
    e/=2;
    emin = -emax;

    fclose(file);

    histSize = (int)((emax-emin)*PRECISION)+1;              // почему резервирование этих массивов происходит именно в этой функции?
    g = (double *) malloc(histSize*sizeof(double));
    visit = (unsigned *) malloc(histSize*sizeof(unsigned));
    hist = (unsigned *) malloc(histSize*sizeof(unsigned));
    nonzero = (int *) malloc(histSize*sizeof(int));

    return 1;
}
/// Переворот спина, подсчет изменения энергии
void rotate(int spin){
    double dE=0;
    spins[spin] *= -1;
    for(unsigned i = sequencies[spin]; i<sequencies[spin]+a_neighbours[spin]; ++i){
        dE += energies[i]*spins[neighbours[i]]*spins[spin]*2;
    }
    e += dE;
}
/// Очистка памяти
void complete(){
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
void mc()
/*
        monte carlo update
*/
{
  unsigned ie,tt; //итераторы
  int check,flag;
  int step, totalstep;
  int count;
  double sum;

/*   initialization  */
  totalstep=0;
  f=1;

  for(ie=0; ie<=histSize; ie++){    //обнуляем массив nonzero
    nonzero[ie]=0;
  }

  for( tt = 0; tt <= nfinal; tt++){    // !!! Кто придумал в качестве итератора WL-шагов использовать n? , я исправил на tt

    flag=0;
    step=0;

    for(ie=0; ie<=histSize; ie++){
      visit[ie]=0;                  // обнуляем массив visit
    }

    while(flag == 0){

      single();

      step++;

      if(step%1000==0){         // каждые 1000 шагов

        for(ie=0; ie<=histSize; ie++){
          if(visit[ie] > 0) {nonzero[ie]=1;}        // проверяем, появились ли новые энергии
        }

        count=0;
        sum=0;
        for(ie=0; ie<=histSize; ie++){
          if(nonzero[ie]==1) {
            count++;                        // подсчитываем количество встреченных энергий
            sum+=visit[ie];                 // сумма посещений всех энергий
          }
        }

        check=1; 
        for(ie=0; ie<=histSize; ie++){                      // проверка на плоскоту
          if(nonzero[ie]==1) {
            if(visit[ie] < factor*(sum/count)){check=0;}    // sum/count = среднее значение количества посещений энергий
          }                                                 // если количество посещений хоть одной энергии меньше среднего значения * factor, то проверка провалилась
        }

        if (false && step%100000) //написать true для дебаг-вывода в файл
            dumpArrays();

        if(check==1){flag++;}
      }
    }

    gupdate();

    totalstep += step;

    printf("# n=%2d    MCS=%9d\n",tt,totalstep);    // ! на самом деле тут totalstep*n MCS, так как в функции single цикл по n

    f = f/2;    // !! НЕОБХОДИМО ИСПРАВИТЬ НА t, при чем t пропорционально шагу моделирования !!  doi:10.1063/1.2803061
  }
  printf("# final   MCS=%9d\n",totalstep);

}

void single()
/*   single spin flip */    // нифига не сингл флип, а n spins flips.
{
    unsigned la,la1;        // итераторы
    double energyOld;       // старая энергия
    double ga,gb;           // g[старой энергии] и g[новой энергии]

    int eoKey, enKey;       // номер столбика гистограммы энергий старой и новой
    
    
    //проверить весь алгоритм!!!!!!!!!!!!!!!!
    for(la1=0; la1 <= n-1; la1++){  //цикл выполняется n раз, не знаю почему
        la=rand()%n;            // выбираем случайный спин
        energyOld = e;          // записываем старую энергию
        rotate(la);             // переворачиваем выбранный спин

        eoKey = (int)((energyOld-emin)*PRECISION); //вычисляем номер столбика гистограммы для старой энергии
        enKey = (int)((e-emin)*PRECISION);         //вычисляем номер столбика гистограммы для новой энергии
        
        ga = g[eoKey];          // g[старой энергии]
        gb = g[enKey];          // g[новой энергии]
        
        if(exp(ga-gb) <= (double)rand()/RAND_MAX){      // условия переворота, если не принимаем, то заходим внутрь цикла
            spins[la] *= -1;        // не принимаем новую конфеигурацию, обратно переворачиваем спин
            e = energyOld;          // обратно записываем старую энергию
            enKey = eoKey;          // берем старый столбик гистограммы
        }
        
        g[enKey]     += f;          // прибавляем f в текущий столбик гистограммы (так как тут хрянятся логарифмы)
        visit[enKey] += 1;
        hist[enKey]  += 1;
    }
}

// нормализация g[E]
void gupdate()
{
    /* set min of g[ie] as 1 */
    double gmin=10000000;              // устанавливаем gmin очень большим, чтобы в системе точно найти минимум.

    for (unsigned ie=0; ie<=histSize; ++ie){
        if (nonzero[ie] == 1) {
            //printf("!# g[%u]=%f\n",ie,g[ie]);
            if(g[ie] < gmin) {
                gmin = g[ie];   // находим наименьшее значение g[E]

            }
        }
    }

    for (unsigned ie=0; ie<=histSize; ++ie){
        if (nonzero[ie] == 1) {
            g[ie] += -gmin;     // !! нормализация g[E], а тут =0, так как тут хранятся логарифмы.
            //printf("!# g[%u]=%f\n",ie,g[ie]);
        }
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

void dumpArrays(){
    FILE *file = fopen("dump.dat", "w");

    fprintf(file,"E=%e; state=",e);
    for(unsigned ie=0; ie<n; ie++){
        fprintf(file,"%d",spins[ie]);
    }
    fprintf(file,"\n");

    fprintf(file,"ie  E  g[ie]  g[ie]/n  hist[ie]  visit[ie]\n");
    for(unsigned ie=0; ie<=histSize; ie++){
      if (nonzero[ie] == 1) {
        fprintf(file,"%d  %e  %e  %e  %d  %d\n",ie,(double)(ie+emin)/PRECISION,g[ie],g[ie]/n,hist[ie],visit[ie]);
      }
    }

    fclose(file);
}

int readCSVintervals(char *filename){
    int numerOfStrings = 0;
    char c;                         //считаный из файла символ
    char symb[100];                 //символ энергии в текстовом файле
    
    //get system sizes
    bool isFirstLine=true;
    FILE *file = fopen(filename, "r");
    
    if (!file)
        return 0;
    
    while(fgetc(file)=='#')     //пропуск комментариев и пустой строки
    {
        fscanf(file,"%[^\n]%*c",symb);
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
    
    for(int i=0; i<numerOfStrings*2; ++i)
        intervalsE[i] = (emax-emin)*intervals[i];
    
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
    
    return 1;
}
