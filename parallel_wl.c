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


int n;                          // количество спинов
signed char *spins;             //массив направления спинов. По умолчанию +1. n - число считанных спинов (число незакомментированных строк в csv-файле).
unsigned short *a_neighbours;   //число соседей каждого спина. Считается как число энергий в соответствующей строке в csv-файле.
unsigned short *neighbours;     // соседи каждого спина
unsigned int *sequencies;       //для каждого спина описывает, с какого ключа в массиве energies[] начинают описываться парные энергии
double *energies;               //сами энергии из файла. Описывается как одномерный массив. Длина массива - число парных энергий в csv-файле.
float emin, emax;               //минимумы и максимумы энергии
float e;                        //текущая энергия системы
unsigned eCount=0;              //число пар энергий

double *g;
double *visit;
double *hist;
int *nonzero;

double f,factor;
int nfinal;



#define PRECISION 3             //Сколько знаков учитывать в энергии после запятой !! НЕ СДЕЛАНО ЕЩЕ


void readCSV(char* filename);
void rotate(int spin);          // Считает энергию системы
void complete();

void mc();
void single();
void gupdate();
void normalize();


int main(void)
{
    int seed=0;                 // Random seed
    printf("# Please, input random number seed:  ");
    scanf("%u",&seed);
    char filename[100];
    printf("# Please, input target filename: ");
    scanf("%s",filename);

    readCSV(filename);

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
    
    
    
    
    
    
    
    printf("# initial energy = %lf\n",e);    // modified
    
    //*
    
    factor = 0.8;
    nfinal = 24; //изменить, не понял что это
    
    int ie;
    for(ie=0; ie<=eCount; ie++){
      g[ie]=0;
      hist[ie]=0;
    }
    
    srand(seed);
    mc();
    normalize();

    for(ie=0; ie<=eCount; ie++){   
      if (nonzero[ie] == 1) {
        printf("%d  %e  %e  %e\n",ie-3*n,g[ie],g[ie]/n,hist[ie]); 
      }
    }
    //*/
    
    
    complete();
}


void readCSV(char *filename){

    char c;                         //считанный из файла символ
    char symb[100];                 //символ энергии в текстовом файле

    //get system sizes
    bool isFirstLine=true;
    n=0;
    FILE *file2 = fopen(filename, "r");
    int fpos = 1, lastFpos=0;

    while(c = fgetc(file2)=='#')     //пропуск комментариев
           {
                fscanf(file2,"%[^\n]%*c",symb);
           }
     fseek(file2,-1,SEEK_CUR);       // сдвиг курсора на один символ назад
     int coursor=ftell(file2);       // положение курсора начала данных

    do{
        c = fgetc(file2);
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
    neighbours=(unsigned short *) malloc(eCount*sizeof(unsigned short));     //поменять размер
    sequencies=(unsigned int *) malloc(n*sizeof(unsigned int));
    energies = (double *) malloc(eCount*sizeof(double));                        //поменять размер
    
    g = (double *) malloc(eCount*sizeof(double)); //не уверен по поводу размеров
    visit = (double *) malloc(eCount*sizeof(double)); //не уверен по поводу размеров
    hist = (double *) malloc(eCount*sizeof(double)); //не уверен по поводу размеров
    nonzero = (double *) malloc(eCount*sizeof(double)); //не уверен по поводу размеров


    // read data

    fseek(file2,coursor,SEEK_SET);      //устанавливаем курсор в начало данных

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
        c = fgetc(file2);

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


        if (c=='\n'){ //if it is newline, mark the flag
            firstSymbolInLine=true;
//            skipFlag=false;
        }
    } while (c != EOF);

    emax/=2;
    e/=2;
    emin = -emax;

    fclose(file2);
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




void mc()
/*
        monte carlo update
*/
{
  int ie,n;
  int check,flag;
  int step, totalstep;
  int count;
  double sum;

/*   initialization  */
  totalstep=0;
  f=1;

  for(ie=0; ie<=eCount; ie++){
    nonzero[ie]=0;
  }

  for( n = 0; n <= nfinal; n++){

    flag=0;
    step=0;

    for(ie=0; ie<=eCount; ie++){ 
      visit[ie]=0;
    }

    while(flag == 0){

      single();

      step++;

      if(step%1000==0){

        for(ie=0; ie<=eCount; ie++){
          if(visit[ie] > 0) {nonzero[ie]=1;}
        }

        count=0;
        sum=0;
        for(ie=0; ie<=eCount; ie++){
          if(nonzero[ie]==1) {
            count++;
            sum+=visit[ie];
          }
        }

        check=1; 
        for(ie=0; ie<=eCount; ie++){
          if(nonzero[ie]==1) {
            if(visit[ie] < factor*(sum/count)){check=0;}
          }
        }

        if(check==1){flag++;}
      }
    }

    gupdate();

    totalstep += step;

    printf("# n=%2d    MCS=%9d\n",n,totalstep);

    f = f/2;
  }
  printf("# final   MCS=%9d\n",totalstep);

}

void single()
/*   single spin flip */
{
    int la,la1;
    double energyOld;
    double ga,gb;
    
    
    //проверить весь алгоритм!!!!!!!!!!!!!!!!
    for(la1=0; la1 <= n-1; la1++){//нужен ли этот цикл?
        la=rand()%n;
        energyOld = e;
        rotate(la);
        
        ga = g[(energyOld-emin)*1ePRECISION]; //проверить ключи
        gb = g[(e-emin)*1ePRECISION]; //проверить ключи
        
        if(exp(ga-gb) > rand()/RAND_MAX){
            rotate(la);
            energyOld  = e;
        }
        
        g[(energyOld-emin)*1ePRECISION]     += f; //проверить ключи
        visit[(energyOld-emin)*1ePRECISION] += 1; //проверить ключи
        hist[(energyOld-emin)*1ePRECISION]  += 1; //проверить ключи
    }
}

void gupdate()
{
    int ie;
    double gmin;
    
    /* set min of g[ie] as 1 */
    gmin=10000000;
    
    for (ie=0; ie<=eCount; ie++){
        if (nonzero[ie] == 1) {
            if(g[ie] < gmin) {
                gmin = g[ie];
            }
        }
    }
    
    for (ie=0; ie<=eCount; ie++){
        if (nonzero[ie] == 1) {
            g[ie] += -gmin;
        }
    }
}

void normalize()
{
    int ie;
    double gmax, sum, a;
    
    gmax = -1000;
    for(ie=0; ie<eCount; ie++){
        if (nonzero[ie] == 1) {
            if(g[ie]>gmax){
                gmax = g[ie];
            }
        }
    }
    
    sum=0;
    for(ie=0; ie<eCount; ie++){
        if (nonzero[ie] == 1) {
            sum += exp(g[ie]-gmax);
        }
    }
    
    a = n*log(2) - gmax - log(sum);
    
    for(ie=0; ie<eCount; ie++){
        if (nonzero[ie] == 1) {
            g[ie] += a;
        }
    }
}

