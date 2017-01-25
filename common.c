/*
 * common.c
 *
 *       Common functions for sequentional and parallel WL method for different magnetic systems
 *
 *
 *  programmed by:
 *           Makarov Aleksandr
 *           Andriushchenko Petr
 *           Shevchenko Yuriy
 *
 */
// нормализация g[E]

void gupdate()
{
    /* set min of g[ie] as 1 */
    double gmin=10000000;              // устанавливаем gmin очень большим, чтобы в системе точно найти минимум.
    unsigned ie;
    for (ie=0; ie<histSize; ++ie){
#ifdef DEBUG
        if (ie>=histSize || ie<0) printf("Error with memory working5");
#endif
        if (nonzero[ie] == 1) {
            //printf("!# g[%u]=%f\n",ie,g[ie]);
            if(g[ie] < gmin) {
                gmin = g[ie];   // находим наименьшее значение g[E]

            }
        }
    }

    for (ie=0; ie<histSize; ++ie){
#ifdef DEBUG
        if (ie>=histSize || ie<0) printf("Error with memory working5");
#endif
        if (nonzero[ie] == 1) {
            g[ie] += -gmin;     // !! нормализация g[E], а тут =0, так как тут хранятся логарифмы.
            //printf("!# g[%u]=%f\n",ie,g[ie]);
        }
    }
}


void dumpArrays(){
    unsigned ie;
    FILE *file = fopen("dump.dat", "w");

    fprintf(file,"E=%e; state=",e);
    for(ie=0; ie<n; ie++){
        fprintf(file,"%d",spins[ie]);
    }
    fprintf(file,"\n");

    fprintf(file,"ie  E  g[ie]  g[ie]/n  hist[ie]  visit[ie]\n");
    for(ie=0; ie<histSize; ie++){
        if (nonzero[ie] == 1) {
            fprintf(file,"%d  %e  %e  %e  %d  %d\n",ie,(double)(ie+emin)/PRECISION,g[ie],g[ie]/n,hist[ie],visit[ie]);
        }
    }

    fclose(file);
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
    int count_n=0;

    while(c = fgetc(file)=='#')     //пропуск комментариев
    {
        if(fscanf(file,"%[^\n]%*c",symb)==1){}
        else{printf("# Error! Failed to read file!\n");}
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
        if (c=='\n'){
            isFirstLine=false;
            count_n++;
        }

        if (c==';' || c=='\n') {
            if (fpos-1 != lastFpos)
                ++eCount;
            lastFpos = fpos;
        }

        fpos++;
    } while (c != EOF);
    ++n;
    if(count_n!=n)
        printf("!!!ERROR with number of element: Number of elements in first line does not correspond with number of lines");

    // reserve memory for arrays
    spins=(signed char *) malloc(n*sizeof(signed char));
    a_neighbours=(unsigned short *) malloc(n*sizeof(unsigned short));
    neighbours=(unsigned short *) malloc(eCount*sizeof(unsigned short));
    sequencies=(unsigned int *) malloc(n*sizeof(unsigned int));
    energies = (double *) malloc(eCount*sizeof(double));


    // read data

    fseek(file,coursor,SEEK_SET);      //устанавливаем курсор в начало данных

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

        if (c==';' || c=='\n' || c == EOF){ //if we found a number, process it
            if (numInSymb!=0){
                sscanf( symb, "%lf", &parsedNumber );
                neighbours[energyNum] = col;
#ifdef DEBUG
                if (energyNum>=eCount || energyNum<0) printf("Error with memory working8");
#endif
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

        if (c=='\n'){
#ifdef DEBUG
            if (row>=n || row<0)
                printf("Error with memory working9");
#endif
            a_neighbours[row] = neighCount;
            sequencies[row] = energyNum-neighCount;
            col=0;
            neighCount=0;
            spins[row]=1;
            ++row;
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
