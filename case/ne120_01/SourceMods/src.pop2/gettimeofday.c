//tao liu @20180508
#include <time.h>   
#include <stdio.h>   
#include<sys/time.h>  
#include <unistd.h>

//interface for C:
void getlooptime(double *mytime)
{
    struct timeval Time;
    gettimeofday(&Time,NULL);
    *mytime = (double)Time.tv_sec*1000000 + (double)Time.tv_usec;
}

//interface for fortran:
void getlooptime_(double* mytime)
{
    struct timeval Time;
    gettimeofday(&Time,NULL);
    *mytime = (double)Time.tv_sec*1000000 + (double)Time.tv_usec;
}
