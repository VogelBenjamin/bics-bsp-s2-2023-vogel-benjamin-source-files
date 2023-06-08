#include "omp.h"
#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<dirent.h>
#include<math.h>
#include<time.h>
#include "extra.h"
void main(int argc, char *argv[]){

  #pragma omp parallel
    {
      printf("WOW: %s", argv[1]);
      printWhat();
      printf("Hello World... from thread = %d\n",
           omp_get_thread_num());
    }
}
