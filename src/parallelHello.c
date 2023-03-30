#include "omp.h"
#include<stdio.h>
void main(){

  #pragma omp parallel
    {
      printf("Hello World... from thread = %d\n",
           omp_get_thread_num());
    }
}