#include<stdio.h>
#include<stdlib.h>
#include "omp.h"

void setVector(int* vector, int vecSize);
void addVector(int *vector1, int *vector2, int vecSize);
void printVector(int *vector, int vecSize);

int main(){
  int *vector1, *vector2;
  int vecSize;
  omp_set_num_threads(3);

  // Takes size of vector as input
  printf("Please input the size of the two vectors: ");
  scanf("%d",&vecSize);
  vector1 = (int*)malloc(vecSize*sizeof(vector1));
  vector2 = (int*)malloc(vecSize*sizeof(vector2));

  // takes vector values
  setVector(vector1, vecSize);
  setVector(vector2, vecSize);

  // adds the vectors together
  addVector(vector1, vector2, vecSize);

  // print vector
  printVector(vector1, vecSize);
}

void setVector(int* vector, int vecSize){
  printf("Please input the values of the vector: ");
  for (int i = 0; i < vecSize; ++i)
  {
    scanf("%d",vector+i);
  }
}

void addVector(int *vector1, int *vector2, int vecSize){
  int i;
  // #pragma omp parallel for shared(vecSize,vector1,vector2) private(i)
  for (i = 0; i < vecSize; ++i)
  {
    vector1[i] += vector2[i];
  }
}

void printVector(int *vector, int vecSize){
  printf("Vector info: ");
  for (int i = 0; i < vecSize; ++i)
  {
    printf("%d ",vector[i]);
  }
  printf("\n");
}