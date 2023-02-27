#include<stdio.h>
#include<stdlib.h>
#include"la.h"


int vectorSize = 3;

double* createVector(){
  double* newVector = (double*)malloc(vectorSize*sizeof(newVector));
  return newVector;
}

void setVector(double *vector){
  printf("Please input the vector values (x,y,z): ");
  for (int i = 0; i < vectorSize; ++i)
  {
    scanf("%lf",vector+i);
  }
}

void printVector(double *vector){
  printf("(");
  for (int i = 0; i < vectorSize; ++i)
  {
    printf("%lf,", vector[i] );
  }
  printf(")\n");
}

double* vectorAddition(double *vector1, double *vector2){
  double* addedVector = (double*)malloc(vectorSize*sizeof(addedVector));
  for (int i = 0; i < vectorSize; ++i)
  {
    addedVector[i] = vector1[i]+vector2[i];
  }
  return addedVector;
}

double dotProduct(double *vector1, double *vector2){
  double dP = 0;
  for (int i = 0; i < vectorSize; ++i)
  {
    dP += vector1[i]*vector2[i];
  }
  return dP;
}

double* crossProduct(double *vector1, double *vector2){
  double *cP = (double*)malloc(vectorSize*sizeof(cP));
  cP[0] = vector1[1]*vector2[2]-vector2[1]*vector1[2];
  cP[1] = -vector1[0]*vector2[2]+vector2[0]*vector1[2];
  cP[2] = vector1[0]*vector2[1]-vector2[0]*vector1[1];;
  return cP;
}