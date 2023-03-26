#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"linearAlgebra.h"
#include "omp.h"


int vectorSize = 3;

/*
  creates and returns a new vector of size 3;
*/
double* createVector(){
  double* newVector = (double*)malloc(vectorSize*sizeof(newVector));
  return newVector;
}

/*
  creates and returns a new vector that connects two points pos1 and pos2 in 3D space
*/
double* connectPoints(double* pos1, double* pos2){
  double* newVector = (double*)malloc(vectorSize*sizeof(newVector));
  for (int i = 0; i < vectorSize; ++i)
  {
    newVector[i] = pos2[i] - pos1[i];
  }
  return newVector;
}

/*
  prompts the user to input values for a vector
*/
void setVector(double *vector){
  printf("Please input the vector values (x,y,z): ");
  for (int i = 0; i < vectorSize; ++i)
  {
    scanf("%lf",vector+i);
  }
}

/*
  modifies the values of an existing vector 'oldVector' tom atch the values of a given vector
*/
void modifyVector(double *oldVector, double *newVector){
  for (int i = 0; i < vectorSize; ++i)
  {
    oldVector[i] = newVector[i];
  }
}

/*
  prints the values of a given vector
*/
void printVector(double *vector){
  printf("(");
  for (int i = 0; i < vectorSize-1; ++i)
  {
    printf("%lf,", vector[i] );
  }
  printf("%lf",vector[2]);
  printf(")\n");
}

/*
  calculates and retuns the Euclidean norm of a given vector
*/
double norm(double* vector){
  double norm = 0;
  norm = sqrt(pow(vector[0],2) + pow(vector[1],2) + pow(vector[2],2));
  return norm;
}

/*
  adds vector2 to vector1
*/
void vectorAddition(double *vector1, double *vector2){
  for (int i = 0; i < vectorSize; ++i)
  {
    vector1[i] += vector2[i];
  }
}

/*
  subtracts vector2 from vector 1
*/
void vectorSubtraction(double *vector1, double *vector2){
  for (int i = 0; i < vectorSize; ++i)
  {
    vector1[i] -= vector2[i];
  }
}

/*
  multiplies each element of a given vector by a scalar value
*/
void scalarMultiplication(double scalar, double *vector){
  for (int i = 0; i < vectorSize; ++i)
  {
    vector[i] *= scalar;
  }
}

/*
  calculates and returns the dot product of two given vectors
*/
double dotProduct(double *vector1, double *vector2){
  double dP = 0;
  for (int i = 0; i < vectorSize; ++i)
  {
    dP += vector1[i]*vector2[i];
  }
  return dP;
}

/*
  calculates the cross product of two given vectors and stores its values in the vector cP
*/
void crossProduct(double *vector1, double *vector2, double* cP){
  cP[0] = vector1[1]*vector2[2]-vector2[1]*vector1[2];
  cP[1] = -vector1[0]*vector2[2]+vector2[0]*vector1[2];
  cP[2] = vector1[0]*vector2[1]-vector2[0]*vector1[1];
}

/*
  normalizes a given vector
*/
void normalize(double *vector){
  double vNorm = norm(vector);
  for (int i = 0; i < vectorSize; ++i)
  {
    vector[i] /= vNorm;
  }
}

/*
  gives the angle between two vectors
*/
double angleBetween(double *v1, double * v2){
  return acos((dotProduct(v1,v2))/(norm(v1)*norm(v2)));
}