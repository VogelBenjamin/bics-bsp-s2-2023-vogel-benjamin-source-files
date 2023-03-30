#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "omp.h"
double* createVector();
double* connectPoints(double* pos1, double* pos2);
void setVector(double *vector);
void modifyVector(double *oldVector, double *newVector);
void printVector(double *vector);
double norm(double* vector);
void vectorAddition(double *vector1, double *vector2);
void vectorSubtraction(double *vector1, double *vector2);
void scalarMultiplication(double scalar, double *vector);
double dotProduct(double *vector1, double *vector2);
void crossProduct(double *vector1, double *vector2,double* cP);
void normalize(double *vector);
double angleBetween(double *v1, double * v2);
double solve(double *equation, int eqSize);
double *getLinearDerivative(double *equation, int eqSize);