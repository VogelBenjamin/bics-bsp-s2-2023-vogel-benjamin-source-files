#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "linearAlgebra.h"
#include "omp.h"


int vectorSize = 3;

/*
  creates and returns a new vector of size 3;
*/
double* createVector(){
  double* newVector = (double*)malloc(vectorSize*sizeof(newVector));
  return newVector;
}

void freeMatrix(double **matrix, int col){
  for (int i = 0; i < col; ++i)
  {
    free(matrix[i]);
  }
  free(matrix);
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
    #pragma omp atomic read
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
  norm = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
  return norm;
}

double distance(double* pos1, double* pos2){
  double* vec = connectPoints(pos1,pos2);
  double distance = norm(vec);
  free(vec);
  return distance;
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
  adds vector2 to vector1
*/
void vectorAtomicAddition(double *vector1, double *vector2){
  for (int i = 0; i < vectorSize; ++i)
  {
    #pragma omp atomic
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
  int i;
  for (i = 0; i < vectorSize; ++i)
  {
    vector[i] *= scalar;
  }
}

/*
  calculates and returns the dot product of two given vectors
*/
double dotProduct(double *vector1, double *vector2){
  double dP = 0;
  int i;
  #pragma omp parallel shared(vectorSize) private(i) reduction(+: dP)
    for (int i; i < vectorSize; ++i)
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
  if (vNorm == 0)
  {
    return;
  }
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

/*
  searches for unique solution of an equation
  assuming equation is linear
*/
double *solveAugMatrix(double **augMatrix, int numEq, int numVar){
  double ratio, sub;
  double *solution = (double*)malloc(sizeof(solution)*numVar);
  double temp;
  int i,j,k;
  // transform matrix into reduced row echelon form
  // int i decides which column to operate on
  for (i = 0; i < numVar; ++i)
  {
    // makes sure that augementmatrix[i][i] is not 0
    if (augMatrix[i][i] == 0)
    {

      for (j = i+1; j < numEq; ++j)
      {

        if (augMatrix[i][j] != 0)
        {

          for (k = 0; k < numVar+1; ++k)
          {

            temp = augMatrix[k][i];
            augMatrix[k][i] = augMatrix[k][j];
            augMatrix[k][j] = temp;
          }

          break;
        }
      }
    }



    if (augMatrix[i][i] != 0)
    {
      // makes sure there is a leading one within the row
      ratio = 1.0/(augMatrix[i][i]);
      for (j = 0; j < numVar+1; ++j)
      {
        augMatrix[j][i] *= ratio;
      }

      // makes sure the values above and below the leading are turned to 0
      // j indicates which row is operated on
      for (j = 0; j < numEq; ++j)
      {
        if (j != i)
        {
          sub = augMatrix[i][j];
          // k indicates on which element of the row is operated on
          for (k = 0; k < numEq; ++k)
          {
            augMatrix[k][j] -= sub*augMatrix[k][i];
          }
        }
      }
    }

  }

  // determine what the solution is

  if (!(abs(augMatrix[numVar][numEq-1]) < 1E-10))
  {
    // no solution -> return an array containing only NAN
    for (int i = 0; i < numVar; ++i)
    {
      solution[i] = NAN;
    }
  } else {
    // the solution is contained within the column containing the equantion constants
    for (int i = 0; i < numVar; ++i)
    {
      solution[i] = augMatrix[numVar][i];
    }
  }

  return solution;
}

