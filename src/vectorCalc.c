#include<stdlib.h>
#include<stdio.h>
#include"linearAlgebra.h"

int main(){

  double **arrayVector;
  double *vector1, *vector2;
  int numberOfVectors;

  printf("Please input the number of vectors: ");
  scanf("%d",&numberOfVectors);

  arrayVector = (double**)malloc(numberOfVectors*sizeof(arrayVector));

  arrayVector[0] = vector1;
  arrayVector[1] = vector2;

  for (int i = 0; i < numberOfVectors; ++i)
  {
    arrayVector[i] = createVector();
    setVector(arrayVector[i]);
    printVector(arrayVector[i]);
  }

  return 0;
}


