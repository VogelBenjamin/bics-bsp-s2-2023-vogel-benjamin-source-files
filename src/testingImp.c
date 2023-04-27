#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>
#include"testing.h"
#include"linearAlgebra.h"


bool compareArray(double *array1, double *array2, int size){
  for (int i = 0; i < size; ++i)
  {
    if (!(abs(array1[i]-array2[i]<0.1)))
    {
      /*
      printf("\n");
      printf("Error: array not equal!\n");
      printVector(array1);
      printVector(array2);
      */
      printf("%d: %lf is not %lf",i,array1[i], array2[i]);
      return false;
    }
  }

  /*
  printf("\n");
  printf("Correct: array equal!\n");
  printVector(array1);
  printVector(array2);
  */
  return true;
}

bool compareValue(double v1, double v2, bool *check){
  if (v1 != v2)
  {
    printf("Error: \n");
    printf("%lf is not equal to %lf \n", v1, v2 );
    *check = false;
    return false;
  }
  return true;
}

bool checkIfAllNAN(double *arr, int size){
  bool result = true;
  for (int i = 0; i < size; ++i)
  {
    if (!isnan(arr[i]))
    {
      result = false;
    }
  }
  return result;
}

bool isWithIn(double value, double min, double max){
  if (value >= min && value < max)
  {
    return true;
  }
  return false;
}