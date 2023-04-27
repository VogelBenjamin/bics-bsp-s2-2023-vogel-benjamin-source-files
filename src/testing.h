#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

bool compareArray(double *array1, double *array2, int size);
bool compareValue(double v1, double v2, bool *check);
bool checkIfAllNAN(double *arr, int size);
bool isWithIn(double value, double min, double max);