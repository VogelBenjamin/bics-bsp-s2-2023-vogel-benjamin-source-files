#include<stdlib.h>
#include<stdio.h>

int vectorSize = 3;

double* createVector();
void setVector(double* vector);
void printVector(double* vector);
double* vectorAddition(double *vector1, double *vector2);
double dotProduct(double *vector1, double *vector2);
double* crossProduct(double *vector1, double *vector2);

struct Particle {
  double* position;
  double* velocity;
  double* acceleration;
  double weight;
  double charge;
  double Hamaker_coef;
};

int main(){
  double **arrayVector;
  double *vector1, *vector2;
  int numberOfVectors;

  printf("Please input the number of vectors: ");
  scanf("%d",&numberOfVectors);

  arrayVector = (double**)malloc(numberOfVectors*sizeof(arrayVector));

  arrayVector[0] = vector1;
  arrayVector[1] = vector2;

  for (int i = 0; i < 2; ++i)
  {
    arrayVector[i] = createVector();
    setVector(arrayVector[i]);
    printVector(arrayVector[i]);
  }

  return 0;
}


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