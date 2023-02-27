#include<stdlib.h>
#include<stdio.h>
#include"linearAlgebra.h"

/*
struct particle {
  double* position;
  double* velocity;
  double* acceleration;
  double weight;
  double charge;
  double Hamaker_coef;
};

struct electricField {
  int energy;
};

struct magneticField {
  int energy;
};


void run(double **allParticles);
void move(struct particle particle);
*/

int main(){
  double **arrayVector;
  double *vector1, *vector2;
  int numberOfVectors;
  /*
  struct particle proton;
  proton.position = createVector();
  proton.position[0] = 0;
  proton.position[1] = 0;
  proton.position[2] = 0;
  proton.velocity = createVector();
  proton.velocity[0] = 0;
  proton.velocity[1] = 0;
  proton.velocity[2] = 0;

  */
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
/*
void move(struct particle particle){
  particle.velocity = vectorAddition(particle.velocity,particle.acceleration);
  particle.position = vectorAddition(particle.position,particle.velocity);
}
*/