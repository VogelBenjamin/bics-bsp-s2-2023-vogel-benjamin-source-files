#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"electromagnetism.h"
#include"linearAlgebra.h"
#include "omp.h"

int vecSize = 3;
const double K = 8.987551792E9;

/*
  allocate memory for a particle struct
*/
struct particle *createParticle(){
  struct particle *p = (struct particle*)malloc(sizeof(struct particle));
  p->position = (double*)malloc(vecSize*sizeof(p->position));
  p->velocity = (double*)malloc(vecSize*sizeof(p->velocity));
  p->acceleration = (double*)malloc(vecSize*sizeof(p->acceleration));
  p->attributes = (struct particleInfo*)malloc(sizeof(struct particleInfo));
  return p;
}
/*
  free memory allocated for a particle struct
*/
void freeParticle(struct particle *p){
  free(p->position);
  free(p->velocity);
  free(p->acceleration);
  free(p->attributes);
  free(p);
}

/*
  free memory allocated for a electricField struct
*/
void freeElField(struct electricField *eF){
  free(eF->fieldVector);
  free(eF);
}

/*
  free memory allocated for a magneticField struct
*/
void freeMagField(struct magneticField *mF){
  free(mF->fieldVector);
  free(mF);
}

/*
  free memory allocated for all particles in an array including the array itself
*/
void freeArrayOfParticles(struct particle **particleArray, int arraySize){
  for (int i = 0; i < arraySize; ++i)
  {
    freeParticle(particleArray[i]);
  }
  free(particleArray);
}

/*
  allocate memory for a electricField struct
*/
struct electricField *createElectricField(){
  struct electricField *eF = (struct electricField*)malloc(sizeof(eF));
  eF->fieldVector = (double*)malloc(sizeof(eF->fieldVector));
}

/*
  allocate memory for a magneticField struct
*/
struct magneticField *createMagneticField(){
  struct magneticField *mF = (struct magneticField*)malloc(sizeof(mF));
  mF->fieldVector = (double*)malloc(sizeof(mF->fieldVector));
}

/*
  reads file containg information on a particle and returns an array containing all important values
*/
double *getParticleAttributes(FILE *particleDoc){
  double *attr = (double*)malloc(4*sizeof(double));
  char output[15];
  char *ptr;

  // Extract mass from a file
  fgets(output, 15, particleDoc);
  fgets(output, 7, particleDoc);
  fgets(output, 5, particleDoc);
  attr[0] = strtod(output,&ptr);

  // Extract charge from file
  fgets(output, 9, particleDoc);
  fgets(output, 5, particleDoc);
  attr[1] = strtod(output,&ptr);

  // Extract radius from file
  fgets(output, 9, particleDoc);
  fgets(output, 5, particleDoc);
  attr[2] = strtod(output,&ptr);

  // Extract epsilon from file
  fgets(output, 10, particleDoc);
  fgets(output, 5, particleDoc);
  attr[3] = strtod(output,&ptr);

  return attr;
}

/*
  takes 4 arrays containing initial position, velcoity and acceleration and attributes and assigns them to a particle
*/
void setParticleInfo(struct particle *p, double *pos, double *vel, double *acc, double * attr){
  for (int i = 0; i < vecSize; ++i)
  {
    p->position[i] = pos[i];
    p->velocity[i] = vel[i];
    p->acceleration[i] = acc[i];
  }
  p->attributes->mass = attr[0];
  p->attributes->charge = attr[1];
  p->attributes->radius = attr[2];
  p->attributes->epsilon = attr[3];

}

/*
  sets the electric field to the given array 'norm'
*/
void setElectricField(struct electricField *eF, double *norm){
  for (int i = 0; i < vecSize; ++i)
  {
    eF->fieldVector[i] = norm[i];
  }
}

/*
  sets the magnetic field to the given array 'norm'
*/
void setMagneticField(struct magneticField *mF, double *norm){
  for (int i = 0; i < vecSize; ++i)
  {
    mF->fieldVector[i] = norm[i];
  }
}

/*
  print the memory addresses of fields of a particle struct
*/
void printDevInfo(struct particle *p){
  printf("Struct address:%p\n\n",p );
  printf("Dynamics: \n");
  printf("Position: %p\n",p->position);
  printf("Velocity: %p\n",p->velocity);
  printf("Acceleration: %p\n",p->acceleration);
  printf("\nAttributes:\n");
  printf("Struct: %p\n",p->attributes);
  printf("Mass: %p\n",&(p->attributes->mass));
  printf("Charge: %p\n",&(p->attributes->charge));
  printf("Radius: %p\n",&(p->attributes->radius));
}

/*
  print the fields of a particle struct
*/
void printAttributes(struct particle *p){
  printf("Position:");
  printVector(p->position);
  printf("Velocity:");
  printVector(p->velocity);
  printf("Acceleration:");
  printVector(p->acceleration);
  printf("\n");
  printf("Attributes: %lf, %lf, %lf, %lf\n\n",p->attributes->mass,p->attributes->charge,p->attributes->radius,p->attributes->epsilon);
}

/*
  calculates the change in position and velocity of a particle and modifies the particles fields accordingly
*/
void move(struct particle *particle,double timeStep){
  double velProg[3];
  double posProg[3];
  double posAcc[3];


  // calculation of the progression in velocity
  modifyVector(velProg, particle->acceleration);
  scalarMultiplication(timeStep,velProg);

  // calculation of the progression in position
  modifyVector(posProg, particle->velocity);
  modifyVector(posAcc, particle->acceleration);

  scalarMultiplication(timeStep,posProg);
  scalarMultiplication(0.5*pow(timeStep,2),posAcc);

  // adding the position and velocity progression to their respective fields
  vectorAddition(particle->position,posProg);
  vectorAddition(particle->position,posAcc);

  vectorAddition(particle->velocity,velProg);
}

/*
  applies the acceleration cause from electric and magnetic fields
*/
void applyAcceleration(struct particle *p, struct electricField *eF, struct magneticField *mF){
  double chargeMassRatio = p->attributes->charge/p->attributes->mass;
  double force;
  double crossP[3];
  crossProduct(p->velocity,mF->fieldVector,crossP);

  // for the 3 directions the acceleration force of both fields are calculated and added to the acceleration of the particle
  for (int i = 0; i < vecSize; ++i)
  {
    force = 0;
    force += chargeMassRatio*eF->fieldVector[i];
    force += chargeMassRatio*crossP[i];
    p->acceleration[i] = force;
  }
}

/*
  checks wheter two particles collide
*/
void checkCollision(struct particle *p1, struct particle *p2){
  double colRadius = p1->attributes->radius + p2->attributes->radius;
  double* direct = connectPoints(p1->position,p2->position);
  double distance = norm(direct);
  if (distance > colRadius)
  {

  }
  else{

  }
  free(direct);
}

/*
  calculates the forces caused by a collision and applies them on the respective particles
*/
void calcCollision(struct particle *p1, struct particle *p2){

  double relativeSpeed[3];
  double normalSpeed[3];
  double *direct;

  // vector between the two centers of the particles
  direct = connectPoints(p1->position,p2->position);

  normalize(direct);

  // calculate the relative velocity
  modifyVector(relativeSpeed,p1->velocity);
  vectorSubtraction(relativeSpeed,p2->velocity);

  // calculate the velocity in the direction of 'direct'
  modifyVector(normalSpeed,direct);
  scalarMultiplication(dotProduct(relativeSpeed,direct),normalSpeed);

  // the respective normal (direct) forces are applied on the other particle
  vectorAddition(p2->velocity,normalSpeed);
  scalarMultiplication(-1,normalSpeed);
  vectorAddition(p1->velocity,normalSpeed);

  free(direct);
}
/*
  approximates the van der waals force between two particles
*/
double *lennardJonesPotentialForce(struct particle *p1, struct particle *p2){
  double *direct;
  double distance;
  double sigma;
  double epsilon;
  double lJForce;

  // get distance between the two particles
  direct = connectPoints(p1->position,p2->position);
  distance = norm(direct);

  // determine sigma
  sigma = 0.5*(p1->attributes->radius + p2->attributes->radius);

  // determine depth of potential well
  epsilon = pow((p1->attributes->epsilon * p1->attributes->epsilon),0.5);

  // calculate potential
  lJForce = 48*epsilon*(pow(sigma,6)/pow(distance,6))*((pow(sigma,6)/pow(distance,7))-0.5*(pow(distance,-1)));

  // apply the intensity of the force to the direction unit vector to obtain the acctual force vector
  scalarMultiplication(lJForce,direct);

  return direct;
}

double *electroStaticForce(struct particle *p1, struct particle *p2){
  double *direct;
  double eSForce;

  // determine direction of the force
  direct = connectPoints(p1->position,p2->position);

  // calculates the intensity of the Force
  eSForce = K*(p1->attributes->charge)*(p2->attributes->charge)/pow(norm(direct),2);

  // normalizes the direction vector and multiplies magnitude of Force to it
  normalize(direct);

  scalarMultiplication(eSForce,direct);

  return direct;
}