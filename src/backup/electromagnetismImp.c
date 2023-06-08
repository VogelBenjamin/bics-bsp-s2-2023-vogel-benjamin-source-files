#include<stdlib.h>
#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#include"electromagnetism.h"
#include"linearAlgebra.h"
#include "testing.h"
#include "omp.h"

int vecSize = 3;
double timeStep = 0.000001;
double electroKonst = 8.987551792E9;
double chargeToCoulombsRatio = 1.6022E-19;
double uToKgRatio = 1.66054E-27;
double avogadroConst = 6.02214076E23;
double scaling = 1E-7;

double getStep(){
  return timeStep;
}

void setStep(double val){
  timeStep = val;
}

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
  eF->fieldVector = (double*)malloc(sizeof(eF->fieldVector)*vecSize);
}

/*
  allocate memory for a magneticField struct
*/
struct magneticField *createMagneticField(){
  struct magneticField *mF = (struct magneticField*)malloc(sizeof(mF));
  mF->fieldVector = (double*)malloc(sizeof(mF->fieldVector)*vecSize);
}

/*
  reads file containg information on a particle and returns an array containing all important values
*/
double *getParticleAttributes(FILE *particleDoc){
  double *attr = (double*)malloc(5*sizeof(double));
  char output[15];
  char *ptr;

  // Extract mass from a file
  fgets(output, 15, particleDoc);
  fgets(output, 7, particleDoc);
  fgets(output, 7, particleDoc);
  attr[0] = strtod(output,&ptr);

  // Extract charge from file
  fgets(output, 9, particleDoc);
  fgets(output, 7, particleDoc);
  attr[1] = strtod(output,&ptr);

  // Extract radius from file
  fgets(output, 9, particleDoc);
  fgets(output, 7, particleDoc);
  attr[2] = strtod(output,&ptr);

  // Extract epsilon from file
  fgets(output, 10, particleDoc);
  fgets(output, 7, particleDoc);
  attr[3] = strtod(output,&ptr);

  // Extract the sigma from the file
  fgets(output, 8, particleDoc);
  fgets(output, 7, particleDoc);
  attr[4] = strtod(output,&ptr);

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
  p->attributes->sigma = attr[4];

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
  print position into a given file for csv usage
*/
void printPos(struct particle *p, FILE *particleDoc, int ID){
  char buffer[15];
  for (int i = 0; i < 2; ++i)
  {
    snprintf(buffer, 15, "%lf, ", p->position[i]);
    fputs(buffer,particleDoc);
  }
  snprintf(buffer, 15, "%lf, %d\n", p->position[2], ID);
  fputs(buffer,particleDoc);
}

/*
  calculates the change in position and velocity of a particle and modifies the particles fields accordingly
*/
void move(struct particle *particle,double timeStep){
  double velProg[3];
  double posProg[3];
  double posAcc[3];


  // calculation of the change in velocity
  modifyVector(velProg, particle->acceleration);
  scalarMultiplication(timeStep,velProg);

  // calculation of the change in position

  // vΔt
  modifyVector(posProg, particle->velocity);
  modifyVector(posAcc, particle->acceleration);

  // (aΔt^2)/2
  scalarMultiplication(timeStep,posProg);
  scalarMultiplication(0.5*pow(timeStep,2),posAcc);

  // adding the position and velocity variation to their respective fields
  vectorAddition(particle->position,posProg);
  vectorAddition(particle->position,posAcc);

  vectorAddition(particle->velocity,velProg);
}

/*
  applies the acceleration cause from electric and magnetic fields
*/
void applyFieldAcceleration(struct particle *p, double *eF, double *mF){
  double massRatio = 1/uToKg(p->attributes->mass);
  double *force = (double*)calloc(vecSize,sizeof(force));
  double *eFo, *mFo;

  // for the 3 directions the acceleration force of both fields are calculated and added to the acceleration of the particle

  eFo = electricForce(p, eF);
  mFo = magneticForce(p, mF);

  vectorAddition(force, eFo);
  vectorAddition(force, mFo);

  scalarMultiplication(massRatio, force);
  scalarMultiplication(scaling, force);

  vectorAddition(p->acceleration, force);

  free(eFo);
  free(mFo);
  free(force);
}

void applyInterAcceleration(struct particle *p1, struct particle *p2){
  double *force = (double*)calloc(vecSize,sizeof(force));
  double *temp = (double*)calloc(vecSize,sizeof(force));
  double *ljFo, *elstFo;

  // apply electrostatic force
  if (p1->attributes->charge != 0 && p2->attributes->charge != 0)
  {
    elstFo = electroStaticForce(p1,p2);
    modifyVector(temp, elstFo);


    scalarMultiplication(1/uToKg(p2->attributes->mass),elstFo);
    scalarMultiplication(-1/uToKg(p1->attributes->mass), temp);
    #pragma omp critical(fazz)
      // p2
      vectorAddition(p2->acceleration, elstFo);
      // p1
      vectorAddition(p1->acceleration, temp);

    free(elstFo);
  }
  else {
    // apply lennard-jones force
    ljFo = lennardJonesPotentialForce(p1, p2);
    modifyVector(temp, ljFo);

    // p2
    scalarMultiplication(1/uToKg(p2->attributes->mass), ljFo);

    vectorAtomicAddition(p2->acceleration, ljFo);

    // p1
    scalarMultiplication(-1/uToKg(p1->attributes->mass), temp);

    vectorAtomicAddition(p1->acceleration, temp);
    free(ljFo);
  }

  free(temp);
  free(force);
}

double *electricForce(struct particle *p, double *eF){
  double charge;
  double *elForce;

  // initialization of values for the final computation
  charge = chargeToCoulombs(p->attributes->charge);
  elForce = (double*)calloc(vecSize,sizeof(elForce));

  // qE
  modifyVector(elForce, eF);
  scalarMultiplication(charge,elForce);

  return elForce;
}

double *magneticForce(struct particle *p, double *mF){
  double charge;
  double *magForce;
  double crossP[3];

  // initialization of values for the final computation
  charge = chargeToCoulombs(p->attributes->charge);
  magForce = (double*)calloc(vecSize,sizeof(magForce));

  // v x B
  crossProduct(p->velocity,mF,crossP);

  // q(v x B)
  modifyVector(magForce,crossP);
  scalarMultiplication(charge,magForce);

  return magForce;
}

/*
  checks a priori wheter two particles collide
*/
double checkPrioriCollision(struct particle *p1, struct particle *p2){
  double colRadius = 1E-12*(p1->attributes->radius + p2->attributes->radius);

  if (nanoToMeter(distance(p1->position,p2->position)) < colRadius)
  {
    return 0;
  }

  double **augementedMatrix = (double**)malloc(sizeof(double*)*3);
  double *c1 = (double*)malloc(sizeof(double)*3);
  double *c2 = (double*)malloc(sizeof(double)*3);
  double *c3 = (double*)malloc(sizeof(double)*3);
  double *solution;
  double ratio;

  // create augemented matrix corresponding to linear system composed of the instantaneous tragectory of two particles
  modifyVector(c1, p1->velocity);
  scalarMultiplication(timeStep, c1);

  modifyVector(c2, p2->velocity);
  scalarMultiplication(-timeStep, c2);

  modifyVector(c3, p2->position);
  vectorSubtraction(c3, p1->position);

  augementedMatrix[0] = c1;
  augementedMatrix[1] = c2;
  augementedMatrix[2] = c3;

  // solve the system to determine, whether the trajectory intersect

  solution = solveAugMatrix(augementedMatrix,3,2);

  if (solution[0] < 0 || solution[1] < 0 || (solution[0] == 0 && solution[1] == 0 && !compareArray(p1->position,p2->position,3)))
  {
    ratio = NAN;
    freeMatrix(augementedMatrix,3);

    return ratio;
  }

  if (solution[1] == 0)
  {
    solution[0] = augementedMatrix[2][0]/(augementedMatrix[0][0]+augementedMatrix[1][0]);
    solution[1] = augementedMatrix[2][0]/(augementedMatrix[0][0]+augementedMatrix[1][0]);
    ratio = evaluateCollision(solution);
    freeMatrix(augementedMatrix,3);
    free(solution);
    return ratio;
  }

  ratio = evaluateCollision(solution);
  freeMatrix(augementedMatrix,3);
  return ratio;
}

double evaluateCollision(double *sol){
  double distance = 0;
  if (!checkIfAllNAN(sol,2) && isWithIn(sol[0],0,1) && isWithIn(sol[1],0,1) && abs(sol[0]-sol[1] < 0.001))
  {
    distance = sol[0];
    free(sol);
    return distance;
  }
  free(sol);
  return NAN;
}

/*
  calculates the forces caused by a collision and applies them on the respective particles
*/
void calcCollision(struct particle *p1, struct particle *p2, double *direct){

  double vCenter1[3];
  double vCenter2[3];
  double vNormal1[3];
  double vNormal2[3];

  if (abs(norm(p1->velocity) < 1E-10) || abs(norm(p2->velocity)) < 1E-10)
  {
    double relativeSpeed[3];
    double normalSpeed[3];

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
  } else{
    double angle1 = angleBetween(direct, p1->velocity);
    scalarMultiplication(-1,direct);
    double angle2 = angleBetween(direct, p2->velocity);

    modifyVector(vCenter1,p1->velocity);
    scalarMultiplication(cos(angle1),vCenter1);

    modifyVector(vCenter2,p2->velocity);
    scalarMultiplication(cos(angle2),vCenter2);

    modifyVector(vNormal1, p1->velocity);
    vectorSubtraction(vNormal1, vCenter1);

    modifyVector(vNormal2, p2->velocity);
    vectorSubtraction(vNormal2, vCenter2);

    modifyVector(p1->velocity, vCenter2);
    modifyVector(p2->velocity, vCenter1);

    vectorSubtraction(p1->velocity, vNormal1);
    vectorSubtraction(p2->velocity, vNormal2);
  }

  free(direct);

}

/*
  executes necessary functions for the checking of collisions
*/
bool handleCollision(struct particle *p1, struct particle *p2, int pID, int *handleArray){
  double eval = checkPrioriCollision(p1,p2);
  double *direct;
  double a[3],b[3],c[3],d[3];
  int check = 1;
  //#pragma omp critical(luzz)
    if (!(handleArray[pID] == 0 && !isnan(eval)))
    {
      check = 0;
    }

  if (check == 0)
  {
    return check;
  }
  if (eval)
  {
    #pragma omp critical(lizz)
      move(p1,eval*timeStep);
      move(p2,eval*timeStep);
  }

  modifyVector(a,p1->position);
  modifyVector(b,p2->position);

  modifyVector(c,p1->velocity);
  modifyVector(d,p2->velocity);

  scalarMultiplication(0.001,c);
  scalarMultiplication(0.001,d);

  vectorSubtraction(a,c);
  vectorSubtraction(b,d);
  // vector between the two centers of the particles
  direct = connectPoints(a,b);
  normalize(direct);

  #pragma omp critical(lizz)
    calcCollision(p1, p2, direct);
    move(p1,(1-eval)*timeStep);
    move(p2,(1-eval)*timeStep);

  return check;
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

  // initialization of values for the final computation
  sigma = p1->attributes->sigma;
  epsilon = p1->attributes->epsilon;

  // get distance between the two particles and the direction for the force
  direct = connectPoints(p1->position,p2->position);
  distance = norm(direct);
  distance = meterToAngstrom(distance);
  normalize(direct);

  // calculate potential
  lJForce = 48*epsilon*(pow(sigma,6)/pow(distance,6))*((pow(sigma,6)/pow(distance,7))-0.5*(pow(distance,-1)));

  // apply the intensity of the force to the direction unit vector to obtain the acctual force vector
  scalarMultiplication(lJForce,direct);

  return direct;
}

double *electroStaticForce(struct particle *p1, struct particle *p2){
  double *direct;
  double distance;
  double eSForce;

  // get distance between the two particles and the direction for the force
  direct = connectPoints(p1->position,p2->position);
  distance = norm(direct);
  normalize(direct);

  // force magnitude
  eSForce = electroKonst*(chargeToCoulombs(p1->attributes->charge))*(chargeToCoulombs(p2->attributes->charge))/pow(distance,2);

  // final vector
  scalarMultiplication(eSForce,direct);

  return direct;
}

double chargeToCoulombs(double charge){
  return chargeToCoulombsRatio*charge;
}

double nanoToMeter(double value){
  return value*1E-9;
}

double uToKg(double weight){
  return uToKgRatio*weight;
}

double meterToAngstrom(double distance){
  return 1E10*distance;
}

