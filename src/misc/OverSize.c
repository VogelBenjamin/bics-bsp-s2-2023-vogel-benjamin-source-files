#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<stdbool.h>
#include<string.h>
#include<time.h>
#include "omp.h"

bool compareArray(double *array1, double *array2, int size);
bool compareValue(double v1, double v2, bool *check);
bool checkIfAllNAN(double *arr, int size);
bool isWithIn(double value, double min, double max);

void setVecSize(int val);
double* createVector();
void freeMatrix(double **matrix, int col);
double* connectPoints(double* pos1, double* pos2);
void setVector(double *vector);
void modifyVector(double *oldVector, double *newVector);
void printVector(double *vector);
double norm(double* vector);
double distance(double* pos1, double* pos2);
void vectorAddition(double *vector1, double *vector2);
void vectorAtomicAddition(double *vector1, double *vector2);
void vectorSubtraction(double *vector1, double *vector2);
void scalarMultiplication(double scalar, double *vector);
double dotProduct(double *vector1, double *vector2);
void crossProduct(double *vector1, double *vector2,double* cP);
void normalize(double *vector);
double angleBetween(double *v1, double * v2);
double solve(double *equation, int eqSize);
double *solveAugMatrix(double **augMatrix, int numEq, int numVar);

struct particleInfo {
  double mass;
  double charge;
  double radius;
  double epsilon;
  double sigma;
};

struct particle {
  double *position;
  double *velocity;
  double *acceleration;
  struct particleInfo *attributes;
};

struct electricField {
  double *fieldVector;
};

struct magneticField {
  double *fieldVector;
};

double getStep();
void setStep(double val);
struct particle *createParticle();
void freeParticle(struct particle *p);
void freeElField(struct electricField *eF);
void freeMagField(struct magneticField *mF);
void freeArrayOfParticles(struct particle **particleArray, int arraySize);
struct electricField *createElectricField();
struct magneticField *createMagneticField();
double *getParticleAttributes(FILE *particleDoc);
void setParticleInfo(struct particle *p, double *pos, double *vel, double *acc, double *attr);
void setElectricField(struct electricField *eF, double *norm);
void setMagneticField(struct magneticField *mF, double *norm);
void printDevInfo(struct particle *p);
void printAttributes(struct particle *p);
void printPos(struct particle *p, FILE *particleDoc, int ID);
void move(struct particle *particle,double timeStep);
void applyFieldAcceleration(struct particle *p, double *eF, double *mF);
void applyInterAcceleration(struct particle *p1, struct particle *p2);
double *electricForce(struct particle *p, double *eF);
double *magneticForce(struct particle *p, double *mF);
double checkPrioriCollision(struct particle *p1, struct particle *p2);
double evaluateCollision(double *sol);
void calcCollision(struct particle *p1, struct particle *p2, double *direct);
bool handleCollision(struct particle *p1, struct particle *p2, int pID, int *handleArray);
double *lennardJonesPotentialForce(struct particle *p1, struct particle *p2);
double *electroStaticForce(struct particle *p1, struct particle *p2);
double chargeToCoulombs(double charge);
double nanoToMeter(double value);
double uToKg(double weight);
double meterToAngstrom(double distance);


void printDirectory();
int **createPartition(int size);

double currTimeStep;

int main(int argc, char *argv[]){

  char *ptr;

  int menuInput = 0;
  int particleArraySize;
  int i = 0;
  int j = 0;
  int k;
  int timePeriod;
  double deltaX;
  double deltaY;
  double deltaZ;


  struct particle **particleArray;
  struct particle *tempParticle;
  struct electricField *elFi = createElectricField();
  struct magneticField *magFi = createMagneticField();

  currTimeStep = getStep();

  double *posVec = createVector();
  double *velVec = createVector();
  double *accVec = createVector();
  double attr[] = {1,1,1,1,1};
  double elNorm[3];
  double magNorm[3];
  double initialVec[3] = {0,0,0};
  double spawnSphere[3];
  double spawnPoint[3];


  particleArraySize = 2;
  timePeriod = 1;

  velVec[0] = 1;
  velVec[0] = 0;
  velVec[0] = 1;

  elNorm[0] = 0;
  elNorm[1] = 1;
  elNorm[2] = 0;
  for (int i = 0; i < 3; ++i)
  {
    spawnSphere[i] = 0.5;
  }
  magNorm[1] = 1;
  magNorm[1] = 0;
  magNorm[1] = 1;
  modifyVector(spawnPoint,spawnSphere);

  particleArray = (struct particle**)malloc(particleArraySize*sizeof(struct particle*));

  modifyVector(posVec, initialVec);
  modifyVector(accVec, initialVec);

  setElectricField(elFi, elNorm);
  setMagneticField(magFi,magNorm);

  srand(time(0));

  for (i = 0; i < particleArraySize; ++i)
  {
    tempParticle = createParticle();


    deltaX = rand();
    deltaY = rand();
    deltaZ = rand();

    spawnPoint[0] *= deltaX;
    spawnPoint[1] *= deltaY;
    spawnPoint[2] *= deltaZ;

    scalarMultiplication(5E-10,spawnPoint);

    setParticleInfo(tempParticle, spawnPoint, velVec, accVec, attr);

    particleArray[i] = tempParticle;
    modifyVector(spawnPoint,spawnSphere);
  }

  int *handled = (int*)malloc(particleArraySize*sizeof(handled));
  double *elF = elFi->fieldVector;
  double *magF = elFi->fieldVector;
  int **particleIdxPartition = createPartition(particleArraySize);
  int chunk = ((particleArraySize*(particleArraySize-1))/2);
  /*
  for (int i = 0; i < particleArraySize; ++i)
  {
    handled[i] = 0;
  }
  */
  for (i = 0; i < timePeriod*(1/currTimeStep); ++i)
  {
    /*
    if (i % 100000 == 0 && argc != 13)
    {
      printf("\n%lfs:\n",i*currTimeStep);
    }
  */
    // apply the acceleration to all particles
    #pragma omp parallel shared(particleArraySize,particleArray,elF,magF, particleIdxPartition, chunk)
      #pragma omp for schedule(static,1) private(j,k)
        for (j = 0; j < particleArraySize; ++j)
        {

          applyFieldAcceleration(particleArray[j],elF,magF);
          int id = omp_get_thread_num();
          int total = omp_get_num_threads();
          int *ptr;
          //printf("id: %d, total: %d, chunk: %d, part: %d\n", id, total, chunk,chunk/total );
          for (k = chunk*id/total; k < chunk*(id+1)/total; ++k)
          {
            ptr = particleIdxPartition[k];
            if (ptr)
            {
              applyInterAcceleration(particleArray[ptr[0]],particleArray[ptr[1]]);
            }
          }


          for (k = chunk*id/total; k < chunk*(id+1)/total; ++k)
          {
            ptr = particleIdxPartition[k];
            if (ptr && handled[j] == 0 && handled[k] == 0 && handleCollision(particleArray[ptr[0]],particleArray[ptr[1]],particleIdxPartition[k][0],handled))
            {
              handled[j] = 1;
              handled[k] = 1;
            }
          }
          // check whether the movement of a particle has already been computed
          if (handled[j] == 0)
          {
            move(particleArray[j],currTimeStep);
          }

          handled[j] = 0;
        }

    if (i % 100000 == 0 && argc != 13)
    {
      for (int j = 0; j < particleArraySize; ++j)
      {
        printAttributes(particleArray[j]);

      }
    }


    for (int j = 0; j < particleArraySize; ++j)
    {
      handled[j] = 0;
      scalarMultiplication(0,particleArray[j]->acceleration);
    }

  }

  freeArrayOfParticles(particleArray,particleArraySize);
  free(handled);
  free(posVec);
  free(velVec);
  free(accVec);
  for (int i = 0; i < chunk; ++i)
  {
    free(particleIdxPartition[i]);
  }
  free(particleIdxPartition);

  freeElField(elFi);
  freeMagField(magFi);

  return 0;
}



int **createPartition(int size){
  int partitionSize = (size*(size-1))/2;
  partitionSize += 8 - (partitionSize % omp_get_max_threads());
  int counter = 0;
  int **partition = (int**)malloc(sizeof(int*)*partitionSize);
  int *pair;
  for (int i = 0; i < size; ++i)
  {
    for (int j = i+1; j < size; ++j)
    {
      pair = (int*)malloc(sizeof(int)*2);
      pair[0] = i;
      pair[1] = j;
      partition[counter] = pair;
      counter++;
    }
  }

  for (int i = counter; i < partitionSize; ++i)
  {
    partition[i] = NULL;
  }
  return partition;
}

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
  norm = sqrt(pow(vector[0],2) + pow(vector[1],2) + pow(vector[2],2));
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

  // transform matrix into reduced row echelon form
  // int i decides which column to operate on
  for (int i = 0; i < numVar; ++i)
  {
    // makes sure that augementmatrix[i][i] is not 0
    if (augMatrix[i][i] == 0)
    {

      for (int j = i+1; j < numEq; ++j)
      {

        if (augMatrix[i][j] != 0)
        {

          for (int k = 0; k < numVar+1; ++k)
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
      for (int j = 0; j < numVar+1; ++j)
      {
        augMatrix[j][i] *= ratio;
      }

      // makes sure the values above and below the leading are turned to 0
      // j indicates which row is operated on
      for (int j = 0; j < numEq; ++j)
      {
        if (j != i)
        {
          sub = augMatrix[i][j];
          // k indicates on which element of the row is operated on
          for (int k = 0; k < numEq; ++k)
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

  // Extract the sigma from the file
  fgets(output, 8, particleDoc);
  fgets(output, 5, particleDoc);
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

    // p2
    scalarMultiplication(1/uToKg(p2->attributes->mass),elstFo);
    #pragma omp critical(fazz)
      vectorAddition(p2->acceleration, elstFo);

    // p1
    scalarMultiplication(-1/uToKg(p1->attributes->mass), temp);
    #pragma omp critical(fazz)
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
    free(solution);

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
  free(solution);
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

  #pragma omp critical(bazz)
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





