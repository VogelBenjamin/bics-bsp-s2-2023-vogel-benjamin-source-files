#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <stdbool.h>
#include "linearAlgebra.h"
#include "electromagnetism.h"
#include "testing.h"
#include "omp.h"



bool testNorm(){
  printf("NORM TESTING\n");
  double v1[] = {1,0,0};
  double v2[] = {0,1,0};
  double v3[] = {0,0,1};
  double v4[] = {-1,0,0};
  double v5[] = {1,2,2};
  double v6[] = {-2,4,-4};
  double v7[] = {0,0,0};

  bool success = true;
  bool *successPtr = &success;

  compareValue(norm(v1), 1, successPtr);

  compareValue(norm(v2), 1, successPtr);

  compareValue(norm(v3), 1, successPtr);

  compareValue(norm(v4), 1, successPtr);

  compareValue(norm(v5), 3, successPtr);

  compareValue(norm(v6), 6, successPtr);

  compareValue(norm(v7), 0, successPtr);

  return success;
}

bool testConnectPoints(){
  printf("\nPOINTCONNECTION TESTING\n");

  bool success = true;

  // TEST 1

  double p1[] = {0,0,0};
  double p2[] = {2,3,-4};
  double *v3 = connectPoints(p1,p2);
  double *v4 = connectPoints(p2,p1);

  double r3[] = {2,3,-4};
  double r4[] = {-2,-3,4};

  bool result = false;

  if (!(compareArray(v3,r3,3) && compareArray(v4,r4,3)))
  {
    success = false;
  }

  // TEST 2

  double p5[] = {5,-3,7};
  double p6[] = {12,4,-4};
  double *v7 = connectPoints(p5,p6);
  double *v8 = connectPoints(p6,p5);

  double r7[] = {7,7,-11};
  double r8[] = {-7,-7,11};

  result = false;

  if (!(compareArray(v7,r7,3) && compareArray(v8,r8,3)))
  {
    success = false;
  }


  free(v3);
  free(v4);
  free(v7);
  free(v8);

  return success;
}

bool testNormalize(){
  printf("\nNORMALIZE TESTING\n");

  bool success = true;

  double v1[] = {1,1,1};
  double v2[] = {2,-5,9};

  normalize(v1);
  normalize(v2);

  bool result = false;

  if (!(norm(v1) == 1.0 && norm(v2) == 1.0))
  {
    success = false;
  }

  return success;
}

bool testDotProduct(){
  printf("\nDOT PRODUCT TESTING\n");

  bool success = true;

  double v1[] = {1,2,3};
  double v2[] = {-3,2,7};

  double v3[] = {3,6,9};
  double v4[] = {10,-2,-2};

  double v5[] = {0,0,0};

  bool result = false;

  if (!(dotProduct(v1,v2) == 22 && dotProduct(v3,v4) == 0 && dotProduct(v1,v3) == dotProduct(v3,v1) && dotProduct(v5,v2) == 0))
  {
    success = false;
  }

  return success;

}

bool testSolveAugMatrix(){
  printf("\nLINEAR SYSTEM SOLVING\n");

  bool success = true;

  double **augMatrix= (double**)malloc(sizeof(double*)*3);
  double *answer;

  // TEST 1
  printf("\nTEST 1 \n");
  double m[3][3] = {{1,1,1},{2,-3,-5},{10,-5,-11}};

  for (int i = 0; i < 3; ++i)
  {
    double *row = (double*)malloc(sizeof(double)*3);
    modifyVector(row,m[i]);
    augMatrix[i] = row;
  }


  answer = solveAugMatrix(augMatrix, 3,2);
  double r1[] = {4,3};

  if (!compareArray(answer,r1,2)){
    printf("Test1: Error\n");
    success = false;
  }

  free(answer);
  for (int i = 0; i < 3; ++i)
  {
    free(augMatrix[i]);
  }

  // TEST 2
  printf("\nTEST 2 \n");
  double m2[3][3] = {{1,1,1},{2,-3,-5},{10,-5,-14}};

  for (int i = 0; i < 3; ++i)
  {
    double *row = (double*)malloc(sizeof(double)*3);
    modifyVector(row,m2[i]);
    augMatrix[i] = row;
  }


  answer = solveAugMatrix(augMatrix, 3, 2);


  if (!checkIfAllNAN(answer,2)){
    printf("Test2: Error\n");
    success = false;
  }


  for (int i = 0; i < 3; ++i)
  {
    free(augMatrix[i]);
  }
  free(answer);

  // TEST 3
  printf("\nTEST 3 \n");
  double m3[3][3] = {{-2,2,0},{2,2,0},{0,1,0}};

  for (int i = 0; i < 3; ++i)
  {
    double *row = (double*)malloc(sizeof(double)*3);
    modifyVector(row,m3[i]);
    augMatrix[i] = row;
  }


  answer = solveAugMatrix(augMatrix, 3, 2);


  double r2[] = {0.25,0.25};

  if (!compareArray(answer,r2,2)){
    printf("Test3: Error\n");
    success = false;
  }


  for (int i = 0; i < 3; ++i)
  {
    free(augMatrix[i]);
  }

  free(answer);

  // TEST 4
  printf("\nTEST 4 \n");
  double m4[3][3] = {{0,2,6},{1,1,3},{1,5,15}};

  for (int i = 0; i < 3; ++i)
  {
    double *row = (double*)malloc(sizeof(double)*3);
    modifyVector(row,m4[i]);
    augMatrix[i] = row;
  }


  answer = solveAugMatrix(augMatrix, 3, 2);


  double r3[] = {2,1};

  if (!compareArray(answer,r3,2)){
    printf("Test4: Error\n");
    success = false;
  }


  for (int i = 0; i < 3; ++i)
  {
    free(augMatrix[i]);
  }

  free(augMatrix);
  free(answer);


  return success;

}

bool testCollisionCheck(){
  printf("\nCOLLISON CHECK TESTING\n");

  // TEST 1
  printf("\nTEST 1 \n");
  bool success = true;

  struct particle *p1 = createParticle();
  struct particle *p2 = createParticle();

  double pos1[] = {-1,2,0};
  double vel1[] = {2,-2,0};
  double pos2[] = {-1,0,0};
  double vel2[] = {2,2,0};
  double acc[] = {0,0,0};
  double att[] = {1,0,1};

  double res = checkPrioriCollision(p1, p2);

  setParticleInfo(p1,pos1,vel1,acc,att);
  setParticleInfo(p2,pos2,vel2,acc,att);

  if(abs(res-0.5) > 0.001){
    printf("Test1: Error\n");
    success = false;
  }


  // TEST 2
  printf("\nTEST 2 \n");
  double pos3[] = {1,1,1};
  double vel3[] = {2,-2,0};
  double pos4[] = {3,5,4};
  double vel4[] = {2,0,3};

  setParticleInfo(p1,pos3,vel3,acc,att);
  setParticleInfo(p2,pos4,vel4,acc,att);

  res = checkPrioriCollision(p1, p2);

  if(!isnan(res)){
    printf("Test2: Error\n");
    success = false;
  }


  // TEST 3
  printf("\nTEST 3 \n");
  double pos5[] = {1,1,1};
  double vel5[] = {2,2,2};
  double pos6[] = {0,0,0};
  double vel6[] = {1,1,1};

  setParticleInfo(p1,pos5,vel5,acc,att);
  setParticleInfo(p2,pos6,vel6,acc,att);

  res = checkPrioriCollision(p1, p2);
  if(!isnan(res)){
    printf("Test3: Error\n");
    success = false;
  }

  freeParticle(p1);
  freeParticle(p2);

  return success;
}

bool testFullCollisionScenario(){
  printf("\nPROPPER COLLISON TESTING\n");

  // TEST 1


  bool success = true;

  setStep(1);

  printf("\nTEST 1 \n");

  struct particle *p1 = createParticle();
  struct particle *p2 = createParticle();

  double pos1[] = {0,0,0};
  double vel1[] = {3,0,0};
  double pos2[] = {2,0,0};
  double vel2[] = {0,0,0};
  double acc[] = {0,0,0};
  double att[] = {1,0,1,0,0};

  setParticleInfo(p1,pos1,vel1,acc,att);
  setParticleInfo(p2,pos2,vel2,acc,att);

  handleCollision(p1, p2);

  double res1[] = {2,0,0};
  double res2[] = {3,0,0};

  if (!(compareArray(p1->position,res1,3) && compareArray(p1->velocity,vel2,3) && compareArray(p2->position,res2,3) && compareArray(p2->velocity,vel1,3)))
  {
    printf("Test1: Error\n");
    success = false;
  }

  // TEST 2
  printf("\nTEST 2 \n");

  double pos3[] = {0,0,0};
  double vel3[] = {4,8,0};
  double pos4[] = {2,4,0};
  double vel4[] = {-1,-2,0};

  setParticleInfo(p1,pos3,vel3,acc,att);
  setParticleInfo(p2,pos4,vel4,acc,att);

  handleCollision(p1, p2);

  double res3[] = {1,2,0};
  double res4[] = {4.89,9.78,0};

  if (!(compareArray(p1->position,res3,3) && compareArray(p1->velocity,vel4,3) && compareArray(p2->position,res4,3) && compareArray(p2->velocity,vel3,3)))
  {
    printf("Test2: Error\n");
    success = false;
  }

  // TEST 3

  printf("\nTEST 3 \n");

  double pos5[] = {0,0,0};
  double vel5[] = {2,3,4};
  double pos6[] = {7,8,9};
  double vel6[] = {-4,5,3};

  setParticleInfo(p1,pos5,vel5,acc,att);
  setParticleInfo(p2,pos6,vel6,acc,att);

  if(!handleCollision(p1, p2)){
    move(p1,1);
    move(p2,1);
  }

  double res5[] = {2,3,4};
  double res6[] = {3,13,12};

  if (!(compareArray(p1->position,res5,3) && compareArray(p1->velocity,vel5,3) && compareArray(p2->position,res6,3) && compareArray(p2->velocity,vel6,3)))
  {
    printf("Test3: Error\n");
    success = false;
  }

  // TEST 4

  printf("\nTEST 4 \n");

  double pos7[] = {9,12,0};
  double vel7[] = {6,-5.5,0};
  double pos8[] = {20,10,0};
  double vel8[] = {-10,-1.5,0};

  setParticleInfo(p1,pos7,vel7,acc,att);
  setParticleInfo(p2,pos8,vel8,acc,att);

  handleCollision(p1, p2);

  double res7[] = {7.231094,9.09,0};
  double res8[] = {15.924,7.652,0};

  if (!(compareArray(p1->position,res7,3) && compareArray(p2->position,res8,3)))
  {
    printf("Test4: Error\n");
    success = false;
  }

  // TEST 5

  printf("\nTEST 5 \n");

  setStep(0.001);

  double pos9[] = {0,0,0};
  double vel9[] = {1,1,1};
  double pos10[] = {1,1,1};
  double vel10[] = {-1,-1,-1};

  setParticleInfo(p1,pos9,vel9,acc,att);
  setParticleInfo(p2,pos10,vel10,acc,att);

  handleCollision(p1, p2);

  double res9[] = {0,0,0};
  double res10[] = {1,1,1};

  if (!(compareArray(p1->position,res9,3) && compareArray(p2->position,res10,3)))
  {
    printf("Test5: Error\n");
    success = false;
  }

  // TEST 6

  printf("\nTEST 6 \n");

  double pos11[] = {0.048,0.032,0.16};
  double vel11[] = {18,-27,21};
  double pos12[] = {2,2,2};
  double vel12[] = {-5838,-5931,-5499};

  setParticleInfo(p1,pos11,vel11,acc,att);
  setParticleInfo(p2,pos12,vel12,acc,att);

  handleCollision(p1, p2);

  double res11[] = {-3.8480,-3.9158,-3.5107};
  double res12[] = {0.0559,0.0201,3.4107};

  if (!(compareArray(p1->position,res11,3) && compareArray(p2->position,res12,3)))
  {
    printf("Test6: Error\n");
    success = false;
  }

  // TEST 7

  printf("\nTEST 7 \n");

  setStep(0.1);

  double pos13[] = {1,0,0};
  double vel13[] = {12,0,0};
  double pos14[] = {8,0,0};
  double vel14[] = {-9,0,0};

  setParticleInfo(p1,pos13,vel13,acc,att);
  setParticleInfo(p2,pos14,vel14,acc,att);

  for (int i = 0; i < 10; ++i)
  {
    if (!handleCollision(p1,p2))
    {
      move(p1,0.1);
      move(p2,0.1);
    }
  }


  double res13[] = {-1,0,0};
  double res14[] = {13,0,0};

  if (!(compareArray(p1->position,res13,3) && compareArray(p2->position,res14,3)))
  {
    printf("Test7: Error\n");
    success = false;
  }

  // TEST 8

  printf("\nTEST 8 \n");

  setStep(0.001);

  double pos15[] = {48,32,160};
  double vel15[] = {18,-27,21};
  double pos16[] = {2000,2000,2000};
  double vel16[] = {-5838,-5931,-5499};

  setParticleInfo(p1,pos15,vel15,acc,att);
  setParticleInfo(p2,pos16,vel16,acc,att);

  for (int i = 0; i < 1000; ++i)
  {
    if (!handleCollision(p1,p2))
    {
      move(p1,0.001);
      move(p2,0.001);
    }
  }

  double res15[] = {-3848,-3915.8,-3510.7};
  double res16[] = {55.95449,20.13987,3410.7};

  if (!(compareArray(p1->position,res15,3) && compareArray(p2->position,res16,3)))
  {
    printf("Test8: Error\n");
    success = false;
  }

  freeParticle(p1);
  freeParticle(p2);

  return success;
}

bool testLennardJones(){
  printf("\nLENNARD-JONES TESTING\n");

  // TEST 1


  bool success = true;

  printf("\nTEST 1 \n");

  struct particle *p1 = createParticle();
  struct particle *p2 = createParticle();

  double pos1[] = {0,0,0};
  double vel1[] = {0,0,0};
  double pos2[] = {4E-10,0,0};
  double vel2[] = {0,0,0};
  double acc[] = {0,0,0};
  double att[] = {40,0,188,0.997,3.4};

  double *response;

  setParticleInfo(p1,pos1,vel1,acc,att);
  setParticleInfo(p2,pos2,vel2,acc,att);

  response = lennardJonesPotentialForce(p1,p2);

  double res1[] = {-0.554328,0,0};

  if(!compareArray(response,res1,1))
  {
    printf("Test1: Error\n");
    success = false;
  }

  free(response);

  // TEST 2

  printf("\nTEST 2 \n");

  double pos3[] = {0,0,0};
  double vel3[] = {0,0,0};
  double pos4[] = {1,1,1};
  double vel4[] = {0,0,0};

  double *response2;

  setParticleInfo(p1,pos3,vel3,acc,att);
  setParticleInfo(p2,pos4,vel4,acc,att);

  response2 = lennardJonesPotentialForce(p1,p2);

  double val = 13718431877.851925;

  double res2[] = {val*pow(3,-2),val*pow(3,-2),val*pow(3,-2)};

  if(!compareArray(response2,res2,1))
  {
    printf("Test2: Error\n");
    success = false;
  }

  free(response2);

  // TEST 3

  setStep(0.001);

  printf("\nTEST 3 \n");

  double pos5[] = {0,0,0};
  double vel5[] = {5E-6,0,0};
  double pos6[] = {1E-5,0,0};
  double vel6[] = {0,0,0};

  double *response3;

  setParticleInfo(p1,pos5,vel5,acc,att);
  setParticleInfo(p2,pos6,vel6,acc,att);

  for (int i = 0; i < 10; ++i)
  {
    response3 = lennardJonesPotentialForce(p1,p2);
    scalarMultiplication(1/uToKg(p1->attributes->mass),response3);
    modifyVector(p2->acceleration,response3);
    scalarMultiplication(-1, response3);
    modifyVector(p1->acceleration, response3);
    if (!handleCollision(p1,p2))
    {
      move(p1,0.001);
      move(p2,0.001);
    }
    free(response3);
  }

  printf("Not supposed to be far appart:\n");
  printVector(p1->position);
  printVector(p2->position);


  freeParticle(p1);
  freeParticle(p2);

  return success;
}

bool testElectroStatic(){
  printf("\nELECTROSTATIC TESTING\n");

  // TEST 1

  bool success = true;

  printf("\nTEST 1 \n");

  struct particle *p1 = createParticle();
  struct particle *p2 = createParticle();

  double pos1[] = {0,0,0};
  double vel1[] = {0,0,0};
  double pos2[] = {1E-14,0,0};
  double vel2[] = {0,0,0};
  double acc[] = {0,0,0};
  double att1[] = {1,1,180,1,1};
  double att2[] = {1,-1,180,1,1};

  double *response;

  setParticleInfo(p1,pos1,vel1,acc,att1);
  setParticleInfo(p2,pos2,vel2,acc,att1);

  response = electroStaticForce(p1, p2);

  double res1[] = {2.307145,0,0};

  if(!compareArray(response,res1,1))
  {
    printf("Test1: Error\n");
    success = false;
  }

  free(response);

  // TEST 2

  printf("\nTEST 2 \n");

  setStep(0.001);

  double pos3[] = {0,0,0};
  double vel3[] = {0,0,0};
  double pos4[] = {1,0,0};
  double vel4[] = {0,0,0};

  double *response2;

  setParticleInfo(p1,pos3,vel3,acc,att1);
  setParticleInfo(p2,pos4,vel4,acc,att2);

  for (int i = 0; i < 1000000; ++i)
  {
    response2 = electroStaticForce(p1,p2);
    scalarMultiplication(1/uToKg(p1->attributes->mass), response2);
    modifyVector(p2->acceleration,response2);
    scalarMultiplication(-1, response2);
    modifyVector(p1->acceleration, response2);
    if (!handleCollision(p1,p2))
    {
      move(p1,0.001);
      move(p2,0.001);
    }
    free(response2);
  }

  printf("Not supposed to be far appart:\n");
  printVector(p1->position);
  printVector(p2->position);

  freeParticle(p1);
  freeParticle(p2);

  return success;
}

bool testLorentzForce(){
  printf("\nMAGNETIC TESTING\n");

  // TEST 1

  bool success = true;

  printf("\nTEST 1 \n");

  setStep(0.001);

  struct particle *p1 = createParticle();
  struct magneticField *magFi = createMagneticField();
  struct electricField *elFi = createElectricField();

  double pos1[] = {0,0,0};
  double vel1[] = {1,0,0};
  double magN[] = {0,1,0};
  double elN[] = {0,0,0};
  double acc[] = {0,0,0};
  double att1[] = {1,1,120,1,1};

  setParticleInfo(p1,pos1,vel1,acc,att1);

  setElectricField(elFi,elN);
  setMagneticField(magFi, magN);

  for (int i = 0; i < 1000; ++i)
  {
    applyFieldAcceleration(p1, elFi, magFi);
    move(p1,0.001);
  }

  freeParticle(p1);
  freeElField(elFi);
  freeMagField(magFi);

  return success;
}

bool testForceApplication(){
  printf("\nFORCE APPLICATION TESTING\n");

  // TEST 1

  bool success = true;

  printf("\nTEST 1 \n");

  struct particle *p1 = createParticle();
  struct particle *p2 = createParticle();
  struct electricField *elFi = createElectricField();
  struct magneticField *magFi = createMagneticField();

  double pos1[] = {0,0,0};
  double vel1[] = {0.01,0,0};
  double pos2[] = {1,0,0};
  double vel2[] = {0,0,0};
  double acc[] = {0,0,0};
  double att1[] = {1,1,120,1,1};
  double att2[] = {1,-1,120,1,1};

  double elForce[] = {1,0,0};
  double magForce[] = {0,1,0};

  setParticleInfo(p1,pos1,vel1,acc,att1);
  setParticleInfo(p2,pos2,vel2,acc,att1);
  setElectricField(elFi, elForce);
  setMagneticField(magFi, magForce);

  applyFieldAcceleration(p1,elFi, magFi);

  double res1[] = {9.648668505E8,0,9.648668505E5};

  if(!compareArray(p1->acceleration,res1,3))
  {
    printf("Test2: Error\n");
    success = false;
  }

  // TEST 2

  printf("\nTEST 2 \n");

  double pos3[] = {0,0,0};
  double vel3[] = {0,0,0};
  double pos4[] = {1,0,0};
  double vel4[] = {0,0,0};

  setParticleInfo(p1,pos3,vel3,acc,att1);
  setParticleInfo(p2,pos4,vel4,acc,att2);

  applyInterAcceleration(p1, p2);

  double res2[] = {0.138939,0,0};

  if(!compareArray(p1->acceleration,res2,3))
  {
    printf("Test2: Error\n");
    success = false;
  }

  freeParticle(p1);
  freeParticle(p2);
  freeElField(elFi);
  freeMagField(magFi);

  return success;
}



int main(){

  int results[12];

  results[0] = testNorm();

  results[1] = testConnectPoints();

  results[2] = testNormalize();

  results[3] = testDotProduct();

  results[4] = testSolveAugMatrix();

  results[6] = testCollisionCheck();

  results[7] = testFullCollisionScenario();

  results[8] = testLennardJones();

  results[9] = testElectroStatic();

  results[10] = testLorentzForce();

  results[11] = testForceApplication();



  printf("\nSummary\n");

  printf("NormTest: %d\n",results[0]);

  printf("ConnectTest: %d\n",results[1]);

  printf("Normalization: %d\n",results[2]);

  printf("DotProduct: %d\n",results[3]);

  printf("SolveAugMatrix: %d\n",results[4]);

  printf("CollisionCheck: %d\n",results[6]);

  printf("FullCollisionScenario: %d\n",results[7]);

  printf("Lennard-Jones Potential: %d\n",results[8]);

  printf("ElectroStaticForce: %d\n",results[9]);

  printf("Magnetic Force: %d\n",results[9]);

  printf("Force application: %d\n",results[11]);

  return 0;
}