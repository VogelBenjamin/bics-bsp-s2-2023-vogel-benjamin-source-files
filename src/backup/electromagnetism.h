#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <stdbool.h>
#include "omp.h"


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
