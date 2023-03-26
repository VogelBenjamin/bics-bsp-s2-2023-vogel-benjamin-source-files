#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "omp.h"

struct particleInfo {
  double mass;
  double charge;
  double radius;
  double epsilon;
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
void move(struct particle *particle,double timeStep);
void applyAcceleration(struct particle *p, struct electricField *eF, struct magneticField *mF);
void checkCollision(struct particle *p1, struct particle *p2);
void calcCollision(struct particle *p1, struct particle *p2);
double *lennardJonesPotential(struct particle *p1, struct particle *p2);
double *electroStaticForce(struct particle *p1, struct particle *p2);