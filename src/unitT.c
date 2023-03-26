#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include <dirent.h>
#include"linearAlgebra.h"
#include"electromagnetism.h"
#include "omp.h"

int main(){
  struct particle *p1 = createParticle();
  struct particle *p2 = createParticle();

  double pos1[] = {0,0,0};
  double vel1[] = {1,0,0};
  double pos2[] = {1,0,0};
  double vel2[] = {0,0,0};
  double acc[] = {0,0,0};
  double att[] = {1,0,1,0};

  setParticleInfo(p1,pos1,vel1,acc,att);
  setParticleInfo(p2,pos2,vel2,acc,att);

  printAttributes(p1);
  printAttributes(p2);

  calcCollision(p1,p2);

  printAttributes(p1);
  printAttributes(p2);



  return 0;
}