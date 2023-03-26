#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include <dirent.h>
#include"linearAlgebra.h"
#include"electromagnetism.h"
#include "omp.h"

void printMenuText();
void printDirectory();

int main(){

  FILE *fp;

  int menuInput = 0;
  int particleArraySize;
  int i;
  int j;

  char fileName[100];
  char path[200] = "./particleCollection/";

  struct particle **particleArray;
  struct particle *tempParticle;
  struct particle *p;
  struct electricField *elFi = createElectricField();
  struct magneticField *magFi = createMagneticField();

  printMenuText();
  while(0==0){
    printf("Please input an instruction (1-3): ");
    scanf("%d",&menuInput);
    printf("\n");
    switch(menuInput){
      case 1:
        printf("Please input particle.txt you want to use: ");
        scanf("%s",fileName);
        strcat(path,fileName);
        break;
      case 2:
        printDirectory();
        break;
      case 3:
        return 0;
        break;
      default:
        return -1;
        break;
    }
    if (menuInput == 1){
      break;
    }
  }


  printf("Please input the number of particles to simulate: ");
  scanf("%d",&particleArraySize);
  particleArray = (struct particle**)malloc(particleArraySize*sizeof(struct particle*));


  double *posVec = createVector();
  double *velVec = createVector();
  double *accVec = createVector();
  double *elNorm = createVector();
  double *magNorm = createVector();

  printf("Please give the initial position of the particle: ");
  setVector(posVec);
  printf("Please give the initial velocity of the particle: ");
  setVector(velVec);
  printf("Please give the initial acceleration of the particle: ");
  setVector(accVec);
  printf("Please give the vector of electric field: ");
  setVector(elNorm);
  printf("Please give the vector of magnetic field: ");
  setVector(magNorm);

  setElectricField(elFi, elNorm);
  setMagneticField(magFi,magNorm);

  double *pos;
  double *vel;
  double *acc;
  double *attr;

  fp = fopen(path,"r");

  attr = getParticleAttributes(fp);
  fclose(fp);

  pos = createVector();
  vel = createVector();
  acc = createVector();

  modifyVector(pos,posVec);
  modifyVector(vel,velVec);
  modifyVector(acc,accVec);


  for (i = 0; i < particleArraySize; ++i)
  {
    tempParticle = createParticle();

    setParticleInfo(tempParticle, pos, vel, acc, attr);

    applyAcceleration(tempParticle,elFi,magFi);
    particleArray[i] = tempParticle;

  }

  free(pos);
  free(vel);
  free(acc);

  //#pragma omp parallel shared(particleArraySize,particleArray,elFi,magFi)
    //#pragma omp for schedule(static,200) private(i,j)
      for (i = 0; i < particleArraySize; ++i)
      {
        for (j = 1; j < 1001; ++j)
        {
          applyAcceleration(particleArray[i],elFi,magFi);
          move(particleArray[i],0.001);
          if (j % 100 == 0)
          {
            printAttributes(particleArray[i]);
          }
        }
      }


  freeArrayOfParticles(particleArray,particleArraySize);

  free(posVec);
  free(velVec);
  free(accVec);
  free(elNorm);
  free(attr);

  free(magNorm);

  freeElField(elFi);
  freeMagField(magFi);
  
  return 0;
}

void printMenuText(){
  printf("Menu:\n");
  printf("---------------------------\n");
  printf("1. Simulation\n");
  printf("2. Particle list\n");
  printf("3. Quit\n");
  printf("---------------------------\n" );
}

void printDirectory(){
  FILE *fp;
  DIR *folder;
  struct dirent *dir;
  char character;
  char directoryName[] = "./particleCollection/";
  char *fileLocation;
  folder = opendir("./particleCollection/");


  while((dir = readdir(folder)) != NULL){
    if (strcmp(dir->d_name,".")!=0 && strcmp(dir->d_name,"..")!=0)
    {
      fileLocation = (char *)malloc(50*sizeof(fileLocation));
      strcpy(fileLocation,directoryName);
      strcat(fileLocation, dir->d_name);
      fp = fopen(fileLocation, "r");
      if (fp != NULL){
        printf("%s:\n\n",dir->d_name);
        character = (char)fgetc(fp);
        while(character != EOF){
          printf("%c",character);
          character = (char)fgetc(fp);
        }
        printf("\n\n");
        fclose(fp);
      }

      free(fileLocation);

    }
  }
  closedir(folder);
}






