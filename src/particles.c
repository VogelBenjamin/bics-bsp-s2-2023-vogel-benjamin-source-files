#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<dirent.h>
#include<time.h>
#include"linearAlgebra.h"
#include"electromagnetism.h"
#include "omp.h"

void printDirectory();

double currTimeStep;

int main(int argc, char *argv[]){

  FILE *fp;
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

  char fileName[100];
  char path[200] = "./particleCollection/";

  struct particle **particleArray;
  struct particle *tempParticle;
  struct particle *p;
  struct electricField *elFi = createElectricField();
  struct magneticField *magFi = createMagneticField();

  currTimeStep = getStep();

  double *posVec = createVector();
  double *velVec = createVector();
  double *accVec = createVector();
  double *attr;
  double elNorm[3];
  double magNorm[3];
  double initialVec[3] = {0,0,0};
  double spawnSphere[3];
  double spawnPoint[3];

  if (!strcmp(argv[1], "-listParticles") && argc == 2)
  {

    printDirectory();
    free(elFi);
    free(magFi);
    free(posVec);
    free(velVec);
    free(accVec);
    return 0;

  }
  else if (argc == 7)
  {

    strcpy(fileName, argv[1]);
    strcat(path,fileName);
    particleArraySize = atoi(argv[2]);
    timePeriod = atoi(argv[3]);
    velVec[0] = strtod(argv[4], &ptr);
    for (int i = 0; i < 3; ++i)
    {

      elNorm[i] = strtod(argv[5], &ptr);
      spawnSphere[i] = atoi(argv[2]);

    }
    magNorm[1] = strtod(argv[6], &ptr);
    modifyVector(spawnPoint,spawnSphere);
  }
  else
  {
    printf("%d\n", argc);
    for (int i = 0; i < argc; ++i)
    {
      printf("%s\n",argv[i]);
    }
    printf("Invalid Input!\n");
    return 1;
  }

  particleArray = (struct particle**)malloc(particleArraySize*sizeof(struct particle*));

  modifyVector(posVec, initialVec);
  modifyVector(accVec, initialVec);

  setElectricField(elFi, elNorm);
  setMagneticField(magFi,magNorm);

  fp = fopen(path,"r");

  attr = getParticleAttributes(fp);
  fclose(fp);

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

    applyFieldAcceleration(tempParticle,elFi,magFi);
    particleArray[i] = tempParticle;
    modifyVector(spawnPoint,spawnSphere);
  }

  int *handled = (int*)malloc(particleArraySize*sizeof(handled));
  for (int i = 0; i < particleArraySize; ++i)
  {
    handled[i] = 0;
  }

  //#pragma omp parallel shared(particleArraySize,particleArray,elFi,magFi,currTimeStep)
    //#pragma omp for schedule(dynamic) private(i,j,k)

        for (i = 0; i < timePeriod*1000; ++i)
        {
          if (i % 100 == 0)
          {
            printf("\n%lfs:\n",i*0.001);
          }
          // apply the acceleration to all particles
          for (j = 0; j < particleArraySize; ++j)
          {
            applyFieldAcceleration(particleArray[j],elFi,magFi);
          }

          // apply interaction forces
          for (j = 0; j < particleArraySize; ++j)
          {

              // checks all the remaining possible collisions and performs particle displacement
              for (k = j+1; k < particleArraySize; ++k)
              {
                applyInterAcceleration(particleArray[j],particleArray[k]);
              }
          }

          // check collisions and move the particles
          for (j = 0; j < particleArraySize; ++j)
          {

              // checks all the remaining possible collisions and performs particle displacement
              for (k = j+1; k < particleArraySize; ++k)
              {

                if (handleCollision(particleArray[j],particleArray[k]))
                {
                  handled[j] = 1;
                  handled[k] = 1;
                  break;
                }

              }

              // check whether the movement of a particle has already been computed
              if (handled[j] == 0)
              {
                move(particleArray[j],currTimeStep);
              }



            if (i % 100 == 0)
            {
              printAttributes(particleArray[j]);
            }
          }

        }




  freeArrayOfParticles(particleArray,particleArraySize);

  free(posVec);
  free(velVec);
  free(accVec);
  free(attr);

  freeElField(elFi);
  freeMagField(magFi);
  
  return 0;
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






