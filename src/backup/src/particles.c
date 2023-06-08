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
int **createPartition(int size);

double currTimeStep;

int main(int argc, char *argv[]){

  FILE *fp;
  FILE *outPutp;
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

  char fileNameInput[100];
  char pathInput[200] = "./particleCollection/";

  char fileNameOutput[150];

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

  if (argc == 2 && !strcmp(argv[1], "-listParticles"))
  {

    printDirectory();
    free(elFi);
    free(magFi);
    free(posVec);
    free(velVec);
    free(accVec);
    return 0;

  }
  else if (argc == 13 || argc == 14 || argc == 15)
  {

    strcpy(fileNameInput, argv[1]);
    strcat(pathInput,fileNameInput);
    particleArraySize = atoi(argv[2]);
    timePeriod = atoi(argv[3]);

    velVec[0] = strtod(argv[4], &ptr);
    velVec[1] = strtod(argv[5], &ptr);
    velVec[2] = strtod(argv[6], &ptr);

    elNorm[0] = strtod(argv[7], &ptr);
    elNorm[1] = strtod(argv[8], &ptr);
    elNorm[2] = strtod(argv[9], &ptr);
    for (int i = 0; i < 3; ++i)
    {
      spawnSphere[i] = atoi(argv[2])/2;
    }
    magNorm[0] = strtod(argv[10], &ptr);
    magNorm[1] = strtod(argv[11], &ptr);
    magNorm[2] = strtod(argv[12], &ptr);
    modifyVector(spawnPoint,spawnSphere);
  }
  else
  {
    printf("Invalid Input!\n");
    strcat(pathInput,"proton.txt\0");

    particleArraySize = 128;
    timePeriod = 1;

    velVec[0] = 1;
    velVec[1] = 0;
    velVec[2] = 1;

    elNorm[0] = 0;
    elNorm[1] = 1;
    elNorm[2] = 0;
    for (int i = 0; i < 3; ++i)
    {
      spawnSphere[i] = atoi(argv[2])/2;
    }
    magNorm[0] = 1;
    magNorm[1] = 0;
    magNorm[2] = 1;
    modifyVector(spawnPoint,spawnSphere);
    argc = 15;
    
  }

  particleArray = (struct particle**)malloc(particleArraySize*sizeof(struct particle*));

  modifyVector(posVec, initialVec);
  modifyVector(accVec, initialVec);

  setElectricField(elFi, elNorm);
  setMagneticField(magFi,magNorm);

  fp = fopen(pathInput,"r");
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

    particleArray[i] = tempParticle;
    modifyVector(spawnPoint,spawnSphere);
  }

  int *handled = (int*)malloc(particleArraySize*sizeof(handled));
  double *elF = elFi->fieldVector;
  double *magF = magFi->fieldVector;
  int **particleIdxPartition = createPartition(particleArraySize);
  int chunk = ((particleArraySize*(particleArraySize-1))/2);

  for (int i = 0; i < particleArraySize; ++i)
  {
    handled[i] = 0;
  }

  for (i = 0; i < timePeriod*(1/currTimeStep); ++i)
  {
/*
    if (i % 100000 == 0 && argc == 15)
    {
      printf("\n%lfs:\n",i*currTimeStep);
    }
*/
    // apply the acceleration to all particles

    #pragma omp parallel shared(particleArraySize,particleArray,elF,magF, particleIdxPartition, chunk, i)
      #pragma omp for schedule(static,1) private(j,k)
        for (j = 0; j < particleArraySize; ++j)
        {
          applyFieldAcceleration(particleArray[j],elF,magF);
          int id = omp_get_thread_num();
          int total = omp_get_num_threads();
          int *ptr;
          //printf("id: %d, total: %d, chunk: %d, part: %d\n", id, total, chunk,chunk/total );

	  applyInterAcceleration(particleArray[j],particleArray,j,particleArraySize);

          for (k = chunk*id/total; k < chunk*(id+1)/total; ++k)
          {
            ptr = particleIdxPartition[k];
            if (ptr && handled[j] == 0 && handleCollision(particleArray[ptr[0]],particleArray[ptr[1]],ptr[0],handled))
            {
              handled[j] = 1;
	      handled[k-chunk*id/total] = 1;
            }
          }
          // check whether the movement of a particle has already been computed
          if (handled[j] == 0)
          {
            move(particleArray[j],currTimeStep);
          }
          handled[j] = 0;
          scalarMultiplication(0,particleArray[j]->acceleration);
        }

    if (i % 100000 == 0 && argc == 15)
    {
      if (argc == 15)
      {
        snprintf(fileNameOutput, 150, "./visualData/%s.csv.%d", argv[14], i/100000);
      } else {
        snprintf(fileNameOutput, 150, "./visualData/timeStep.csv.%d", i/100000);
      }

      outPutp = fopen(fileNameOutput, "w");
      fputs("x coord, y coord, z coord, scalar\n", outPutp);
      for (j = 0; j < particleArraySize; ++j)
      {
        //printAttributes(particleArray[j]);
        printPos(particleArray[j], outPutp, j);
      }
      fclose(outPutp);
    }

  }

  freeArrayOfParticles(particleArray,particleArraySize);
  free(handled);
  free(posVec);
  free(velVec);
  free(accVec);
  free(attr);
  for (int i = 0; i < chunk; ++i)
  {
    free(particleIdxPartition[i]);
  }
  free(particleIdxPartition);

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

int **createPartition(int size){
  int partitionSize = (size*(size-1))/2;
  if (partitionSize % omp_get_max_threads() != 0){
    partitionSize += omp_get_max_threads() - (partitionSize % omp_get_max_threads());
  }
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






