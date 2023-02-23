#include<stdio.h>
#include<stdbool.h>
#include<stdlib.h>
#include<math.h>

double** createMatrix(int num_row, int num_col){
    double** matrix;
    double* row;

    matrix = (double**)malloc(num_row*sizeof(matrix));

    for (size_t i = 0; i < num_row; i++)
    {
        row = (double*)malloc(num_col*sizeof(row));
        matrix[i] = row;
    }

    return matrix;
}

void fillMatrix(double** empty_matrix, int num_row, int num_col){
    for (size_t i = 0; i < num_row; i++)
    {
        for (size_t j = 0; j < num_col; j++)
        {
            scanf("%lf",&empty_matrix[i][j]);
        }
        
    }
}

void printMatrix(double** empty_matrix, int num_row, int num_col){
    for (size_t i = 0; i < num_row; i++)
    {
        for (size_t j = 0; j < num_col; j++)
        {
            printf("%d ",empty_matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

double** sumMatrix(double** matrixA, double** matrixB, int num_row, int num_col){
    double** matrixC;
    matrixC = createMatrix(num_row, num_col);
    for (size_t i = 0; i < num_row; i++)
    {
        for (size_t j = 0; j < num_col; j++)
        {
            matrixC[i][j] = matrixA[i][j]+matrixB[i][j];
        }
        
    }
    return matrixC;
}

double** multiMatrix(double** matrixA, double** matrixB, int* sizeA, int* sizeB){
    
    double** matrixC, sum;
    matrixC = createMatrix(sizeA[0],sizeB[1]);
    for (size_t i = 0; i < sizeA[0]; i++)
    {
        for (size_t j = 0; j < sizeB[1]; j++)
        {
            sum = 0;
            for (size_t k = 0; k < sizeA[1]; k++)
            {
                sum += matrixA[i][k]*matrixB[k][j];
            }
            matrixC[i][j] = sum;
        }   
    }

    return matrixC;
}

double dotProduct(double** vectorA, double** vectorB){
    double result=0;
    for (size_t i = 0; i < 3; i++)
    {
        result += vectorA[i][0] + vectorB[i][0];
    }
    return result;
}

long double norm(double** vector){
    return sqrt(pow(vector[0][0],2) + pow(vector[1][0],2) + pow(vector[2][0],2));
}

double** crossProduct(double** vectorA, double** vectorB){
    double** new_vector = malloc(sizeof(new_vector));
    double *x,*y,*z = malloc(3*sizeof(new_vector));

    *x = vectorA[1][0]*vectorB[2][0]-vectorA[2][0]*vectorB[1][0];
    *y = -(vectorA[0][0]*vectorB[2][0]-vectorA[2][0]*vectorB[0][0]);
    *z = vectorA[0][0]*vectorB[1][0]-vectorA[1][0]*vectorB[0][0];

    new_vector[0] = x;
    new_vector[1] = y;
    new_vector[2] = z;

    return new_vector;
}


double** copyMatrix(double** matrix, int size){
    double** newMatrix = (double**)malloc(size*sizeof(matrix));
    double* row;
    for (size_t i = 0; i < size; i++)
    {
        row = (double*)malloc(size*sizeof(row));
        newMatrix[i] = row;
        for (size_t j = 0; j < size; j++)
        {
            newMatrix[i][j] = matrix[i][j];
        }
    }
    return newMatrix;
}

double det(double** matrix, int size){
    double determinant;
    determinant = matrix[0][0] * (matrix[1][1]*matrix[2][2]-matrix[2][1]*matrix[1][2]) - 
    matrix[0][1] * (matrix[1][0]*matrix[2][2]-matrix[2][0]*matrix[1][2]) + 
    matrix[0][2] * (matrix[1][0]*matrix[2][1]-matrix[2][0]*matrix[1][1]);
    return determinant;
}

struct particle
{
    double** position;
    double** movement;
    double weight;
    double charge;
    double radius;
    double Hamaker_coef;
};



int main(){
    
    
    double** matrix_A, d;
    int sA[2] = {3,3};
    matrix_A = createMatrix(3,3);
    fillMatrix(matrix_A,3,3);
    printMatrix(matrix_A, 3,3);

    d = det(matrix_A, 3);

    printf("\n%lf",d);
    

    
    return 0; 
}