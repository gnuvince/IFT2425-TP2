/*  IFT2425 - TP2 */
/* Vincent Foley-Bourgon (FOLV08078309) */
/* Eric Thivierge (THIE09016601) */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 10
#define M 0.36
#define EPSILON 0.0000000001

/*
  Structure pour représenter une matrice:
  - elems: pointeurs vers les rangées
  - start: pointeur vers le début des données
           (pour libérer la mémoire et copier les données)
  - rows, cols: dimensions de la matrice.
 */
typedef struct {
    float** elems;
    float*  start;
    int     rows, cols;
} matrix_t;


/* Allocation dynamique d'une nouvelle matrice.  Tous les éléments
 * sont initialisés à 0. */
matrix_t* NewMatrix(int rows, int cols) {
    matrix_t* mat = malloc(sizeof(matrix_t));
    if (mat == NULL) {
        fprintf(stderr, "memoire insuffisante\n");
        abort();
    }

    float** elems = malloc(sizeof(float*) * rows);
    if (elems == NULL) {
        fprintf(stderr, "memoire insuffisante\n");
        abort();
    }

    float* data = calloc(rows*cols, sizeof(float));
    if (data == NULL) {
        fprintf(stderr, "memoire insuffisante\n");
        abort();
    }

    for (int i = 0; i < rows; ++i) {
        elems[i] = data + (i*cols);
    }

    mat->elems = elems;
    mat->start = data;
    mat->rows = rows;
    mat->cols = cols;

    return mat;
}

/* Libération de la mémoire utilisée par une matrice. */
void FreeMatrix(matrix_t* mat) {
    free(mat->start);
    free(mat->elems);
    free(mat);
}


void FillMatrix(matrix_t* mat, float val) {
  for (int i = 0; i < mat->rows; ++i)
    for (int j = 0; j < mat->cols; ++j)
        mat->elems[i][j] = val;
}


/* Créer une matrice ayant le format penta-diagonal:
  1    m   m^2 0             ... 0
  m    1   m   m^2 0         ... 0
  m^2  m   1   m   m^2 0     ... 0
  0    m^2 m   1   m   m^2 0 ... 0
  .
  .
  .
  0 ...                 0  m^2 m 1
 */
void MakePentadiagonalMatrix(matrix_t* matrix, float m) {
    float m2 = m*m;

    // Do nothing if the matrix isn't square or if m not in [0..1].
    if ((matrix->rows != matrix->cols) || (m < 0.0) || (m > 1.0))
        return;

    matrix->elems[0][0] = 1.0;
    matrix->elems[0][1] = m;
    matrix->elems[0][2] = m2;
    matrix->elems[1][0] = m;
    matrix->elems[1][1] = 1.0;
    matrix->elems[1][2] = m;
    matrix->elems[1][3] = m2;
    if ( matrix->rows > 4) {
        for (int i = 2; i < matrix->rows - 2; ++i) {
            matrix->elems[i][i-2] = m2;
            matrix->elems[i][i-1] = m;
            matrix->elems[i][i] = 1.0;
            matrix->elems[i][i+1] = m;
            matrix->elems[i][i+2] = m2;
        }
    }
    matrix->elems[matrix->rows-2][matrix->cols-4] = m2;
    matrix->elems[matrix->rows-2][matrix->cols-3] = m;
    matrix->elems[matrix->rows-2][matrix->cols-2] = 1.0;
    matrix->elems[matrix->rows-2][matrix->cols-1] = m;
    matrix->elems[matrix->rows-1][matrix->cols-3] = m2;
    matrix->elems[matrix->rows-1][matrix->cols-2] = m;
    matrix->elems[matrix->rows-1][matrix->cols-1] = 1.0;
}


/* Returns true if mat is diagonally dominant */
int DiagonallyDominant(matrix_t* mat) {
    float linetotal;
    
    for (int i = 0; i < mat->rows; ++i) {
        linetotal = 0.0;
        for (int j = 0; j < mat->cols; ++j) {
            if (i != j)
                linetotal += fabsf(mat->elems[i][j]);
        }
        if (fabs(mat->elems[i][i]) <= linetotal)
            return 0;
    }
    return 1;
}

float VectorNorm1(matrix_t* vector) {
    float norm = 0.0f;

    for (int i = 0; i < vector->rows; ++i)
        norm += fabs(vector->elems[i][0]);

    return norm;
}

float MatrixNorm1(matrix_t* mat) {
    float norm = 0.0f;
    float max = 0.0f;

    for (int j = 0; j < mat->cols; ++j) {
        norm = 0.0;
        for (int i = 0; i < mat->rows; ++i) {
            norm += fabs(mat->elems[i][j]);
        }
        if (norm > max)
            max = norm;
    }
    return max;
}

int MatrixAddCompatible(matrix_t* A, matrix_t* B, matrix_t* C) {
    if ((A->rows != B->rows) || (A->rows != C->rows) ||
        (A->cols != B->cols) || (A->cols != C->cols) ||
        (B->rows != C->rows) || (B->cols != C->cols))
        return 0;
    
    return 1;
}

void MatrixSub(matrix_t* A, matrix_t* B, matrix_t* C) {
    if (!MatrixAddCompatible(A, B, C))
        return;
    
    for (int i = 0; i < C->rows; ++i)
        for (int j = 0; j < C->cols; ++j)
            C->elems[i][j] = A->elems[i][j] - B->elems[i][j];    
}

void MatrixAdd(matrix_t* A, matrix_t* B, matrix_t* C) {
    if (!MatrixAddCompatible(A, B, C))
        return;

    for (int i = 0; i < C->rows; ++i)
        for (int j = 0; j < C->cols; ++j)
            C->elems[i][j] = A->elems[i][j] + B->elems[i][j];
}

void MatrixScalarMult(matrix_t* A, float scal, matrix_t* C) {
    for (int i = 0; i < C->rows; ++i)
        for (int j = 0; j < C->cols; ++j)
            C->elems[i][j] = A->elems[i][j] * scal;    
}

void MatrixMult(matrix_t* A, matrix_t* B, matrix_t* C) {
    if ((A->cols != B->rows) || (C->rows != A->rows) || (C->cols != B->cols))
        return;

    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < B->cols; j++) {
            float temp = 0;
            for (int k = 0; k < B->rows; k++)
                temp += A->elems[i][k] * B->elems[k][j];
            C->elems[i][j] = temp;
        }
}


/*
 * vector n <- A * [1, 1, ..., 1]^T
 */
void MakeBVector(matrix_t* A, matrix_t* vector) {
    matrix_t* v1 = NewMatrix(N, 1);

    FillMatrix(v1, 1.0);
    MatrixMult(A, v1, vector);
    
    FreeMatrix(v1);
}


/* Afficher une matrice. */
void PrintMatrix(matrix_t* matrix) {
    for (int i = 0; i < matrix->rows; ++i) {
        for (int j = 0; j < matrix->cols; ++j) {
            printf("%+.2f ", matrix->elems[i][j]);
        }
        putchar('\n');
    }
    putchar('\n');
}


/* Copier les données d'une matrice dans une autre. */
void CopyMatrix(matrix_t* dst, matrix_t* src) {
    if (src->rows != dst->rows || src->cols != dst->cols)
        return;

    memcpy(dst->start, src->start, sizeof(float) * src->rows * src->cols);
}



matrix_t* MakeI(int n) {
  matrix_t* I = NewMatrix(n, n);

  for (int i = 0; i < n; i++)
    I->elems[i][i] = 1;

  return I;
}

int MatrixEq(matrix_t* A, matrix_t* B) {
  if ((A->rows != B-> rows) || (A->cols != B->cols))
    return 0;

  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->cols; j++)
      if (A->elems[i][j] != B->elems[i][j])
        return 0;

  return 1;
}


/* x contains x0 for the iterative evaluation */
void SolveJacobi(matrix_t* A, matrix_t* b, matrix_t* x) {
    matrix_t* x1 = NewMatrix(N, 1);
    matrix_t* x2 = NewMatrix(N, 1);
    matrix_t* x3 = NewMatrix(N, 1);
    matrix_t* temp;
    int iter = 0;

    CopyMatrix(x1, x);
    do {
        FillMatrix(x2, 0.0f);        
        iter += 1;
        for (int i = 0; i < x->rows; ++i) {
            for (int j = 0; j < A->cols; ++j) {
                if (i != j)
                    x2->elems[i][0] += (A->elems[i][j] / A->elems[i][i]) * x1->elems[j][0]; 
            }
            x2->elems[i][0] = (b->elems[i][0] / A->elems[i][i]) -  x2->elems[i][0];
        }
        MatrixSub(x2, x1, x3);
        temp = x1;
        x1 = x2;
        x2 = temp;
    } while (VectorNorm1(x3) > EPSILON);

    printf("iter: %d\n", iter);
    CopyMatrix(x, x1);
    FreeMatrix(x3);
    FreeMatrix(x2);
    FreeMatrix(x1);
}

/* x contains x0 for the iterative evaluation */
void SolveGaussSeidel(matrix_t* A, matrix_t* b, matrix_t* x) {
    matrix_t* x1 = NewMatrix(N, 1);
    matrix_t* x2 = NewMatrix(N, 1);
    matrix_t* x3 = NewMatrix(N, 1);
    matrix_t* temp;
    float accumul = 0.0f;
    int iter = 0;

    CopyMatrix(x1, x);
    do {
        FillMatrix(x2, 0.0f);        
        iter += 1;
        for (int i = 0; i < x->rows; ++i) {
            accumul = 0.0f;
            for (int j = 0; j < i; ++j) {
                accumul += (A->elems[i][j] / A->elems[i][i]) * x2->elems[j][0]; 
            }
            for (int j = (i + 1); j < A->cols; ++j) {
                accumul += (A->elems[i][j] / A->elems[i][i]) * x1->elems[j][0]; 
            }
            x2->elems[i][0] = (b->elems[i][0] / A->elems[i][i]) - accumul;
        }
        MatrixSub(x2, x1, x3);
        temp = x1;
        x1 = x2;
        x2 = temp;
    } while (VectorNorm1(x3) > EPSILON);

    printf("iter: %d\n", iter);    
    CopyMatrix(x, x1);
    FreeMatrix(x3);
    FreeMatrix(x2);
    FreeMatrix(x1);
}

int main(void) {
    matrix_t* A = NewMatrix(N, N);
    matrix_t* B = NewMatrix(N, N);
    matrix_t* C = NewMatrix(N, N);
    matrix_t* x = NewMatrix(N, 1);
    matrix_t* b = NewMatrix(N, 1);

    MakePentadiagonalMatrix(A, M);
    MakeBVector(A, b);

    printf("A:\n");
    PrintMatrix(A);
    if (DiagonallyDominant(A))
        printf("Matrice diagonalement dominante\n\n");
    else
        printf("Matrice non diagonalement dominante\n\n");
    printf("b:\n");
    PrintMatrix(b);

    FillMatrix(x, 0.0f);    // x contains x0    
    SolveJacobi(A, b, x);
    printf("x Jacobi:\n");
    PrintMatrix(x);    

    FillMatrix(x, 0.0f);    // x contains x0    
    SolveGaussSeidel(A, b, x);
    printf("x GaussSeidel:\n");
    PrintMatrix(x);    
    
    FreeMatrix(b);
    FreeMatrix(x);
    FreeMatrix(C);
    FreeMatrix(B);
    FreeMatrix(A);

    return 0;
}
