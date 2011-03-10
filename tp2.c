/*  IFT2425 - TP2 */
/* Vincent Foley-Bourgon (FOLV08078309) */
/* Eric Thivierge (THIE09016601) */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Paramètres donnés dans l'énoncé. */
#define N 5
#define M 0.3
#define EPSILON 10e-10

/*
  Structure pour représenter une matrice:
  - elems: pointeurs vers les rangées
  - start: pointeur vers le début des données
           (pour libérer la mémoire et copier les données)
  - rows, cols: dimensions de la matrice.
 */
typedef struct {
    double** elems;
    double*  start;
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

    double** elems = malloc(sizeof(double*) * rows);
    if (elems == NULL) {
        fprintf(stderr, "memoire insuffisante\n");
        abort();
    }

    double* data = calloc(rows*cols, sizeof(double));
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


void FillMatrix(matrix_t* mat, double val) {
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

void MakePentadiagonalMatrix(matrix_t* matrix, double m) {
    // Do nothing if the matrix isn't square or if m not in [0..1].
    if ((matrix->rows != matrix->cols) || (m < 0.0) || (m > 1.0))
        return;

    double m2 = m*m;

    for (int i = 0; i < matrix->rows; ++i) {
        for (int j = 0; j < matrix->cols; ++j) {
            if (i == j) matrix->elems[i][j] = 1;
            else if (fabs(i - j) == 1) matrix->elems[i][j] = m;
            else if (fabs(i - j) == 2) matrix->elems[i][j] = m2;
        }
    }
}



/* Returns true if mat is diagonally dominant */
int DiagonallyDominant(matrix_t* mat) {
    double linetotal;

    for (int i = 0; i < mat->rows; ++i) {
        linetotal = 0.0;
        for (int j = 0; j < mat->cols; ++j) {
            linetotal += fabsf(mat->elems[i][j]) * (i != j);
        }
        if (fabs(mat->elems[i][i]) <= linetotal)
            return 0;
    }
    return 1;
}

double VectorNorm1(matrix_t* vector) {
    double norm = 0.0f;

    for (int i = 0; i < vector->rows; ++i)
        norm += fabs(vector->elems[i][0]);

    return norm;
}

double MatrixNorm1(matrix_t* mat) {
    double norm = 0.0f;
    double max = 0.0f;

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
    return ((A->rows == B->rows) || (A->rows == C->rows) ||
            (A->cols == B->cols) || (A->cols == C->cols) ||
            (B->rows == C->rows) || (B->cols == C->cols));
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

void MatrixScalarMult(matrix_t* A, double scal, matrix_t* C) {
    for (int i = 0; i < C->rows; ++i)
        for (int j = 0; j < C->cols; ++j)
            C->elems[i][j] = A->elems[i][j] * scal;
}

void MatrixMult(matrix_t* A, matrix_t* B, matrix_t* C) {
    if ((A->cols != B->rows) || (C->rows != A->rows) || (C->cols != B->cols))
        return;

    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < B->cols; j++) {
            double temp = 0;
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

    memcpy(dst->start, src->start, sizeof(double) * src->rows * src->cols);
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


int SolveJacobi(matrix_t* A, matrix_t* b, matrix_t* x) {
    matrix_t* x_prev = NewMatrix(x->rows, 1);
    matrix_t* vect_difference = NewMatrix(x->rows, 1);

    double norm_difference;
    int iterations = 1;

    FillMatrix(x_prev, 0);
    do {
        for (int i = 0; i < A->rows; ++i) {
            double val = b->elems[i][0] / A->elems[i][i];
            for (int j = 0; j < A->cols; ++j) {
                if (i != j) {
                    val -= (A->elems[i][j] / A->elems[i][i]) * x_prev->elems[j][0];
                }
            }
            x->elems[i][0] = val;
        }
        MatrixSub(x, x_prev, vect_difference);
        norm_difference = fabsf(VectorNorm1(vect_difference));
        CopyMatrix(x_prev, x);
        iterations++;
    } while (norm_difference > EPSILON);

    FreeMatrix(x_prev);
    FreeMatrix(vect_difference);
    return iterations;
}


int SolveGaussSeidel(matrix_t* A, matrix_t* b, matrix_t* x) {
    double norm_difference;
    double prev_norm;
    int iterations = 1;

    FillMatrix(x, 0);
    prev_norm = VectorNorm1(x);
    do {
        for (int i = 0; i < A->rows; ++i) {
            double val = b->elems[i][0] / A->elems[i][i];
            for (int j = 0; j < A->cols; ++j) {
                if (i != j) {
                    val -= (A->elems[i][j] / A->elems[i][i]) * x->elems[j][0];
                }
            }
            x->elems[i][0] = val;
        }
        double norm = VectorNorm1(x);
        norm_difference = fabsf(norm - prev_norm);
        prev_norm = norm;
        iterations++;
    } while (norm_difference > EPSILON);

    return iterations;
}


int main(void) {
    matrix_t* A = NewMatrix(N, N);
    matrix_t* B = NewMatrix(N, N);
    matrix_t* C = NewMatrix(N, N);
    matrix_t* x = NewMatrix(N, 1);
    matrix_t* b = NewMatrix(N, 1);

    puts("QUESTION 1");
    puts("==========");
    MakePentadiagonalMatrix(A, M);
    MakeBVector(A, b);
    puts("A:");
    PrintMatrix(A);
    puts("b:");
    PrintMatrix(b);


    puts("QUESTION 2");
    puts("==========");
    puts("Les méthodes de Jacobi et Gauss-Seidel convergent si A est diagonalement dominante.");
    if (DiagonallyDominant(A))
        printf("La matrice A est diagonalement dominante\n\n");
    else
        printf("La matrice A n'est pas diagonalement dominante\n\n");

    puts("QUESTION 3");
    puts("==========");
    FillMatrix(x, 0.0f);    // x contains x0
    int iter_jacobi = SolveJacobi(A, b, x);
    printf("Solution x par la méthode de Jacobi:\n");
    printf("Nombre d'itérations: %d\n", iter_jacobi);
    PrintMatrix(x);

    FillMatrix(x, 0.0f);    // x contains x0
    int iter_gauss_seidel = SolveGaussSeidel(A, b, x);
    printf("Solution x par la méthode de Gauss-Seidel:\n");
    printf("Nombre d'itérations: %d\n", iter_gauss_seidel);
    PrintMatrix(x);

    FreeMatrix(b);
    FreeMatrix(x);
    FreeMatrix(C);
    FreeMatrix(B);
    FreeMatrix(A);

    return 0;
}
