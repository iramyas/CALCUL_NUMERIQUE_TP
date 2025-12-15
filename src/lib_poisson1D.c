/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
    int n = *la;
    int lab_ = *lab;

    /* mettre tout à 0 */
    for(int j = 0; j < n; j++){
        for(int i = 0; i < lab_; i++){
            AB[i + j*lab_] = 0.0;
        }
    }
    /* remplir la tridiagonale */
    for(int j = 0; j < n; j++){
        /* diag principale */
        int i = j;
        AB[indexABCol(i, j, lab)] = 2.0;

        /* sous-diagonale (i = j-1) */
        if(j > 0){
            i = j - 1;
            AB[indexABCol(i, j, lab)] = -1.0;
        }

        /* sur-diagonale (i = j+1) */
        if(j < n - 1){
            i = j + 1;
            AB[indexABCol(i, j, lab)] = -1.0;
        }
    }
}


void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
    int n = *la;
    int lab_ = *lab;

    for(int j = 0; j < n; j++){
        for(int i = 0; i < lab_; i++){
            AB[i + j*lab_] = 0.0;
        }
    }

    for(int j = 0; j < n; j++){
        int i = j;
        AB[indexABCol(i, j, lab)] = 1.0;
    }
}


void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
    int n = *la;
    int i;

    /* on part de RHS = 0 */
    for (i = 0; i < n; i++){
        RHS[i] = 0.0;
    }
    double h = 1.0 / (double)(n + 1);

    /* contributions des bords (schéma -u'' ≈ (2u_i - u_{i-1} - u_{i+1}) / h^2) */
    RHS[0]     += (*BC0) / (h*h);      /* i = 1 -> terme en u_0 */
    RHS[n-1]   += (*BC1) / (h*h);      /* i = n -> terme en u_{n+1} */
}


void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
    int n = *la;
    int i;
    double T0 = *BC0;
    double T1 = *BC1;

    for (i = 0; i < n; i++){
        double x = X[i];
        EX_SOL[i] = T0 + (T1 - T0) * x;
    }
}
 

void set_grid_points_1D(double* x, int* la){
  int n = *la;
  double h = 1.0 / (double)(n + 1);

  for(int i = 0; i < n; i++){
    x[i] = (i + 1) * h;
  }
}

double relative_forward_error(double* x, double* y, int* la){
    int n = *la;
    int i;
    double num = 0.0;
    double den = 0.0;

    for (i = 0; i < n; i++){
        double diff = x[i] - y[i];
        num += diff * diff;
        den += y[i] * y[i];
    }

    if (den == 0.0){
        return 0.0;
    }
    return sqrt(num / den);
}

int indexABCol(int i, int j, int *lab){
  int kv = *lab - 2;          
  int k = i - j + kv;         
  return k + j * (*lab);

}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  return *info;
}