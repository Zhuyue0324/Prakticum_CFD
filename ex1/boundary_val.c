#include "boundary_val.h"
#include "helper.h"
#include <math.h>

void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V)
{
    for (int j=1;j<=jmax; j++){
        U[0][j]=0;
        V[0][j]=-V[1][j];
        U[imax][j]=0;
        V[imax+1][j]=-V[imax][j];
    } 
    for (int i=1;i<=imax; i++){
        U[i][0]=-U[i][1];
        V[i][0]=0;
        U[i][jmax+1]=2.0-U[i][jmax];
        V[i][jmax]=0;
    }
}