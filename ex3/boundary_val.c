#include "boundary_val.h"
#include "helper.h"
#include <math.h>

//set the boundary conditions
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  int il,
  int ir,
  int jb,
  int jt)
{
    if (il==1){
        V[0][jb-1]=-V[1][jb-1];
        for (int j=jb;j<=jt; j++){
            U[0][j]=0;
            V[0][j]=-V[1][j];
        }
    }
    if (ir==imax){
        V[imax+1][jb-1]=-V[imax][jb-1];
        for (int j=jb;j<=jt; j++){
            U[imax][j]=0;
            V[imax+1][j]=-V[imax][j];
        }
    }
    if (jb==1){
        U[il-1][0]=-U[il-1][1];
        for (int i=il;i<=ir; i++){
            U[i][0]=-U[i][1];
            V[i][0]=0;
        }
    }
    if (jt==jmax){
        U[il-1][jmax+1]=2.0-U[il-1][jmax];
        for (int i=il;i<=ir;i++){
            U[i][jmax+1]=2.0-U[i][jmax];
            V[i][jmax]=0;
        }
    }
}