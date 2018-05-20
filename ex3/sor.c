#include "sor.h"
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int il,
  int ir,
  int jb,
  int jt
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* SOR iteration */
  for(j = jb; j <= jt; j++) {
    for(i = il; i<=ir; i++) {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }

  /* compute the residual */
  rloc = 0;
  for(j = jb; j <= jt; j++) {
    for(i = il; i <= ir; i++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  rloc = rloc/(imax*jmax);
  //rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;


  /* set boundary values */
  if (il==1){
    for(j = jb; j <= jt; j++) {
      P[0][j] = P[1][j];
    }
  }
  if (ir==imax){
    for(j = jb; j <= jt; j++) {
      P[imax+1][j] = P[imax][j];
    }
  }
  if (jb==1){
    for(i = il; i <= ir; i++) {
      P[i][0] = P[i][1];
    }
  }
  if (jt==jmax){
    for(i = il; i <= ir; i++) {
      P[i][jmax+1] = P[i][jmax];
    }
  }
  
}

