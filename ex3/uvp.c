#include "uvp.h"
#include "helper.h"
#include <math.h>


/**
 * Determines the value of U and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */
void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  int il,
  int ir,
  int jb,
  int jt
)
{
    for (int j=jb;j<=jt;j++){
        for (int i=il;i<ir;i++){
            F[i][j]=U[i][j]+dt*(
                ((U[i+1][j]-2.0*U[i][j]+U[i-1][j])/(dx*dx)+(U[i][j+1]-2.0*U[i][j]+U[i][j-1])/(dy*dy))/Re
                -(((U[i][j]+U[i+1][j])/2.0)*((U[i][j]+U[i+1][j])/2.0)-((U[i-1][j]+U[i][j])/2.0)*((U[i-1][j]+U[i][j])/2.0)
                    +alpha*((fabs(U[i][j]+U[i+1][j])/2.0)*((U[i][j]-U[i+1][j])/2.0)-(fabs(U[i-1][j]+U[i][j])/2.0)*((U[i-1][j]-U[i][j])/2.0))
                )/dx
                -(((V[i][j]+V[i+1][j])/2.0)*((U[i][j]+U[i][j+1])/2.0)-((V[i][j-1]+V[i+1][j-1])/2.0)*((U[i][j-1]+U[i][j])/2.0)
                    +alpha*((fabs(V[i][j]+V[i+1][j])/2.0)*((U[i][j]-U[i][j+1])/2.0)-(fabs(V[i][j-1]+V[i+1][j-1])/2.0)*((U[i][j-1]-U[i][j])/2.0))
                )/dy
                +GX
            );
        }
    }
    for (int j=jb;j<jt;j++){
        for (int i=il;i<=ir;i++){
            G[i][j]=V[i][j]+dt*(
                ((V[i+1][j]-2.0*V[i][j]+V[i-1][j])/(dx*dx)+(V[i][j+1]-2.0*V[i][j]+V[i][j-1])/(dy*dy))/Re
                -(((U[i][j]+U[i][j+1])/2.0)*((V[i][j]+V[i+1][j])/2.0)-((U[i-1][j]+U[i-1][j+1])/2.0)*((V[i-1][j]+V[i][j])/2.0)
                    +alpha*((fabs(U[i][j]+U[i][j+1])/2.0)*((V[i][j]-V[i+1][j])/2.0)-(fabs(U[i-1][j]+U[i-1][j+1])/2.0)*((V[i-1][j]-V[i][j])/2.0))
                )/dx
                -(((V[i][j]+V[i][j+1])/2.0)*((V[i][j]+V[i][j+1])/2.0)-((V[i][j-1]+V[i][j])/2.0)*((V[i][j-1]+V[i][j])/2.0)
                    +alpha*((fabs(V[i][j]+V[i][j+1])/2.0)*((V[i][j]-V[i][j+1])/2.0)-(fabs(V[i][j-1]+V[i][j])/2.0)*((V[i][j-1]-V[i][j])/2.0))
                )/dy
                +GY
            );
        }
    }
    if (il==1){
        for (int j=jb;j<=jt; j++){
            F[0][j]=U[0][j];
        }
    }
    if (ir==imax){
        for (int j=jb;j<=jt; j++){
            F[imax][j]=U[imax][j];
        }
    }
    if (jb==1){
        for (int i=il;i<=ir; i++){
            G[i][0]=V[i][0];
        }    
    }
    if (jt==jmax){
        for (int i=il;i<=ir; i++){
            G[i][jmax]=V[i][jmax];
        }
    }
}


/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  double **F,
  double **G,
  double **RS,
  int il,
  int ir,
  int jb,
  int jt
)
{
    for (int j=jb;j<=jt;j++){
        for (int i=il;i<=ir;i++){
            RS[i][j]=((F[i][j]-F[i-1][j])/dx+(G[i][j]-G[i][j-1])/dy)/dt;
        }
    }
}


/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  double *umax,
  double *vmax)
{
    *dt=tau*(fmin(Re/(2.0*(1.0/(dx*dx)+1.0/(dy*dy))),fmin(dx/(*umax),dy/(*vmax))));
}

void calculate_local_uv_max(
  double **U,
  double **V,
  double *umax,
  double *vmax,
  int il,
  int ir,
  int jb,
  int jt
){
    *umax=0;
    *vmax=0;
    for(int j = jb; j <= jt; j++) {
        for(int i = il; i<=ir; i++) {
            if (fabs(U[i][j])>*umax){
                *umax=fabs(U[i][j]);
            }
            if (fabs(V[i][j])>*vmax){
                *vmax=fabs(V[i][j]);
            }
        }
    }
}


/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */
void calculate_uv(
  double dt,
  double dx,
  double dy,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P,
  int il,
  int ir,
  int jb,
  int jt
)
{
    for (int j=jb;j<=jt;j++){
        for (int i=il;i<ir;i++){
            U[i][j]=F[i][j]-dt*(P[i+1][j]-P[i][j])/dx;
        }
    }
    for (int j=jb;j<jt;j++){
        for (int i=il;i<=ir;i++){
            V[i][j]=G[i][j]-dt*(P[i][j+1]-P[i][j])/dy;
        }
    }
}


