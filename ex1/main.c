#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){
  double valRe=1;
  double valUI=1;
  double valVI=1;
  double valPI=1;
  double valGX=1;
  double valGY=1;
  double valt_end=1;
  double valxlength=1;
  double valylength=1;
  double valdt=1;
  double valdx=1;
  double valdy=1;
  int  valimax=1;
  int  valjmax=1;
  double valalpha=1;
  double valomg=1;
  double valtau=1;
  int  valitermax=1;
  double valeps=1;
	double valdt_value=1;

  double *Re=&valRe;
  double *UI=&valUI;
  double *VI=&valVI;
  double *PI=&valPI;
  double *GX=&valGX;
  double *GY=&valGY;
  double *t_end=&valt_end;
  double *xlength=&valxlength;
  double *ylength=&valylength;
  double *dt=&valdt;
  double *dx=&valdx;
  double *dy=&valdy;
  int  *imax=&valimax;
  int  *jmax=&valjmax;
  double *alpha=&valalpha;
  double *omg=&valomg;
  double *tau=&valtau;
  int  *itermax=&valitermax;
  double *eps=&valeps;
	double *dt_value=&valdt_value;
  read_parameters("cavity100.dat", Re,UI,VI,PI,GX,GY,t_end,xlength,ylength,dt,dx,dy,imax,jmax,alpha,omg,tau,itermax,eps,dt_value);

  //Q4 for changing omega
  //*omg=2.0;
  //Q6 for U tends to infinity
  *imax=512;
  *jmax=512;
  *t_end=100;
  //Q7 for changing Re
  //*Re=10000;
  
  double **U = matrix ( 0 , *imax+1 , 0 , *jmax+1 );
  double **V = matrix ( 0 , *imax+1 , 0 , *jmax+1 );
  double **P = matrix ( 0 , *imax+1 , 0 , *jmax+1 );
  double **RS = matrix ( 0 , *imax+1 , 0 , *jmax+1 );
  double **F = matrix ( 0 , *imax+1 , 0 , *jmax+1 );
  double **G = matrix ( 0 , *imax+1 , 0 , *jmax+1 );

  init_matrix(F,0 , *imax+1 , 0 , *jmax+1, 0.0);
  init_matrix(G,0 , *imax+1 , 0 , *jmax+1, 0.0);
  init_matrix(RS,0 , *imax+1 , 0 , *jmax+1, 0.0);

  double valt=0;
  int valn=0;
  int valit=0;
  double valres= valeps+1;

  double *t=&valt;
  int *n=&valn;
  int *it=&valit;
  double *res= &valres;
  init_uvp(*UI,*VI,*PI,*imax,*jmax,U,V,P);

  while (*t<*t_end){
    if (*tau<=0){//(1>0){
      *dt=*dt_value;//0.005;
    }
    else{
      calculate_dt(*Re, *tau, dt, *dx, *dy, *imax, *jmax, U, V);
    }
    boundaryvalues(*imax, *jmax, U, V);
    calculate_fg(*Re, *GX, *GY, *alpha, *dt, *dx, *dy, *imax, *jmax, U, V, F, G);
    calculate_rs(*dt, *dx, *dy, *imax, *jmax, F, G, RS);
    *it=0;
    *res=valeps+1;
    while ((*it<*itermax) && (*res>*eps)){
      sor(*omg, *dx, *dy, *imax, *jmax, P, RS, res);
      *it=*it+1;
    }
    calculate_uv(*dt, *dx, *dy, *imax, *jmax, U, V, F, G, P);
    *t+=*dt; 
    *n=*n+1;
    if (1>0){
      //write_vtkFile("output170/output", *n, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P);
      printf("n=");
      printf("%i",*n);
      printf(", t=");
      printf("%f",*t);
      printf(", U[imax/2][7*jmax/8]=");
      printf("%f \n",U[*imax/2][*jmax/8*7]);
    }
  }
  free_matrix (U, 0 , *imax+1 , 0 , *jmax+1 );
  free_matrix (V, 0 , *imax+1 , 0 , *jmax+1 );
  free_matrix (P, 0 , *imax+1 , 0 , *jmax+1 );
  free_matrix (RS, 0 , *imax+1 , 0 , *jmax+1 );
  free_matrix (F, 0 , *imax+1 , 0 , *jmax+1 );
  free_matrix (G, 0 , *imax+1 , 0 , *jmax+1 );
  return 0;
}
