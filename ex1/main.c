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
  double valRe=0;
  double valUI=0;
  double valVI=0;
  double valPI=0;
  double valGX=0;
  double valGY=0;
  double valt_end=0;
  double valxlength=0;
  double valylength=0;
  double valdt=0;
  double valdx=0;
  double valdy=0;
  int  valimax=0;
  int  valjmax=0;
  double valalpha=0;
  double valomg=0;
  double valtau=0;
  int  valitermax=0;
  double valeps=0;
	double valdt_value=0;

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
  //*imax=512;
  //*jmax=512;
  //*t_end=500;
  //Q7 for changing Re
  //*Re=10000;
  
  double **U = matrix ( 0 , *imax , 0 , *jmax+1 );
  double **V = matrix ( 0 , *imax+1 , 0 , *jmax );
  double **P = matrix ( 0 , *imax+1 , 0 , *jmax+1 );
  double **RS = matrix ( 1 , *imax , 1 , *jmax );
  double **F = matrix ( 0 , *imax , 1 , *jmax );
  double **G = matrix ( 1 , *imax , 0 , *jmax );

  init_matrix(F, 0, *imax , 1 , *jmax, 0.0);
  init_matrix(G, 1, *imax , 0 , *jmax, 0.0);
  init_matrix(RS, 1, *imax , 1, *jmax, 0.0);

  double valt=0;
  int valn=0;
  int valit=0;
  double valres= valeps+1;
  double restTime=0;

  double *t=&valt;
  int *n=&valn;
  int *it=&valit;
  double *res= &valres;
  int solbeh=0;
  init_uvp(*UI,*VI,*PI,*imax,*jmax,U,V,P);
  boundaryvalues(*imax, *jmax, U, V);
  //write_vtkFile("output/output", *n, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P);

  while (*t<*t_end){
    if (*tau<=0){
      *dt=*dt_value;//0.005;
    }
    else{
      calculate_dt(*Re, *tau, dt, *dx, *dy, *imax, *jmax, U, V);
    }
    boundaryvalues(*imax, *jmax, U, V);
    calculate_fg(*Re, *GX, *GY, *alpha, *dt, *dx, *dy, *imax, *jmax, U, V, F, G);
    calculate_rs(*dt, *dx, *dy, *imax, *jmax, F, G, RS);
    *it=0;
    *res=*eps+1;
    while ((*it<*itermax) && (*res>*eps)){
      sor(*omg, *dx, *dy, *imax, *jmax, P, RS, res);
      *it=*it+1;
      solbeh=solbeh+1;
    }
    if (*it==*itermax){
      printf("WARNING! The SOR have not converged in this iteration!\n");
    }
    calculate_uv(*dt, *dx, *dy, *imax, *jmax, U, V, F, G, P);
    //Q5 help check if stable
    /*for (int i=1;i<=*imax;i++){
      for (int j=1;j<=*jmax;j++){
        if ((isnan(U[i][j]))||(isnan(V[i][j]))||(isnan(F[i][j]))||(isnan(G[i][j]))||(isnan(P[i][j]))||(isnan(RS[i][j]))){
          printf("Unstable, exist NaN in U,V,P,F,G or RS \n");
          free_matrix (U, 0 , *imax , 0 , *jmax+1 );
          free_matrix (V, 0 , *imax+1 , 0 , *jmax );
          free_matrix (P, 0 , *imax+1 , 0 , *jmax+1 );
          free_matrix (RS, 1 , *imax , 1 , *jmax );
          free_matrix (F, 0 , *imax , 1 , *jmax );
          free_matrix (G, 1 , *imax , 0 , *jmax );
          return 0;
        }
      }
    }*/
    *t+=*dt;
    restTime+=*dt; 
    *n=*n+1;
    if (restTime>=*dt_value){
      //write_vtkFile("output/output", *n, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P);
      restTime=restTime-*dt_value;
    }
    printf("n=");
    printf("%i",*n);
    printf(", dt=");
    printf("%f",*dt);
    printf(", t=");
    printf("%f ,\n   ",*t);
    printf("#SOR=");
    printf("%i",*it);
    printf(", Res=");
    printf("%f",*res);
    printf(", U[imax/2][7*jmax/8]=");
    printf("%f \n",U[*imax/2][(int)(((double)*jmax)/8.0*7.0)]);
    //stop if is NaN
    if (isnan(U[*imax/2][(int)(((double)*jmax)/8.0*7.0)])){
      free_matrix (U, 0 , *imax , 0 , *jmax+1 );
      free_matrix (V, 0 , *imax+1 , 0 , *jmax );
      free_matrix (P, 0 , *imax+1 , 0 , *jmax+1 );
      free_matrix (RS, 1 , *imax , 1 , *jmax );
      free_matrix (F, 0 , *imax , 1 , *jmax );
      free_matrix (G, 1 , *imax , 0 , *jmax ); 
      return 0;
    }
  }
  printf("Total #SOR=");
  printf("%i \n", solbeh);
  free_matrix (U, 0 , *imax , 0 , *jmax+1 );
  free_matrix (V, 0 , *imax+1 , 0 , *jmax );
  free_matrix (P, 0 , *imax+1 , 0 , *jmax+1 );
  free_matrix (RS, 1 , *imax , 1 , *jmax );
  free_matrix (F, 0 , *imax , 1 , *jmax );
  free_matrix (G, 1 , *imax , 0 , *jmax );
  return 0;
}
