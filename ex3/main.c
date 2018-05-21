#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include "mpi.h"
#include "parallel.h"
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

  //Definit variables and read them from file
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

  int valiproc=1;
  int valjproc=1;
  int *iproc=&valiproc;
  int *jproc=&valjproc;

  read_parameters("cavity100.dat", Re,UI,VI,PI,GX,GY,t_end,xlength,ylength,dt,dx,dy,imax,jmax,alpha,omg,tau,itermax,eps,dt_value,iproc,jproc);
  printf("-------------------------------------------------------------\n");

  int valil=1;
  int valir=1;
  int valjb=1;
  int valjt=1;
  int valrank_l=0;
  int valrank_r=0;
  int valrank_b=0;
  int valrank_t=0;
  int valomg_i=1;
  int valomg_j=1;
  int *il=&valil;
  int *ir=&valir;
  int *jb=&valjb;
  int *jt=&valjt;
  int *rank_l=&valrank_l;
  int *rank_r=&valrank_r;
  int *rank_b=&valrank_b;
  int *rank_t=&valrank_t;
  int *omg_i=&valomg_i;
  int *omg_j=&valomg_j;
  
  int myrank;
  int num_proc;
  MPI_Status status;

  //Initialization
  MPI_Init( &argn, &args );
  MPI_Comm_size( MPI_COMM_WORLD, &num_proc );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  init_parallel(*iproc, *jproc, *imax, *jmax, &myrank, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, omg_i, omg_j, &num_proc);

  double **U = matrix ( *il-2 , *ir+1 , *jb-1 , *jt+1 );
  double **V = matrix ( *il-1 , *ir+1 , *jb-2 , *jt+1 );
  double **P = matrix ( *il-1 , *ir+1 , *jb-1 , *jt+1 );
  double **RS = matrix ( *il , *ir , *jb , *jt);
  double **F = matrix ( *il-2 , *ir+1 , *jb-1 , *jt+1 );
  double **G = matrix ( *il-1 , *ir+1 , *jb-2 , *jt+1 );

  init_matrix(F, *il-2 , *ir+1 , *jb-1 , *jt+1 , 0.0);
  init_matrix(G, *il-1 , *ir+1 , *jb-2 , *jt+1, 0.0);
  init_matrix(RS, *il , *ir , *jb , *jt, 0.0);

  double valt=0;
  int valn=0;
  int valit=0;
  double valres= valeps+1;
  double valreslocal=valeps+1;
  double restTime=0;
  double valumax=0;
  double valumaxlocal=0;
  double valvmax=0;
  double valvmaxlocal=0;


  double *t=&valt;
  int *n=&valn;
  int *it=&valit;
  double *res= &valres;
  double *reslocal=&valreslocal;
  double *umax=&valumax;
  double *umaxlocal=&valumaxlocal;
  double *vmax=&valvmax;
  double *vmaxlocal=&valvmaxlocal;
  int chunk=2*(max((*imax)/(*iproc),(*jmax)/(*jproc))+1)+1;

  double bufSend[chunk];
  double bufRecv[chunk];

  //to test function MPI_Allreduce
  /*double a=myrank*0.1;
  printf("a %i, %f\n", myrank, a);
  double b;
  MPI_Allreduce(&a, &b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  printf("b %i, %f\n", myrank, b);*/


  //int solbeh=0;
  init_uvp(*UI,*VI,*PI,U,V,P,*il, *ir, *jb, *jt);
  boundaryvalues(*imax, *jmax, U, V, *il, *ir, *jb, *jt);
  
  while (*t<*t_end){

    //compute dt
    MPI_Barrier(MPI_COMM_WORLD); 
    if (*tau<=0){
      *dt=*dt_value;//0.005;
    }
    else{
      //Compute the maximum values of u (n+1) and v (n+1) for each process
      calculate_local_uv_max(U,V, umaxlocal, vmaxlocal, *il, *ir, *jb, *jt);
      MPI_Barrier(MPI_COMM_WORLD); 
      //computes the global maximum values of u (n+1) and v (n+1) and broadcasts them to all other processes
      MPI_Allreduce(umaxlocal, umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(vmaxlocal, vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      //Compute the common dt according to Worksheet 1
      calculate_dt(*Re, *tau, dt, *dx, *dy, umax, vmax);
    }

    //Set boundary values for u und v
    boundaryvalues(*imax, *jmax, U, V, *il, *ir, *jb, *jt);

    //Compute F (n) and G (n) according to Worksheet 1
    calculate_fg(*Re, *GX, *GY, *alpha, *dt, *dx, *dy, *imax, *jmax, U, V, F, G, *il, *ir, *jb, *jt);
    MPI_Barrier(MPI_COMM_WORLD);
    //uv_comm(F,G,*il,*ir,*jb,*jt,*rank_l,*rank_r,*rank_b,*rank_t,bufSend,bufRecv,&status,chunk);
    //MPI_Barrier(MPI_COMM_WORLD);

    //Compute the right-hand side of the pressure equation (Worksheet 1)
    calculate_rs(*dt, *dx, *dy, F, G, RS, *il, *ir, *jb, *jt);

    *it=0;
    *res=*eps+1;
    MPI_Barrier(MPI_COMM_WORLD);

    while ((*it<*itermax) && (sqrt(*res)>*eps)){

      //Set the appropriate pressure boundary values (Worksheet 1) and perform an SOR cycle according to Worksheet 1
      MPI_Barrier(MPI_COMM_WORLD);
      sor(*omg, *dx, *dy, *imax, *jmax, P, RS, reslocal, *il, *ir, *jb, *jt);
      MPI_Barrier(MPI_COMM_WORLD);

      //Exchange the pressure values in the boundary strips
      pressure_comm(P, *il, *ir, *jb, *jt, *rank_l, *rank_r, *rank_b, *rank_t, bufSend, bufRecv, &status, chunk);
      *it=*it+1;

      //Compute the partial residual sum and send it to the master process
      //The master process computes the residual norm of the pressure equation and broadcasts it to all the other processes
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(reslocal, res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      //solbeh=solbeh+1;
    }
    if (*it==*itermax){
      if (myrank==0){
        printf("WARNING! The SOR have not converged in this iteration!\n");
      }
    }
    //Compute u (n+1) and v (n+1) according to Worksheet 1
    calculate_uv(*dt, *dx, *dy, U, V, F, G, P, *il,*ir,*jb,*jt);
    //Exchange the velocity values in the boundary strips
    MPI_Barrier(MPI_COMM_WORLD);
    uv_comm(U,V,*il,*ir,*jb,*jt,*rank_l,*rank_r,*rank_b,*rank_t,bufSend,bufRecv,&status,chunk);
    boundaryvalues(*imax, *jmax, U, V, *il, *ir, *jb, *jt);
    
    *t+=*dt;
    restTime+=*dt; 
    *n=*n+1;
    /*if (restTime>=*dt_value){
      write_vtkFile("output/output", *n, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P);
      restTime=restTime-*dt_value;
    }*/
    if(myrank==0){
      printf("n=");
      printf("%i",*n);
      printf(", dt=");
      printf("%f",*dt);
      printf(", t=");
      printf("%f ,\n   ",*t);
      printf("#SOR=");
      printf("%i",*it);
      printf(", Res=");
      printf("%f\n",sqrt(*res));
      /*printf("   U=");
      printf("%f",U[2][2]);
      printf(", V=");
      printf("%f",V[2][2]);
      printf(", P=");
      printf("%f\n",P[2][2]);*/
      
    }
    
  }
  //printf("Total #SOR=");
  //printf("%i \n", solbeh);
  //output the result
  output_uvp(U, V, P, *il, *ir, *jb, *jt, *omg_i, *omg_j, "output/MPI", *dx, *dy);

  free_matrix (U, *il-2 , *ir+1 , *jb-1 , *jt+1  );
  free_matrix (V, *il-1 , *ir+1 , *jb-2 , *jt+1 );
  free_matrix (P, *il-1 , *ir+1 , *jb-1 , *jt+1);
  free_matrix (RS, *il , *ir , *jb , *jt );
  free_matrix (F, *il-2 , *ir+1 , *jb-1 , *jt+1);
  free_matrix (G, *il-1 , *ir+1 , *jb-2 , *jt+1);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
