#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>
#include "precice_adapter.h"
//please include your path to the references correctly, like "Your-path-to-precice/src/precice/adapters/c/SolverInterfaceC.h"
#include "adapters/c/SolverInterfaceC.h"
#include "adapters/c/Constants.h"



int main(int argn, char** args){
  
  //The initial variables read from the .dat files, together with their pointers

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

  double valINUI=0.0;
  double valINVI=0.0;
  double valINTI=0.0;
  double *INUI=&valINUI;
  double *INVI=&valINVI;
  double *INTI=&valINTI;
  //The equivalent default value for a non Heat problem so that they can use the same equations
  double valTI=0.0;
  double valBeta=0.0;
  double valPr=1.0;
  double *TI=&valTI;
  double *beta=&valBeta;
  double *Pr=&valPr;

  double valx_origin=0.0;
  double valy_origin=0.0;
  double *x_origin=&valx_origin;
  double *y_origin=&valy_origin;

  char* precice_config = (char*) malloc(sizeof(char) * 100); //the path to the precice-config.xml file,
  char* participant_name= (char*) malloc(sizeof(char) * 20); //which should typically be Fluid,
  char* mesh_name= (char*) malloc(sizeof(char) * 20); //which should typically be Fluid-Mesh,
  char* read_data_name= (char*) malloc(sizeof(char) * 20); //which should typically be Heat-Flux, and
  char* write_data_name= (char*) malloc(sizeof(char) * 20); //which should typically be Temperature.


  //Set flags for the question: is this a fluid or heat problem? is it karman or step or natural convection or heat trap or rbc?

  char* DAT;//address of the .dat file
  char* problem= (char*) malloc(sizeof(char) * 50);//type of problem: 1.heated-plate, 2.convection, 3.F1-heat-exchange, 4.F2-heat-exchange
  char* geometry= (char*) malloc(sizeof(char) * 50);//address of the .pgm file
  char* Outputfilename;//where to save output

  double Th=323.0;
  double Tc=283.0;

  int NameQuestion=0;//1=plate, 2=convection, 3=exchange1, 4=exchange2

  if (strcmp(args[1],"plate")==0){
    DAT="configs/heated-plate.dat";
    Outputfilename="output/heated-plate";
    NameQuestion=1;
  }
  else if (strcmp(args[1],"convection")==0){
    DAT="configs/convection.dat";
    Outputfilename="output/convection";
    NameQuestion=2;
  }
  else if (strcmp(args[1],"exchange1")==0){
    DAT="configs/F1-heat-exchange.dat";
    Outputfilename="output/F1-heat-exchange";
    NameQuestion=3;
  }
  else if (strcmp(args[1],"exchange2")==0){
    DAT="configs/F2-heat-exchange.dat";
    Outputfilename="output/F2-heat-exchange";
    NameQuestion=4;
  }
  

  //read parameters
  read_parameters(DAT, Re,UI,VI,PI,GX,GY,t_end,xlength,ylength,dt,dx,dy,imax,jmax,alpha,omg,tau,itermax,eps,dt_value, TI, Pr, beta, x_origin, y_origin, problem, geometry, precice_config, participant_name, mesh_name, read_data_name, write_data_name);

  valINUI=valUI;
  valINVI=valVI;
  valINTI=valTI;

  //reserve memory for the global matrices

  double **U = matrix ( 0 , *imax , 0 , *jmax+1 );
  double **V = matrix ( 0 , *imax+1 , 0 , *jmax );
  double **P = matrix ( 0 , *imax+1 , 0 , *jmax+1 );
  double **RS = matrix ( 1 , *imax , 1 , *jmax );
  double **F = matrix ( 0 , *imax , 1 , *jmax );
  double **G = matrix ( 1 , *imax , 0 , *jmax );
  double **Temp = matrix (0, *imax+1, 0, *jmax+1);

  double **U_cp = matrix ( 0 , *imax , 0 , *jmax+1 );
  double **V_cp = matrix ( 0 , *imax+1 , 0 , *jmax );
  double **Temp_cp = matrix (0, *imax+1, 0, *jmax+1);
  double valt_cp=0.0;
  double *t_cp=&valt_cp;

  init_matrix(F, 0, *imax , 1 , *jmax, 0.0);
  init_matrix(G, 1, *imax , 0 , *jmax, 0.0);
  init_matrix(RS, 1, *imax , 1, *jmax, 0.0);
  init_matrix(U_cp, 0, *imax , 0 , *jmax+1, 0.0);
  init_matrix(V_cp, 0, *imax+1 , 0 , *jmax, 0.0);
  init_matrix(Temp_cp, 0, *imax+1 , 0, *jmax+1, 0.0);

  int **Flag=imatrix(0,*imax+1,0,*jmax+1);
  init_imatrix(Flag, 0,*imax+1,0,*jmax+1, 0);
  double **Geo=matrix(0,*imax+1,0,*jmax+1);
  init_matrix(Geo, 0,*imax+1,0,*jmax+1, 0.0);

  //set global variables  to compute timestamps, residual of sor, and restTime variable judges when to print out a vtk file

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

  //initialize U,V,P, and T

  init_uvp(*UI,*VI,*PI,*imax,*jmax,U,V,P, *TI, Temp);

  //initialize flags together with its geometry structure (Geo, which can be print out later)

  init_flag(problem,geometry, *imax, *jmax, Flag, Geo);
  
  //Compute the number of fluid cells in the problem, can be used in the sor iteration as a global variable

  double NumberFluidCell=0.0;
  int num_coupling_cells=0;
  for (int i=0;i<=(*imax+1);i++){
    for (int j=0;j<=(*jmax+1);j++){
      if (Flag[i][j]>511){
        NumberFluidCell+=1.0;
      }
      if ((Flag[i][j]&256)==256){
        num_coupling_cells+=1;
      }
    }
  }
  printf("num_coupling_cell=%i\n", num_coupling_cells);

  
  
  // initialize preCICE
  precicec_createSolverInterface(participant_name, precice_config, 0, 1);

  const char* coric = precicec_actionReadIterationCheckpoint();
  const char* cowic = precicec_actionWriteIterationCheckpoint();

  int dim = precicec_getDimensions();
  //printf("dim= %i\n",dim);
  // define coupling mesh
  int meshID = precicec_getMeshID(mesh_name);
  //printf("meshID= %i\n",meshID);

  //int num_coupling_cells = ... // determine number of coupling cells

  int* vertexIDs = precice_set_interface_vertices(*imax, *jmax, *dx, *dy, *x_origin, *y_origin, num_coupling_cells, meshID, Flag);// get coupling cell ids
  //printf("vertexIDs[0:2]= %i, %i, %i\n",vertexIDs[0],vertexIDs[1],vertexIDs[2]);

  // define Dirichlet part of coupling written by this solver
  int temperatureID = precicec_getDataID(write_data_name, meshID);
  //printf("temperatureID= %i\n",temperatureID);
  double* temperatureCoupled = (double*) malloc(sizeof(double) * num_coupling_cells);

  // define Neumann part of coupling read by this solver
  int heatFluxID = precicec_getDataID(read_data_name, meshID);
  //printf("heatfluxID= %i\n",heatFluxID);
  double* heatfluxCoupled = (double*) malloc(sizeof(double) * num_coupling_cells);
  // call precicec_initialize()
  double precice_dt=precicec_initialize();

  
  // initialize data at coupling interface
  precice_write_temperature(*imax, *jmax, num_coupling_cells, temperatureCoupled, vertexIDs, temperatureID, Temp, Flag);
  precicec_initialize_data(); // synchronize with OpenFOAM
  precicec_readBlockScalarData(heatFluxID, num_coupling_cells, vertexIDs, heatfluxCoupled); // read heatfluxCoupled

  //write a vtk file at the beginning

  boundaryvalues(*imax, *jmax, U, V, Flag, NameQuestion, Th, Tc, Temp, *INUI, *INVI, *INTI);
  set_coupling_boundary(*imax, *jmax, *dx, *dy, heatfluxCoupled, Temp, Flag);
  write_vtkFile(Outputfilename, *n, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P, Geo, Temp, *x_origin, *y_origin);
  

  while (precicec_isCouplingOngoing()) { // time loop

    if (precicec_isActionRequired(cowic)){
      printf("DUMMY: Writing iteration checkpoint %i\n", 0);
      write_checkpoint(*t, U, V, Temp, t_cp, U_cp, V_cp, Temp_cp, *imax, *jmax);
      precicec_fulfilledAction(cowic);
    }

    //1. calculate time step
    // use dt = min(solver_dt, precice_dt)
    if (*tau<=0){
      *dt=*dt_value;//0.005;
    }
    else{
      calculate_dt(*Re, *tau, dt, *dx, *dy, *imax, *jmax, U, V, Flag,*Pr);
    }
    *dt=fmin(*dt, precice_dt);

    //2. set boundary values
    //compute boundary for U,V,T
    boundaryvalues(*imax, *jmax, U, V, Flag, NameQuestion, Th, Tc, Temp, *INUI, *INVI, *INTI);
    set_coupling_boundary(*imax, *jmax, *dx, *dy, heatfluxCoupled, Temp, Flag);

    printf("heatfluxCoupled[0]=%f\n", heatfluxCoupled[0]);
    printf("temperature[0]=%f\n", temperatureCoupled[0]);
    
    //3 - 6. calculate temp, F and G | RHS of P eq. | pressure | new U,V

    
    
    //compute T^(n+1) if this is a heat question

    calculate_temp(*Re, *Pr, *dt, *dx, *dy, *imax, *jmax, U, V, Temp, *alpha, Flag, *TI);
    

    //compute F, G

    calculate_fg(*Re, *GX, *GY, *alpha, *dt, *dx, *dy, *imax, *jmax, U, V, F, G, Flag,*beta, Temp);

    //compute right hand side

    calculate_rs(*dt, *dx, *dy, *imax, *jmax, F, G, RS, Flag);

    //reinitialize the counter of iteration for the sor solver

    *it=0;
    *res=*eps+1;

    //compute sor, either converge or not converge but arrive at itermax

    while ((*it<*itermax) && (*res>*eps)){
      sor(*omg, *dx, *dy, *imax, *jmax, P, RS, res, Flag,NumberFluidCell);
      *it=*it+1;
      solbeh=solbeh+1;
    }
    if (*it==*itermax){
      printf("WARNING! The SOR have not converged in this iteration! %i\n", 0);
    }

    //compute new U, V

    calculate_uv(*dt, *dx, *dy, *imax, *jmax, U, V, F, G, P, Flag);
    
    //7. coupling
    precice_write_temperature(*imax, *jmax, num_coupling_cells, temperatureCoupled, vertexIDs, temperatureID, Temp, Flag); // write new temperature to preCICE buffers
    
    precice_dt = precicec_advance(*dt);

    // advance coupling
    precicec_readBlockScalarData(heatFluxID, num_coupling_cells, vertexIDs, heatfluxCoupled);// read new heatflux from preCICE buffers

    if(precicec_isActionRequired(coric)){ // timestep not converged
		  restore_checkpoint(t, U, V, Temp, *t_cp, U_cp, V_cp, Temp_cp, *imax, *jmax);
      printf("DUMMY: Reading iteration checkpoint %i\n", 0);
		  precicec_fulfilledAction(coric);
    }
    else{
      //8. output U, V, P for visualization and update iteration values

      *t+=*dt;
      restTime+=*dt; 
      if (restTime>=*dt_value){
        boundaryvalues(*imax, *jmax, U, V, Flag, NameQuestion, Th, Tc, Temp, *INUI, *INVI, *INTI);
        set_coupling_boundary(*imax, *jmax, *dx, *dy, heatfluxCoupled, Temp, Flag);
        write_vtkFile(Outputfilename, *n, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P, Geo, Temp, *x_origin, *y_origin);
        restTime=restTime-*dt_value;
      }
      *n=*n+1;
      

      //print out n, dt, t, #sor, Residual in the terminal
      printf("n=%i",*n);
      printf(", dt=%f",*dt);
      printf(", t=%f ,\n   ",*t);
      printf("#SOR=%i",*it);
      printf(", Res=%f\n",*res);
    }
  }
  precicec_finalize();

  boundaryvalues(*imax, *jmax, U, V, Flag, NameQuestion, Th, Tc, Temp, *INUI, *INVI, *INTI);
  set_coupling_boundary(*imax, *jmax, *dx, *dy, heatfluxCoupled, Temp, Flag);
  write_vtkFile(Outputfilename, *n, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P, Geo, Temp, *x_origin, *y_origin);

  
  //print how many iterations in sor in all 

  printf("Total #SOR=%i \n", solbeh);

  //free memory space

  free_matrix (U, 0 , *imax , 0 , *jmax+1 );
  free_matrix (V, 0 , *imax+1 , 0 , *jmax );
  free_matrix (P, 0 , *imax+1 , 0 , *jmax+1 );
  free_matrix (RS, 1 , *imax , 1 , *jmax );
  free_matrix (F, 0 , *imax , 1 , *jmax );
  free_matrix (G, 1 , *imax , 0 , *jmax );
  free_imatrix(Flag, 0, *imax+1, 0, *jmax+1);
  free_matrix(Temp, 0, *imax+1, 0, *jmax+1);
  return 0;
}
