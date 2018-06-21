#include "helper.h"
#include "init.h"

int read_parameters( const char *szFileName,       /* name of the file */
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                               /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
		    double *dt_value,           /* time for output */
        double *TI,
        double *Pr,
        double *beta,
        double *x_origin,
        double *y_origin,
        char* problem,
        char* geometry,
        char* precice_config, //the path to the precice-config.xml file,
        char* participant_name, //which should typically be Fluid,
        char* mesh_name, //which should typically be Fluid-Mesh,
        char* read_data_name, //which should typically be Heat-Flux, and
        char* write_data_name //which should typically be Temperature.

        )
{
   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );

   READ_INT   ( szFileName, *itermax );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   READ_DOUBLE( szFileName, *TI );
   READ_DOUBLE( szFileName, *Pr );
   READ_DOUBLE( szFileName, *beta );

   //READ_DOUBLE( szFileName, *INUI );
   //READ_DOUBLE( szFileName, *INVI );
   //READ_DOUBLE( szFileName, *INTI );

   READ_DOUBLE( szFileName, *x_origin);
   READ_DOUBLE( szFileName, *y_origin);

   READ_STRING( szFileName, problem);
   READ_STRING( szFileName, geometry);

   READ_STRING( szFileName, precice_config);
   READ_STRING( szFileName, participant_name);
   READ_STRING( szFileName, mesh_name);
   READ_STRING( szFileName, read_data_name);
   READ_STRING( szFileName, write_data_name);

   return 1;
}

//initialize the U, V, P, T matrices
void init_uvp(
  double UI,
  double VI,
  double PI,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **P,
  double TI,
  double **T)
{
    init_matrix(U,0,imax,0,jmax+1, UI);
    init_matrix(V,0,imax+1,0,jmax, VI);
    init_matrix(P,0,imax+1,0,jmax+1, PI);
    init_matrix(T,0,imax+1,0,jmax+1, TI);
}

//initialize Flags. each Flag case is a int type between 0 and 512+15=527. Its bit representation represent:
//fluid/ coupling/ no-slip/ free-slip/ outflow/ inflow/ B_N/ B_S/ B_W/ B_E
void init_flag(
  char* problem,
  char* geometry,
  int imax,
  int jmax,
  int **Flag,
  double **Geo
){
  int **aux=imatrix(0,imax+1,0,jmax+1);
  aux=read_pgm(geometry);
  for (int j=0;j<=jmax+1;j++){
    for (int i=0;i<=imax+1;i++){
      Geo[i][j]=(double)aux[i][j];
      Flag[i][j]=0;
      if (aux[i][j]==6){//fluid
        Flag[i][j]+=512;
      }
      if (aux[i][j]==4){//coupling
        Flag[i][j]+=256;
      }
      if (aux[i][j]==0){//no-slip
        Flag[i][j]+=128;
      }
      if (aux[i][j]==1){//free-slip
        Flag[i][j]+=64;
      }
      if (aux[i][j]==2){//outflow
        Flag[i][j]+=32;
      }
      if (aux[i][j]==3){//inflow
        Flag[i][j]+=16;
      }
      if ((j<(jmax+1))&&(aux[i][j+1]==6)){//BN
        Flag[i][j]+=8;
      }
      if ((j>0)&&(aux[i][j-1]==6)){//BS
        Flag[i][j]+=4;
      }
      if ((i>0)&&(aux[i-1][j]==6)){//BW
        Flag[i][j]+=2;
      }
      if ((i<(imax+1))&&(aux[i+1][j]==6)){//BE
        Flag[i][j]+=1;
      }
      //assert only one field from 5 cases. This assertion is useless according to how we build Flag[i][j]
      //assert((Flag[i][j]&240==128)||(Flag[i][j]&240==64)||(Flag[i][j]&240==32)||(Flag[i][j]&240==16)||(Flag[i][j]&240==0));//&& strcat(strcat(strcat("Error! Fields set unproperly at point ",itoa(i)),", "),itoa(j)));
      //assert fluid or obstacle without 2 opposite fluid cells
      assert((Flag[i][j]>511)||((Flag[i][j]&12)!=12));//NS not same time
      assert((Flag[i][j]>511)||((Flag[i][j]&3)!=3));//WE not same time
      //The above 2 assetions can also make sure that we have not an obstacle cell with more than 2 adjacent fluid cells
    }
  }
  free_imatrix(aux, 0, imax+1, 0, jmax+1);
}
