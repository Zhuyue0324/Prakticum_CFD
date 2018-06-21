#include "precice_adapter.h"
#include <stdlib.h>
#include <stdio.h>
#include "adapters/c/SolverInterfaceC.h"

void write_checkpoint(double time, double **U, double **V, double **TEMP, double *time_cp, double **U_cp, double **V_cp,
                      double **TEMP_cp, int imax, int jmax)
    {
        *time_cp=time;
        for (int j=0;j<=(jmax+1);j++){
            for (int i=0;i<=imax;i++){
                U_cp[i][j]=U[i][j];
            }
        }
        for (int j=0;j<=jmax;j++){
            for (int i=0;i<=(imax+1);i++){
                V_cp[i][j]=V[i][j];
            }
        }
        for (int j=0;j<=(jmax+1);j++){
            for (int i=0;i<=(imax+1);i++){
                TEMP_cp[i][j]=TEMP[i][j];
            }
        }
    }

void restore_checkpoint(double *time, double **U, double **V, double **TEMP, double time_cp, double **U_cp,
                        double **V_cp,
                        double **TEMP_cp, int imax, int jmax)
    {
        *time=time_cp;
        for (int j=0;j<=(jmax+1);j++){
            for (int i=0;i<=imax;i++){
                U[i][j]=U_cp[i][j];
            }
        }
        for (int j=0;j<=jmax;j++){
            for (int i=0;i<=(imax+1);i++){
                V[i][j]=V_cp[i][j];
            }
        }
        for (int j=0;j<=(jmax+1);j++){
            for (int i=0;i<=(imax+1);i++){
                TEMP[i][j]=TEMP_cp[i][j];
            }
        }
    }


int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin,
                                    int num_coupling_cells, int meshID, int **FLAG)
    {
        double* vertices = (double*) malloc(sizeof(double) * num_coupling_cells *3);
        int* vertexIDs = (int*) malloc(sizeof(int) * num_coupling_cells);

        //double vertices[num_coupling_cells*3];
        //int vertexIDs[num_coupling_cells];
        for (int i=0;i<(3*num_coupling_cells);i++){
            vertices[i]=0.0;
        }
        int idx=0;
        for (int j=1;j<=jmax;j++){
            if ((FLAG[0][j]&256)==256){
                vertices[idx]=x_origin;
                vertices[idx+1]=y_origin +((double)j-0.5)*dy;
                vertices[idx+2]=0.0;
                idx=idx+3;
            }
        }
        for (int j=1;j<=jmax;j++){
            if ((FLAG[imax+1][j]&256)==256){
                vertices[idx]=x_origin+((double)imax)*dx;
                vertices[idx+1]=y_origin +((double)j-0.5)*dy;
                vertices[idx+2]=0.0;
                idx=idx+3;
            }
        }
        for (int i=1;i<=imax;i++){
            if ((FLAG[i][jmax+1]&256)==256){
                vertices[idx]=x_origin +((double)i-0.5)*dx;
                vertices[idx+1]=y_origin+((double)jmax)*dy;
                vertices[idx+2]=0.0;
                idx=idx+3;
            }
        }
        for (int i=1;i<=imax;i++){
            if ((FLAG[i][0]&256)==256){
                vertices[idx]=x_origin +((double)i-0.5)*dx;
                vertices[idx+1]=y_origin;
                vertices[idx+2]=0.0;
                idx=idx+3;
            }
        }
        for (int j=1;j<=jmax;j++){
            for (int i=1;i<=imax;i++){
                if ((FLAG[i][j]&264)==264){//coupling + north fluid = 256 + 8
                    vertices[idx]=x_origin +((double)i-0.5)*dx;
                    vertices[idx+1]=y_origin +((double)j)*dy;
                    vertices[idx+2]=0.0;
                    idx=idx+3;
                }
                if ((FLAG[i][j]&260)==260){//coupling + south fluid = 256 + 4
                    vertices[idx]=x_origin +((double)i-0.5)*dx;
                    vertices[idx+1]=y_origin +((double)j-1.0)*dy;
                    vertices[idx+2]=0.0;
                    idx=idx+3;
                }
            }
        }
        precicec_setMeshVertices(meshID, num_coupling_cells, vertices, vertexIDs);
        //free(vertices);
        return vertexIDs;
    }


void precice_write_temperature(int imax, int jmax, int num_coupling_cells, double *temperature, int *vertexIDs,
                               int temperatureID, double **TEMP, int **FLAG)
    {
        //for (int i=0;i<num_coupling_cells;i++){
        //    temperature[i]=0.0;
        //}
        int idx=0;
        for (int j=1;j<=jmax;j++){
            if (((FLAG[0][j]&256)==256)&&(FLAG[1][j]>511)){
                temperature[idx]=TEMP[1][j];
                idx=idx+1;
            }
        }
        for (int j=1;j<=jmax;j++){
            if (((FLAG[imax+1][j]&256)==256)&&(FLAG[imax][j]>511)){
                temperature[idx]=TEMP[imax][j];
                idx=idx+1;
            }
        }
        for (int i=1;i<=imax;i++){
            if (((FLAG[i][jmax+1]&256)==256)&&(FLAG[i][jmax]>511)){
                temperature[idx]=TEMP[i][jmax];
                idx=idx+1;
            }
        }
        for (int i=1;i<=imax;i++){
            if (((FLAG[i][0]&256)==256)&&(FLAG[i][1]>511)){
                temperature[idx]=TEMP[i][1];
                idx=idx+1;
            }
        }
        for (int j=1;j<=jmax;j++){
            for (int i=1;i<=imax;i++){
                if ((FLAG[i][j]&264)==264){//coupling + north fluid = 256 + 8
                    temperature[idx]=TEMP[i][j+1];
                    idx=idx+1;
                }
                else if ((FLAG[i][j]&260)==260){//coupling + north fluid = 256 + 4
                    temperature[idx]=TEMP[i][j-1];
                    idx=idx+1;
                }
            }
        }
        precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);
    }
    
void set_coupling_boundary(int imax, int jmax, double dx, double dy, double *heatflux, double **TEMP, int **FLAG)
    {
        int idx=0;
        for (int j=1;j<=jmax;j++){
            if (((FLAG[0][j]&256)==256)&&(FLAG[1][j]>511)){
                TEMP[0][j]=TEMP[1][j]+dx*heatflux[idx];
                idx=idx+1;
            }
        }
        for (int j=1;j<=jmax;j++){
            if (((FLAG[imax+1][j]&256)==256)&&(FLAG[imax][j]>511)){
                TEMP[imax+1][j]=TEMP[imax][j]+dx*heatflux[idx];
                idx=idx+1;
            }
        }
        for (int i=1;i<=imax;i++){
            if (((FLAG[i][jmax+1]&256)==256)&&(FLAG[i][jmax]>511)){
                TEMP[i][jmax+1]=TEMP[i][jmax]+dy*heatflux[idx];
                idx=idx+1;
            }
        }
        for (int i=1;i<=imax;i++){
            if (((FLAG[i][0]&256)==256)&&(FLAG[i][1]>511)){
                TEMP[i][0]=TEMP[i][1]+dy*heatflux[idx];
                idx=idx+1;
            }
        }
        for (int j=1;j<=jmax;j++){
            for (int i=1;i<=imax;i++){
                if ((FLAG[i][j]&264)==264){//coupling + north fluid = 256 + 8
                    TEMP[i][j]=TEMP[i][j+1]+dy*heatflux[idx];
                    idx=idx+1;
                }
                if ((FLAG[i][j]&260)==260){//coupling + south fluid = 256 + 4
                    TEMP[i][j]=TEMP[i][j-1]+dy*heatflux[idx];
                    idx=idx+1;
                }
            }
        }
    }

