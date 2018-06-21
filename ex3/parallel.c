#include "parallel.h"
#include "helper.h"
#include <math.h>
#include "mpi.h"

void init_parallel (int iproc,
                    int jproc,
                    int imax,
                    int jmax,
                    int *myrank,
                    int *il,
                    int *ir,
                    int *jb,
                    int *jt,
                    int *rank_l,
                    int *rank_r,
                    int *rank_b,
                    int *rank_t,
                    int *omg_i,
                    int *omg_j,
                    int *num_proc)
{
    *num_proc=iproc*jproc;
    *omg_j=(*myrank/iproc)+1;
    *omg_i=(*myrank%iproc)+1;

    int quot_i=imax/iproc;
    int res_i=imax%iproc;
    int quot_j=jmax/jproc;
    int res_j=jmax%jproc;

    //we try as equi-distributed as possible
    //We must deal with the case imax is not divisable by iproc,or jmax is not divisable by jproc
    //for example, if imax=10, iproc=4, then we distribute like 3 3 2 2
    if (*omg_i<=res_i){
        *il=(*omg_i-1)*(quot_i+1)+1;
        *ir=(*omg_i)*(quot_i+1);
    }
    else{
        *il=(*omg_i-1)*quot_i+res_i+1;
        *ir=(*omg_i)*quot_i+res_i;
    }
    if (*omg_j<=res_j){
        *jb=(*omg_j-1)*(quot_j+1)+1;
        *jt=(*omg_j)*(quot_j+1);
    }
    else{
        *jb=(*omg_j-1)*quot_j+res_j+1;
        *jt=(*omg_j)*quot_j+res_j;
    }

    //neighbor ranks
    if (*omg_i==1){
        *rank_l=MPI_PROC_NULL;
    }
    else{
        *rank_l=*myrank-1;
    }

    if (*omg_i==iproc){
        *rank_r=MPI_PROC_NULL;
    }
    else{
        *rank_r=*myrank+1;
    }

    if (*omg_j==1){
        *rank_b=MPI_PROC_NULL;
    }
    else{
        *rank_b=*myrank-iproc;
    }

    if (*omg_j==jproc){
        *rank_t=MPI_PROC_NULL;
    }
    else{
        *rank_t=*myrank+iproc;
    }
}

void pressure_comm(double **P,
                    int il,
                    int ir,
                    int jb,
                    int jt,
                    int rank_l,
                    int rank_r,
                    int rank_b,
                    int rank_t,
                    double *bufSend,
                    double *bufRecv,
                    MPI_Status *status,
                    int chunk)
{
    //to left
    MPI_Barrier(MPI_COMM_WORLD);
    if (MPI_PROC_NULL!=rank_l){
        for (int j=jb;j<=jt;j++){
            bufSend[j-jb]=P[il][j];
        }
    }
    if (MPI_PROC_NULL==rank_r){
        MPI_Send(bufSend, chunk, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);
    }
    else if(MPI_PROC_NULL!=rank_l){
        MPI_Sendrecv(bufSend,chunk, MPI_DOUBLE, rank_l, 1, bufRecv, chunk, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL==rank_l){
        MPI_Recv(bufRecv, chunk, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL!=rank_r){
        for (int j=jb;j<=jt;j++){
            P[ir+1][j]=bufRecv[j-jb];
        }
    }

    //to right
    MPI_Barrier(MPI_COMM_WORLD);
    if (MPI_PROC_NULL!=rank_r){
        for (int j=jb;j<=jt;j++){
            bufSend[j-jb]=P[ir][j];
        }
    }
    if (MPI_PROC_NULL==rank_l){
        MPI_Send(bufSend, chunk, MPI_DOUBLE, rank_r, 2, MPI_COMM_WORLD);
    }
    else if(MPI_PROC_NULL!=rank_r){
        MPI_Sendrecv(bufSend,chunk, MPI_DOUBLE, rank_r, 2, bufRecv, chunk, MPI_DOUBLE, rank_l, 2, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL==rank_r){
        MPI_Recv(bufRecv, chunk, MPI_DOUBLE, rank_l, 2, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL!=rank_l){
        for (int j=jb;j<=jt;j++){
            P[il-1][j]=bufRecv[j-jb];
        }
    }

    //to top
    MPI_Barrier(MPI_COMM_WORLD);
    if (MPI_PROC_NULL!=rank_t){
        for (int i=il;i<=ir;i++){
            bufSend[i-il]=P[i][jt];
        }
    }
    if (MPI_PROC_NULL==rank_b){
        MPI_Send(bufSend, chunk, MPI_DOUBLE, rank_t, 3, MPI_COMM_WORLD);
    }
    else if(MPI_PROC_NULL!=rank_t){
        MPI_Sendrecv(bufSend,chunk, MPI_DOUBLE, rank_t, 3, bufRecv, chunk, MPI_DOUBLE, rank_b, 3, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL==rank_t){
        MPI_Recv(bufRecv, chunk, MPI_DOUBLE, rank_b, 3, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL!=rank_b){
        for (int i=il;i<=ir;i++){
            P[i][jb-1]=bufRecv[i-il];
        }
    }

    //to bottom
    MPI_Barrier(MPI_COMM_WORLD);
    if (MPI_PROC_NULL!=rank_b){
        for (int i=il;i<=ir;i++){
            bufSend[i-il]=P[i][jb];
        }
    }
    if (MPI_PROC_NULL==rank_t){
        MPI_Send(bufSend, chunk, MPI_DOUBLE, rank_b, 4, MPI_COMM_WORLD);
    }
    else if(MPI_PROC_NULL!=rank_b){
        MPI_Sendrecv(bufSend,chunk, MPI_DOUBLE, rank_b, 4, bufRecv, chunk, MPI_DOUBLE, rank_t, 4, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL==rank_b){
        MPI_Recv(bufRecv, chunk, MPI_DOUBLE, rank_t, 4, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL!=rank_t){
        for (int i=il;i<=ir;i++){
            P[i][jt+1]=bufRecv[i-il];
        }
    }
    //sync
    MPI_Barrier(MPI_COMM_WORLD);
}

void uv_comm (double**U,
                double**V,
                int il,
                int ir,
                int jb,
                int jt,
                int rank_l,
                int rank_r,
                int rank_b,
                int rank_t,
                double *bufSend,
                double *bufRecv,
                MPI_Status *status,
                int chunk)
{
    //to left
    MPI_Barrier(MPI_COMM_WORLD);
    if (MPI_PROC_NULL!=rank_l){
        for (int j=jb;j<=jt;j++){
            bufSend[2*(j-jb)]=V[il][j-1];
            bufSend[2*(j-jb)+1]=U[il][j];
        }
        bufSend[2*(jt-jb+1)]=V[il][jt];
    }
    if (MPI_PROC_NULL==rank_r){
        MPI_Send(bufSend, chunk, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);
    }
    else if(MPI_PROC_NULL!=rank_l){
        MPI_Sendrecv(bufSend,chunk, MPI_DOUBLE, rank_l, 1, bufRecv, chunk, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL==rank_l){
        MPI_Recv(bufRecv, chunk, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL!=rank_r){
        for (int j=jb;j<=jt;j++){
            V[ir+1][j-1]=bufRecv[2*(j-jb)];
            U[ir+1][j]=bufRecv[2*(j-jb)+1];
        }
        V[ir+1][jt]=bufRecv[2*(jt-jb+1)];
    }

    //to right
    MPI_Barrier(MPI_COMM_WORLD);
    if (MPI_PROC_NULL!=rank_r){
        for (int j=jb;j<=jt;j++){
            bufSend[2*(j-jb)]=V[ir][j-1];
            bufSend[2*(j-jb)+1]=U[ir-1][j];
        }
        bufSend[2*(jt-jb+1)]=V[ir][jt];
    }
    if (MPI_PROC_NULL==rank_l){
        MPI_Send(bufSend, chunk, MPI_DOUBLE, rank_r, 2, MPI_COMM_WORLD);
    }
    else if(MPI_PROC_NULL!=rank_r){
        MPI_Sendrecv(bufSend,chunk, MPI_DOUBLE, rank_r, 2, bufRecv, chunk, MPI_DOUBLE, rank_l, 2, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL==rank_r){
        MPI_Recv(bufRecv, chunk, MPI_DOUBLE, rank_l, 2, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL!=rank_l){
        for (int j=jb;j<=jt;j++){
            V[il-1][j-1]=bufRecv[2*(j-jb)];
            U[il-2][j]=bufRecv[2*(j-jb)+1];
        }
        V[il-1][jt]=bufRecv[2*(jt-jb+1)];
    }

    //to top
    MPI_Barrier(MPI_COMM_WORLD);
    if (MPI_PROC_NULL!=rank_t){
        for (int i=il;i<=ir;i++){
            bufSend[2*(i-il)]=U[i-1][jt];
            bufSend[2*(i-il)+1]=V[i][jt-1];
        }
        bufSend[2*(ir-il+1)]=U[ir][jt];
    }
    if (MPI_PROC_NULL==rank_b){
        MPI_Send(bufSend, chunk, MPI_DOUBLE, rank_t, 3, MPI_COMM_WORLD);
    }
    else if(MPI_PROC_NULL!=rank_t){
        MPI_Sendrecv(bufSend,chunk, MPI_DOUBLE, rank_t, 3, bufRecv, chunk, MPI_DOUBLE, rank_b, 3, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL==rank_t){
        MPI_Recv(bufRecv, chunk, MPI_DOUBLE, rank_b, 3, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL!=rank_b){
        for (int i=il;i<=ir;i++){
            U[i-1][jb-1]=bufRecv[2*(i-il)];
            V[i][jb-2]=bufRecv[2*(i-il)+1];
        }
        U[ir][jb-1]=bufRecv[2*(ir-il+1)];
    }

    //to bottom
    MPI_Barrier(MPI_COMM_WORLD);
    if (MPI_PROC_NULL!=rank_b){
        for (int i=il;i<=ir;i++){
            bufSend[2*(i-il)]=U[i-1][jb];
            bufSend[2*(i-il)+1]=V[i][jb];
        }
        bufSend[2*(ir-il+1)]=U[ir][jb];
    }
    if (MPI_PROC_NULL==rank_t){
        MPI_Send(bufSend, chunk, MPI_DOUBLE, rank_b, 4, MPI_COMM_WORLD);
    }
    else if(MPI_PROC_NULL!=rank_b){
        MPI_Sendrecv(bufSend,chunk, MPI_DOUBLE, rank_b, 4, bufRecv, chunk, MPI_DOUBLE, rank_t, 4, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL==rank_b){
        MPI_Recv(bufRecv, chunk, MPI_DOUBLE, rank_t, 4, MPI_COMM_WORLD, status);
    }
    if (MPI_PROC_NULL!=rank_t){
        for (int i=il;i<=ir;i++){
            U[i-1][jt+1]=bufRecv[2*(i-il)];
            V[i][jt+1]=bufRecv[2*(i-il)+1];
        }
        U[ir][jt+1]=bufRecv[2*(ir-il+1)];
    }
    //sync
    MPI_Barrier(MPI_COMM_WORLD);
}

void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}




