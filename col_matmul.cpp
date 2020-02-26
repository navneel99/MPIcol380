#include <iostream>
#include "mpi.h"
#include <time.h>
#include <chrono>
using namespace std;
using namespace chrono;
int IsEqual(float *A, float *B, int m, int n){
    int i,j;
    for (i = 0; i < m; i++ ){
       for (j =0; j < n; j++){
           if (A[i*n + j] != B[i*n +j]){
               return 0;
           }
       } 
    }
    return 1;
}


void Matrix_Multiply(float *A, float *B, float *C, int m, int n, int p){
	int i, j, k;
	for (i = 0; i < m; i++){
		for (j = 0; j < p; j++){
			C[i*p + j] = 0;
			for (k = 0; k < n; k++)
				C[i*p + j] += A[i*n + k] * B[k*p + j];
		}
	}
}	

void print_matrix(float *A, int n, int m){
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            cout<<A[i*m+j]<<" ";
        }
        cout<<"\n";
    }
}

int main(int argc, char *argv[])
{
    srand(time(0));
    int id;
    int p;
    int n=8000;
    int k=32;
    int ierr = MPI_Init(&argc, &argv);

    ierr = MPI_Comm_size(MPI_COMM_WORLD,&p);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&id);
    MPI_Status* status;
    double wtime;

    float *A, *B, *C, *C2,*Ans, *Part, *e;
    A = (float *)malloc(n*k*sizeof(float));
    B = (float *)malloc(n*k*sizeof(float));
    C = (float *)malloc(n*n*sizeof(float));
    C2 = (float *)malloc(n*n*sizeof(float));
    if (id == 0){ //PID 0 sets the matrices 
        int i,j;
        for (i=0; i<n;i++){
            for(j=0;j<k;j++){
                A[i*k+j] = ((float)rand()/RAND_MAX);
            }
        }
        for (i=0; i<k;i++){
            for(j=0;j<n;j++){
                B[i*n+j] = ((float)rand()/RAND_MAX);
            }
        }
    }
   
    Part = (float *)malloc((n/p) * k * sizeof(float));
    Ans = (float *)malloc((n/p) * n * sizeof(float));
    // // }
    if (id==0){
        wtime=MPI_Wtime();
    }
    MPI_Scatter(A,(n/p)*k,MPI_FLOAT,Part,(n/p)*k,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(B,n*k,MPI_FLOAT,0,MPI_COMM_WORLD);
    Matrix_Multiply(Part,B,Ans,(n/p),k,n);
    MPI_Gather(Ans,(n/p)*n,MPI_FLOAT,C,(n/p)*n,MPI_FLOAT,0,MPI_COMM_WORLD);
    
    if (id == 0){
        Matrix_Multiply(A+((n/p)*p*k),B,C+((n/p)*p*n),(n%p),k,n); //extra part
        // print_matrix(C,n,n);
    }

    if (id ==0){
        wtime=MPI_Wtime()-wtime;
        auto sertime = high_resolution_clock::now();
        Matrix_Multiply(A,B,C2,n,k,n);
        auto st = high_resolution_clock::now();
        auto du = duration_cast<milliseconds>(st-sertime);
        cout<<"Using COllective MPI: "<<wtime<<endl;
        cout<<"Using Serial Code: ";
        cout<<((float)du.count())/1000<<endl;
        cout <<"Sanity Check: "<< IsEqual(C,C2,n,n)<<endl;
        cout<<"-x-x-x-x-x-x-x-"<<endl;
        // Matrix_Multiply(A,B,C2,n,k,n);
    }
    // cout<< IsEqual()    
    ierr=MPI_Finalize();
    return 0;
}
