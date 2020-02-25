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
    int n=1000;
    int k=32;
    int ierr = MPI_Init(&argc, &argv);

    ierr = MPI_Comm_size(MPI_COMM_WORLD,&p);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&id);
    MPI_Status* status;
    double wtime;

    float *A, *B, *C, *C2,*Ans, *Part;
    A = (float *)malloc(n*k*sizeof(float));
    B = (float *)malloc(n*k*sizeof(float));
    C = (float *)malloc(n*n*sizeof(float));
    C2 = (float *)malloc(n*n*sizeof(float));
    if (id == p-1){
        Part = (float *)malloc( ((n/p) + (n%p)) * k * sizeof(float));
        Ans = (float *)malloc( ((n/p)+(n%p)) * n * sizeof(float));
    }else{
        Part = (float *)malloc((n/p) * k * sizeof(float));
        Ans = (float *)malloc((n/p) * n * sizeof(float));
    }

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
        wtime = MPI_Wtime();
        for(i=1;i<p;i++){
            if (i==p-1){
                MPI_Send(A +(p-1)*(n/p)*k,((n/p) + n%p)*k, MPI_FLOAT,i,1,MPI_COMM_WORLD);
            }else{
                MPI_Send(A + i*(n/p)*k,(n/p)*k,MPI_FLOAT,i,1,MPI_COMM_WORLD);
            }
            MPI_Send(B,n*k,MPI_FLOAT,i,1,MPI_COMM_WORLD);
        }
        Matrix_Multiply(A,B,C,(n/p),k,n);
        for(i=1;i<p;i++){
            if (i==p-1){
                MPI_Recv(C+((p-1)*(n/p)*n),((n/p)+(n%p))*n,MPI_FLOAT,i,MPI_ANY_TAG,MPI_COMM_WORLD,status);
            }else{
                MPI_Recv(C+(i*(n/p)*n),(n/p)*n,MPI_FLOAT,i,MPI_ANY_TAG,MPI_COMM_WORLD,status);
            }
        }
        wtime = MPI_Wtime() - wtime;
    }else{
        if (id == p-1){
            MPI_Recv(Part,((n/p)+(n%p))*k,MPI_FLOAT,0,MPI_ANY_TAG,MPI_COMM_WORLD,status);
        }else{
            MPI_Recv(Part,(n/p)*k,MPI_FLOAT,0,MPI_ANY_TAG,MPI_COMM_WORLD,status);
        }
        MPI_Recv(B,n*k,MPI_FLOAT,0,MPI_ANY_TAG,MPI_COMM_WORLD, status);
        if (id==p-1){
            Matrix_Multiply(Part,B,Ans,(n/p)+(n%p) ,k,n);
            MPI_Send(Ans,((n/p)+(n%p))*n,MPI_FLOAT,0,1,MPI_COMM_WORLD);
        }else{
            Matrix_Multiply(Part,B,Ans,(n/p),k,n);
            MPI_Send(Ans,(n/p)*n,MPI_FLOAT,0,1,MPI_COMM_WORLD);
        }
    }
    if (id == 0){
        auto sertime = high_resolution_clock::now();
        Matrix_Multiply(A,B,C2,n,k,n);
        auto st = high_resolution_clock::now();
        auto du = duration_cast<milliseconds>(st-sertime);
        cout<<"Using MPI: "<<wtime<<endl;
        cout<<"Using Serial Code: ";
        cout<<((float)du.count())/1000<<endl;
        cout << IsEqual(C,C2,n,n)<<endl;
    }
    ierr = MPI_Finalize();

    return 0;
}
