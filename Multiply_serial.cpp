#include <iostream>
#include <mpi.h>
#include <time.h>

using namespace std;

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

int main(){
    
    srand(time(0));
    int r = 10;
    int c = 10;
    float *arr = (float *)malloc(r*c*sizeof(float));
    // for (int i=0;i<r;i++){
    //     arr[i] = (double *)malloc(c * sizeof(c * sizeof(double)));
    // }
    for (int i = 0; i <r;i++){
        for(int j=0;j<c;j++){
            arr[i*c+j] = rand();
        }
    }
    cout << IsEqual(arr,arr,r,c)<<endl;

    return 0;

}