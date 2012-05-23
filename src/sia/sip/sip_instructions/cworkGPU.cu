#ifndef __CWORKGPU_CU__
#define __CWORKGPU_CU__
#include <cublas_v2.h>
#include "f77_name.h"
#include "f_types.h"
#include <stdio.h>

//**************************************************
// Configuration
//**************************************************
#define MAX_DIMS 20
#define SCRATCH_BUFFER_SIZE_MB 256
#define REORDER_TRANSFER_SIZE_MB 40
#define REORDER_STREAMS 32
#define REORDER_BLOCKS 48
#define REORDER_THREADS 512


//**************************************************
// Global Pointers
//**************************************************
double* scratch1;
double* scratch2;
double* scratch3;

__constant__ int ldDev[MAX_DIMS];
__constant__ int stepDev[MAX_DIMS];

cublasHandle_t cublasHandle;
cudaStream_t reorderStreams[REORDER_STREAMS];

//**************************************************
// Procedures
//**************************************************

#define imin(a, b) ((a < b)? a : b)

#ifdef __cplusplus
extern "C" {
#endif
void F77_NAME(cwork_gpu_alloc, CWORK_GPU_ALLOC)() {
  cudaMalloc(&scratch1, SCRATCH_BUFFER_SIZE_MB * 1024 * 1024);
  cudaMalloc(&scratch2, SCRATCH_BUFFER_SIZE_MB * 1024 * 1024);
  cudaMalloc(&scratch3, SCRATCH_BUFFER_SIZE_MB * 1024 * 1024);
  cublasCreate(&cublasHandle);

  for(int i = 0; i < REORDER_STREAMS; i++)
    cudaStreamCreate(&reorderStreams[i]);
}
}
#ifdef __cplusplus
extern "C" {
#endif
void F77_NAME(cwork_gpu_free, CWORK_GPU_FREE)() {
  cudaFree(scratch1);
  cudaFree(scratch2);
  cudaFree(scratch3);
  cublasDestroy(cublasHandle);

  for(int i = 0; i < REORDER_STREAMS; i++)
    cudaStreamDestroy(reorderStreams[i]);
}
}
__global__ void reorderScatter(double* newX, double* oldX, int ndims, int size, int offset = 0) {  
  int step = gridDim.x * blockDim.x;
  int oldIndex = blockIdx.x * blockDim.x + threadIdx.x + offset;
  int newIndex;
  int i;

  while(oldIndex < size) {
    newIndex = oldIndex / ldDev[0] * stepDev[0];
    for(i = 1; i < ndims; i++)
      newIndex += (oldIndex % ldDev[i-1]) / ldDev[i] * stepDev[i];    

    newX[newIndex] = oldX[oldIndex];
    oldIndex += step;
  }
}

__global__ void reorderGather(double* newX, double* oldX, int ndims, int size, int offset = 0) {
  int step = gridDim.x * blockDim.x;
  int newIndex = blockIdx.x * blockDim.x + threadIdx.x + offset;
  int oldIndex;
  int i;

  while(newIndex < size) {
    oldIndex = newIndex / ldDev[0] * stepDev[0];
    for(i = 1; i < ndims; i++)
      oldIndex += (newIndex % ldDev[i-1]) / ldDev[i] * stepDev[i];

    newX[newIndex] = oldX[oldIndex];
    newIndex += step;
  }
}

#ifdef __cplusplus
extern "C" {
#endif
void F77_NAME(cwork_gpu, CWORK_GPU) (f_double *y, f_int yOrder, f_int *yDims, f_int *yInds,
				     f_double *x1, f_int x1Order, f_int *x1Dims, f_int *x1Inds,
				     f_double *x2, f_int x2Order, f_int *x2Dims, f_int *x2Inds) {
  int cInds[MAX_DIMS];
  int xLds[MAX_DIMS];
  int xSteps[MAX_DIMS];
  int yLds[MAX_DIMS];
  int ySteps[MAX_DIMS];
  int step, lda, ldb, i, j;
  int x1Length, x2Length, yLength = 1;
  int numTransfers, transferLength;
  int cOrder = 0;

  printf("Y order: %d\nY dims: ", yOrder);
  for(i = 0; i < yOrder; i++)
    printf("%d ", yDims[i]);
  printf("\nY inds: ");
  for(i = 0; i < yOrder; i++)
    printf("%d", yInds[i]);

  printf("\nX1 order: %d\nX1 dims: ", x1Order);
  for(i = 0; i < x1Order; i++)
    printf("%d ", x1Dims[i]);
  printf("X1 inds: ");
  for(i = 0; i < x1Order; i++)
    printf("%d", x1Inds[i]);
  printf("\n");


  // Reverse y, x1, and x2 ordering to match row major
  for(i = 0; i < yOrder/2; i++) {
    j = yDims[i];
    yDims[i] = yDims[yOrder - 1 - i];
    yDims[yOrder - 1 - i] = j;
    j = yInds[i];
    yInds[i] = yInds[yOrder - 1 - i];
    yInds[yOrder - 1 - i] = j;
  }

  for(i = 0; i < x1Order/2; i++) {
    j = x1Dims[i];
    x1Dims[i] = x1Dims[x1Order - 1 - i];
    x1Dims[x1Order - 1 - i] = j;
    j = x1Inds[i];
    x1Inds[i] = x1Inds[x1Order - 1 - i];
    x1Inds[x1Order - 1 - i] = j;
  }


  for(i = 0; i < x2Order/2; i++) {
    j = x2Dims[i];
    x2Dims[i] = x2Dims[x2Order - 1 - i];
    x2Dims[x2Order - 1 - i] = j;
    j = x2Inds[i];
    x2Inds[i] = x2Inds[x2Order - 1 - i];
    x2Inds[x2Order - 1 - i] = j;
  }

  // determine which indices to contract
  for(i = 0; i < x1Order; i++)
    for(j = 0; j < x2Order; j++)
      if(x1Inds[i] == x2Inds[j])
	cInds[cOrder++] = x1Inds[i];

  // copy x1 into scratch3 and then reorder into scratch1
  {
    int cMask[MAX_DIMS] = {0};
    step = 1;
    for(i = x1Order - 1; i >= 0; i--) {
      xLds[i] = step;
      step *= x1Dims[i];
    }
    x1Length = step;
  
    for(i = 0; i < cOrder; i++) {
      for(j = 0; j < x1Order; j++) {
	if(cInds[i] == x1Inds[j]) {
	  step /= x1Dims[j];
	  xSteps[j] = step;
	  cMask[j] = 1;
	}
      }
    }
    lda = step;

    for(i = 0; i < x1Order; i++) {
      if(cMask[i])
	continue;

      step /= x1Dims[i];
      xSteps[i] = step;
      for(j = 0; j < yOrder; j++)
	if(yInds[j] == x1Inds[i])
	  ySteps[j] = step;
    }

    cudaMemcpyToSymbol(ldDev, xLds, sizeof(int) * x1Order);
    cudaMemcpyToSymbol(stepDev, xSteps, sizeof(int) * x1Order);

    numTransfers = (x1Length * sizeof(double) + REORDER_TRANSFER_SIZE_MB * 1024 * 1024 - 1) / (REORDER_TRANSFER_SIZE_MB * 1024 * 1024);
    transferLength = REORDER_TRANSFER_SIZE_MB * 1024 * 1024 / sizeof(double);

    for(i = 0; i < numTransfers; i++) {
      if(transferLength * i < x1Length) {
	cudaMemcpyAsync(scratch3 + transferLength * i, x1 + transferLength * i, 
			imin(transferLength, x1Length - transferLength * i) * sizeof(double), 
			cudaMemcpyHostToDevice, reorderStreams[i % REORDER_STREAMS]);
	reorderScatter<<<REORDER_BLOCKS, REORDER_THREADS, 0, reorderStreams[i % REORDER_STREAMS]>>>
	  (scratch1, scratch3, x1Order, imin(transferLength * (i + 1), x1Length), transferLength * i);
      }
    }

    cudaDeviceSynchronize();
  }

  // copy x2 into scratch3 and then reorder into scratch2
  {
    int cMask[MAX_DIMS] = {0};

    for(i = x2Order - 1; i >= 0; i--) {
      xLds[i] = step;
      step *= x2Dims[i];
    }
    x2Length = step;
     
    step = 1;
    for(i = cOrder - 1; i >= 0; i--) {
      for(j = 0; j < x2Order; j++) {
	if(x2Inds[j] == cInds[i]) {
	  xSteps[j] = step;
	  step *= x2Dims[j];
	  cMask[j] = 1;
	  break;
	}
      }
    }
    ldb = step;
  
    int ys = lda;
    for(i = x2Order - 1; i >= 0; i--) {
      if(cMask[i])
	continue;

      xSteps[i] = step;
      step *= x2Dims[i];
      for(j = 0; j < yOrder; j++)
	if(yInds[j] == x2Inds[i]) {
	  ySteps[j] = ys;
	  ys *= x2Dims[i];
	}
    }

    cudaMemcpyToSymbol(ldDev, xLds, sizeof(int) * x2Order);
    cudaMemcpyToSymbol(stepDev, xSteps, sizeof(int) * x2Order);

    numTransfers = (x2Length * sizeof(double) + REORDER_TRANSFER_SIZE_MB * 1024 * 1024 - 1) / (REORDER_TRANSFER_SIZE_MB * 1024 * 1024);

    for(i = 0; i < numTransfers; i++) {
      if(transferLength * i < x2Length) {
	cudaMemcpyAsync(scratch3 + transferLength * i, x2 + transferLength * i, 
			imin(transferLength, x2Length - transferLength * i) * sizeof(double), 
			cudaMemcpyHostToDevice, reorderStreams[i % REORDER_STREAMS]);
	reorderScatter<<<REORDER_BLOCKS, REORDER_THREADS, 0, reorderStreams[i % REORDER_STREAMS]>>>
	  (scratch2, scratch3, x2Order, imin(transferLength * (i + 1), x2Length), transferLength * i);
      }
    }

    cudaDeviceSynchronize();
  }

  // dGemm scratch1 and scratch2 into scratch 3
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, lda, x2Length / ldb, ldb, 
		&alpha, scratch1, lda, scratch2, ldb, &beta, scratch3, lda);
  cudaDeviceSynchronize();

  // reorder y from scratch3 to scratch1 and copy back from GPU
  for(i = yOrder - 1; i >= 0; i--) {
    yLds[i] = yLength;
    yLength *= yDims[i];
  }

  cudaMemcpyToSymbol(ldDev, yLds, sizeof(int) * x2Order);
  cudaMemcpyToSymbol(stepDev, ySteps, sizeof(int) * x2Order);

  numTransfers = (yLength * sizeof(double) + REORDER_TRANSFER_SIZE_MB * 1024 * 1024 - 1) / (REORDER_TRANSFER_SIZE_MB * 1024 * 1024);
  
  for(i = 0; i < numTransfers; i++) {
    if(transferLength * i < yLength) {
      reorderGather<<<REORDER_BLOCKS, REORDER_THREADS, 0, reorderStreams[i % REORDER_STREAMS]>>>
	(scratch1, scratch3, yOrder, imin(transferLength * (i + 1), yLength), transferLength * i);
      cudaMemcpyAsync(y + transferLength * i, scratch1 + transferLength * i, 
		      imin(transferLength, yLength - transferLength * i) * sizeof(double), 
		      cudaMemcpyDeviceToHost, reorderStreams[i % REORDER_STREAMS]);
    }
  }
  cudaDeviceSynchronize();
}
}
void cworkGPU(double* y, int yOrder, int* yDims, int* yInds,
	      double* x1, int x1Order, int* x1Dims, int* x1Inds,
	      double* x2, int x2Order, int* x2Dims, int* x2Inds) {
  F77_NAME(cwork_gpu, CWORK_GPU)(y, yOrder, yDims, yInds,
				 x1, x1Order, x1Dims, x1Inds,
				 x2, x2Order, x2Dims, x2Inds);
}

#endif // __CWORKGPU_CU__
