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
#define SCRATCH_BUFFER_SIZE_MB 100
#define REORDER_BLOCKS 48
#define REORDER_THREADS 512

void _alloc();
void _free();
void _cwork(f_double*  y, f_int ny, f_int*  yDims, f_int*  yInds,
	    f_double* x1, f_int n1, f_int* x1Dims, f_int* x1Inds,
	    f_double* x2, f_int n2, f_int* x2Dims, f_int* x2Inds); 
//**************************************************
// Global Pointers
//**************************************************
f_double* scratch1;
f_double* scratch2;
f_double* scratch3;

__constant__ f_int dimsDev[MAX_DIMS];
__constant__ f_int stepsDev[MAX_DIMS];

cublasHandle_t cublasHandle;

extern "C" {
  void F77_NAME(cwork_gpu_alloc, CWORK_GPU_ALLOC)() {
    _alloc();
  }
}

extern "C" {
  void F77_NAME(cwork_gpu_free, CWORK_GPU_FREE)() {
    _free();
  }
}

extern "C" {
  void F77_NAME(cwork_gpu, CWORK_GPU)
       (f_double*  y, f_int* ny, f_int*  nya, f_int* nyb, f_int*  yInds,
	f_double* x1, f_int* n1, f_int* nx1a, f_int* nx1b, f_int* x1Inds,
	f_double* x2, f_int* n2, f_int* nx2a, f_int* nx2b, f_int *x2Inds) {
    f_int yDims[MAX_DIMS], x1Dims[MAX_DIMS], x2Dims[MAX_DIMS];
    int i;
    for(i = 0; i < MAX_DIMS; i++) {
      yDims[i] = nyb[i] - nya[i] + 1;
      x1Dims[i] = nx1b[i] - nx1a[i] + 1;
      x2Dims[i] = nx2b[i] - nx2a[i] + 1;
    }
    _cwork( y, *ny,  yDims,  yInds,
	   x1, *n1, x1Dims, x1Inds,
	   x2, *n2, x2Dims, x2Inds);
  }
}

void _alloc() {
  int dev;
  if(cudaGetDevice(&dev))
    printf("GPU ERROR: cworkGPU: alloc: getDevice\n");
  printf("Allocating GPU memory on device %d\n", dev);
  
  if(cudaMalloc(&scratch1, SCRATCH_BUFFER_SIZE_MB * 1024 * 1024) ||
     cudaMalloc(&scratch2, SCRATCH_BUFFER_SIZE_MB * 1024 * 1024) ||
     cudaMalloc(&scratch3, SCRATCH_BUFFER_SIZE_MB * 1024 * 1024))
    printf("GPU ERROR: cworkGPU: alloc: cudaMalloc\n");
  
  cublasCreate(&cublasHandle);
}

void _free() {
  cudaFree(scratch1);
  cudaFree(scratch2);
  cudaFree(scratch3);
  cublasDestroy(cublasHandle);
}

__global__ void reorderScatter(double* newX, double* oldX, int ndims, int size) {  
  int blockstep = gridDim.x * blockDim.x;
  int oldIndex = blockIdx.x * blockDim.x + threadIdx.x;
  int newIndex;
  int t;
  int i;

  while(oldIndex < size) {
    t = oldIndex;

    newIndex = t % dimsDev[0] * stepsDev[0];
    t /= dimsDev[0];

    for(i = 1; i < ndims; i++) {
      newIndex += t % dimsDev[i] * stepsDev[i];
      t /= dimsDev[i];
    }

    newX[newIndex] = oldX[oldIndex];
    oldIndex += blockstep;
  }
}

__global__ void reorderGather(double* newX, double* oldX, int ndims, int size) {
  int blockstep = gridDim.x * blockDim.x;
  int newIndex = blockIdx.x * blockDim.x + threadIdx.x;
  int oldIndex;
  int t;
  int i;

  while(newIndex < size) {
    t = newIndex;
    oldIndex = t % dimsDev[0] * stepsDev[0];
    t /= dimsDev[0];

    for(i = 1; i < ndims; i++) {
      oldIndex += t % dimsDev[i] * stepsDev[i];
      t /= dimsDev[i];
    }

    newX[newIndex] = oldX[oldIndex];
    newIndex += blockstep;
  }
}

void _cwork(f_double*  y, f_int ny, f_int*  yDims, f_int*  yInds,
	    f_double* x1, f_int n1, f_int* x1Dims, f_int* x1Inds,
	    f_double* x2, f_int n2, f_int* x2Dims, f_int* x2Inds) {
  int steps[MAX_DIMS];
  int yIndsP[MAX_DIMS], yDimsP[MAX_DIMS];
  int x1IndsP[MAX_DIMS], x1DimsP[MAX_DIMS];
  int x2IndsP[MAX_DIMS], x2DimsP[MAX_DIMS];
  int step;
  int lda, ldb;
  int i, j;
  int c, k;
  int size;
  int nc = (n1 + n2 - ny) / 2;
  bool isContractedIndex;

  // determine permutations of x1, x2, and y
  c = 0;
  k = 0;
  for(i = 0; i < n1; i++) {
    isContractedIndex = false;

    for(j = 0; j < n2; j++) {
      if(x1Inds[i] == x2Inds[j]) {
	isContractedIndex = true;
	x1IndsP[n1 - nc + c] = x1Inds[i];
	x1DimsP[n1 - nc + c] = x1Dims[i];
	x2IndsP[c] = x2Inds[j];
	x2DimsP[c] = x2Dims[j];
	c++;
	break;
      }
    }

    if(!isContractedIndex) {
      x1IndsP[k] = x1Inds[i];
      x1DimsP[k] = x1Dims[i];
      yIndsP[k] = x1Inds[i];
      yDimsP[k] = x1Dims[i];
      k++;
    }
  }

  c = 0;
  for(i = 0; i < n2; i++) {
    for(j = 0; j < ny; j++) {
      if(x2Inds[i] == yInds[j]) {
	x2IndsP[nc + c] = x2Inds[i];
	x2DimsP[nc + c] = x2Dims[i];
	yIndsP[k] = yInds[j];
	yDimsP[k] = yDims[j];
	k++;
	c++;
      }
    }
  }
  
  // copy x1 into scratch3 and then reorder into scratch1
  step = 1;
  for(i = 0; i < n1; i++) {
    for(j = 0; j < n1; j++)
      if(x1Inds[j] == x1IndsP[i]) {
	steps[j] = step;
	break;
      }
    step *= x1DimsP[i];
  }
  size = step;      

  cudaMemcpyToSymbol(dimsDev, x1Dims, sizeof(f_int) * n1);
  cudaMemcpyToSymbol(stepsDev, steps, sizeof(f_int) * n1);

  cudaMemcpy(scratch3, x1, size * sizeof(f_double), cudaMemcpyHostToDevice);
  reorderScatter<<<REORDER_BLOCKS, REORDER_THREADS>>>(scratch1, scratch3, n1, size);
  cudaDeviceSynchronize();

  // copy x2 into scratch3 and then reorder into scratch2
  step = 1;
  for(i = 0; i < n2; i++) {
    for(j = 0; j < n2; j++)
      if(x2Inds[j] == x2IndsP[i]) {
	steps[j] = step;
	break;
      }
    step *= x2DimsP[i];
  }
  size = step;      
  
  cudaMemcpyToSymbol(dimsDev, x2Dims, sizeof(f_int) * n2);
  cudaMemcpyToSymbol(stepsDev, steps, sizeof(f_int) * n2);
  
  cudaMemcpy(scratch3, x2, size * sizeof(f_double), cudaMemcpyHostToDevice);
  reorderScatter<<<REORDER_BLOCKS, REORDER_THREADS>>>(scratch2, scratch3, n2, size);
  cudaDeviceSynchronize();
  
  // dGemm scratch1 and scratch2 into scratch 3
  double alpha = 1.0;
  double beta = 0.0;

  lda = 1;
  for(i = 0; i < n1 - nc; i++)
    lda *= x1DimsP[i];

  ldb = 1;
  for(i = 0; i < nc; i++)
    ldb *= x2DimsP[i];

  cublasDgemm(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, lda, size / ldb, ldb, 
	      &alpha, scratch1, lda, scratch2, ldb, &beta, scratch3, lda);
  cudaDeviceSynchronize();

  // reorder y from scratch3 to scratch1 and copy back from GPU
  step = 1;
  for(i = 0; i < ny; i++) {
    for(j = 0; j < ny; j++)
      if(yInds[j] == yIndsP[i]) {
	steps[j] = step;
	break;
      }
    step *= yDimsP[i];
  }
  size = step;      

  cudaMemcpyToSymbol(dimsDev, yDims, sizeof(f_int) * ny);
  cudaMemcpyToSymbol(stepsDev, steps, sizeof(f_int) * ny);

  reorderGather<<<REORDER_BLOCKS, REORDER_THREADS>>>(scratch1, scratch3, ny, size);
  cudaMemcpy(y, scratch1, size * sizeof(f_double), cudaMemcpyDeviceToHost);
  
  cudaDeviceSynchronize();
}

#endif // __CWORKGPU_CU__
