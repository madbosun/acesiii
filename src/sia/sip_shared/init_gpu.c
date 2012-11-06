#ifndef _GPU_INIT_C
#define _GPU_INIT_C

#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include "f_types.h"
#include "f77_name.h"

void F77_NAME(init_gpu, INIT_GPU) (f_int* devid) {
  char *filename;
  FILE* file;

  filename = getenv("PBS_GPUFILE");
  file = fopen(filename, "r");
  (*devid) = -1;

  if(file != NULL) {
    printf("gpufile opened\n");
    while(fscanf(file, "gpu%d", devid) == 0) {
      fscanf(file, "%*c");
    }
    printf("Found device id = %d\n", *devid);
    if(*devid >= 0);
//      cudaSetDevice(*devid);
  }
}

#endif //_GPU_INIT_C
