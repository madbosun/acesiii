
c This routine defines the default values for the common block in
c include/ric_heap.com.

      blockdata bd_ric_heap
      implicit none
#define NO_EXTERNAL
#include "ric_heap.com"
      data dRICHeap /0.d0/
      data z_RICHeap, z_DerBMat, z_BMat, z_GMat, z_BTGInv /5*1/
      end

