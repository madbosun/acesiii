#ifndef _RIC_HEAP_COM_
#define _RIC_HEAP_COM_
c ric_heap.com : begin

c This common block contains the heap address and array indices for
c processing RICs. New additions must be initialized to 1 in
c bd_ric_heap.F and set in init_ric_heap.F.

#ifndef NO_EXTERNAL
      external bd_ric_heap
#endif

      double precision dRICHeap(1)
      integer          z_RICHeap, z_DerBMat, z_BMat, z_GMat, z_BTGInv

      common /ric_heap_com/ dRICHeap,
     &                      z_RICHeap, z_DerBMat, z_BMat, z_GMat,
     &                      z_BTGInv
      save   /ric_heap_com/

c ric_heap.com : end
#endif /* _RIC_HEAP_COM_ */
