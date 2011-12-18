c--------------------------------------------------------------------------
c   This file has the parameters about the shared memory between PEs
c--------------------------------------------------------------------------
 
      integer shared_buf(1)
      integer shared_size
      integer*8 shared_mem_offset
      common/shared_mem_params/shared_mem_offset,
     * shared_size
      common/shared_mem_block/shared_buf


