/*All function names prefixed with bas_ or BAS_ to allow simultaneous execution with original */

#ifndef SIPSHARED__C_MEM_ALLOC__H
#define SIPSHARED__C_MEM_ALLOC__H


#include <stddef.h>
#include "f77_name.h"
#include "f_types.h"

/** initializes dynamic memory system by obtaining requested memory from (via malloc)
  *
  *  @param[in] megabytes amount of memory to be allocate and managed
  *  @param[in] sheap_flag  not currently implemented
  *  @param[out] ierr 0 if successful, -1 if not successful
  *
  *  If successful, this routine makes the memory management system active.  If this subroutine is called
  *  again while the memory system is active, it will simply return.
*/
void F77_NAME(mem_alloc_init,MEM_ALLOC_INIT)(f_int *megabytes, f_int *sheap_flag, f_int *ierr);
//void bas_mem_alloc_init(f_int *megabytes, f_int *sheap_flag, f_int *ierr);

/** allocates nwords elements of size element_size from managed array.  It is assumed that this subroutine
 * will be invoked from Fortran, so the parameters are described from the point of view of the Fortran code.
 *  After return, the allocated memory is in base(ixmem):base(ixmem+nwords-1)
 *
 *  @param[in] base  the first element of a 1-element array declared and allocated in fortran
 *  @param[in] element_size the size of the elements in the array
 *  @param[out] ixmem integer*8 index into array indicating start of allocated memory.
 *  @param[in] heap  not currently implemented
 *  @param[out] ierr  value is 0 if successful
 */
void F77_NAME(mem_alloc,MEM_ALLOC)
      (char *base, f_int *nwords, f_int *element_size, long long *ixmem, f_int *heap, f_int *ierr);


/** allocates nbytes bytes from the managed array.  It is assumed that this subroutine
 * will be invoked from C and has an interface similar to malloc
 *
 *  @param nbytes[n] number of bytes
 *  @return  address of allocated memory, or 0 if not successful
 */
char * mem_alloc_c(size_t nbytes);

/**  frees all memory above addr so that it can be reallocated.
 *
 *  @param addr[in] starting address of memory to be free. Must be in allocated memory
 *  @param ierr[out] return code
 *
 *  Note:  mem_alloc_free will fail if called with invalid input.  There is no scenario where there is a normal
 *  return and ierr != 0.  We may want to change this interface.
*/

void F77_NAME(mem_alloc_free,MEM_ALLOC_FREE)(char *addr, f_int *ierr);

/**  mem_alloc_free_all returns all managed memory to OS and makes memory management inactive.
*  After calling this subroutine, the mem_alloc_init may be called to reinitialize the memory
*  management system.
*/
void F77_NAME(mem_alloc_free_all,MEM_ALLOC_FREE_ALL)();

/**  Frees all memory allocated by mem_alloc so that it can be reallocated.  The
 *   memory system remains active.  Freed memory is not reinitialized to 0.
*/
void F77_NAME(mem_alloc_reset,MEM_ALLOC_RESET)();

/** returns the amount of memory allocated
 *
 *     *param nbused[out] amount of memory already allocated by mem_alloc
 */
void F77_NAME(mem_alloc_query, MEM_ALLOC_QUERY)(long long *nbused);



/** returns the base address of the managed memory
 */
long long F77_NAME(get_mem_base_addr,GET_MEM_BASE_ADDR)();


/** returns address of base(ix) where base is an array with 'size'-byte elements
 *
 *   @param[in] base first element of array
 *   @param[in] ix   index of element whose address is to be returned
 *   @param[in] size size of elements of array
 *
 *   @returns address of base(ix)
 */
long long F77_NAME(c_loc64, C_LOC64)(char *base, long long *ix, f_int *size);

/*for debugging*/
void print_mem_values();

#endif /*SIPSHARED__C_MEM_ALLOC__H*/
