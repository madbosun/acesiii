//TODO unit tests

#include "mem_alloc.h"

#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <limits.h>
#include "f77_name.h"
#include "f_types.h"

#define BYTES_PER_MB 1048576
#define ALIGN 8

char * base_ptr = NULL;      //start of allocated memory
char * nxt_ptr = NULL;       //address of next free memory
long long top_int = 0LL;     //highest address (as long long int) of allocated memory + 1
long long total_bytes = 0LL; //total number of bytes allocated


/** 
  *  initializes dynamic memory system by obtaining requested memory via malloc
  *
  *  @param[in] megabytes amount of memory to be allocate and managed
  *  @param[in] sheap_flag  not currently implemented
  *  @param[out] ierr 0 if successful, -1 if not successful
  *
  *  If successful, this routine makes the memory management system active.
  *  If this subroutine is called again while the memory system is active, it will simply return.
*/
void F77_NAME(mem_alloc_init,MEM_ALLOC_INIT)(f_int *megabytes, f_int *sheap_flag, f_int *ierr){
#ifdef ALTIX
	assert (false); //altix memory management not implemented
#endif
//         printf("mem_alloc_init called  with %d megabytes\n", *megabytes);
	 assert( sizeof(long long) == sizeof(char *));  //verify assumption about size of pointers
	 *ierr = 0;
	 if (total_bytes == 0){  //only do this once unless free has been called to reset.
		total_bytes = (long long)(*megabytes) * BYTES_PER_MB;
		// total_bytes = (*megabytes) * BYTES_PER_MB;
		base_ptr = (char *)malloc(total_bytes);
		if (base_ptr == 0)
		{
			printf("malloc failed: %d megabytes",*megabytes);
			total_bytes = 0;
			*ierr = -1;
			return;
		}
		nxt_ptr = base_ptr;
		top_int = (long long)base_ptr + total_bytes;
		assert (top_int > 0); //check that using signed type (long long) for addresses OK
		return;
	 }
	return;
}

/** 
 * Allocates nwords elements of size element_size from managed array.  It is assumed that this subroutine
 * will be invoked from Fortran, so the parameters are described from the point of view of the Fortran code.
 * After return, the allocated memory is in base(ixmem):base(ixmem+nwords-1)
 *
 *  @param[in] base  the first element of a 1-element array declared and allocated in fortran
 *  @param[in] element_size the size of the elements in the array
 *  @param[out] ixmem integer*8 index into array indicating start of allocated memory.
 *  @param[in] heap  not currently implemented
 *  @param[out] ierr  value is 0 if successful
 */
void F77_NAME(mem_alloc,MEM_ALLOC)
      (char *base, f_int *nwords, f_int *element_size, long long *ixmem, f_int *heap, f_int *ierr){
	assert (total_bytes != 0);
	*ierr = 0;
    f_int esize = *element_size;
	long long nbytes = (long long)(*nwords) * esize;
	long long fiddle = ALIGN-1;
	long long nxt_int = (long long)nxt_ptr;
//        printf("mem_alloc called with nwords %d and esize %d, previous nxt_int is %lld.\n",*nwords,esize,nxt_int);
	long long alloc_int = (nxt_int + fiddle) & ~(fiddle); //align on ALIGN byte boundary
	nxt_int = alloc_int + nbytes;
//        printf("mem_alloc calculated nxt_int is %lld.\n",nxt_int);
	if (nxt_int > top_int){
		printf("Error base_ptr=%lld, nxt_int= %lld, top_int = %lld",base_ptr, nxt_int,top_int);
		*ierr = -1;
		return;
	}
	nxt_ptr = (char *)nxt_int;
	long long base_int = (long long)base;
 //       printf("mem_alloc done with nwords %d.base_ptr is %lld and nxt_ptr is %lld\n",*nwords, base_ptr, nxt_ptr);
 //       printf("size of f_int is %d\n.",sizeof(f_int));
 //      printf ("result of adding 1 to esize is %d\n", 1+esize);fflush(stdout); 
	*ixmem = (long long)((alloc_int - base_int)/esize + 1);
}


/** allocates nbytes bytes from the managed array.  It is assumed that this subroutine
 * will be invoked from C and has an interface similar to malloc
 *
 *  @param nbytes[n] number of bytes
 *  @return  address of allocated memory, or 0 if not successful
 */
char * mem_alloc_c(size_t nbytes){
	assert (total_bytes != 0);
	long long fiddle = ALIGN-1;
	long long nxt_int = (long long)nxt_ptr;
	long long alloc_int = (nxt_int + fiddle) & ~(fiddle); //align on ALIGN byte boundary
	nxt_int = alloc_int + nbytes;
        printf("mem_alloc_c called with %d bytes\n",nbytes);
	if (nxt_int >= top_int || nxt_int <= 0 ){  //check that enough memory available
		printf("nxt_int= %lld, top_int = %lld\n",nxt_int,top_int);
		alloc_int = 0;
		return NULL;
	}
	nxt_ptr = (char *)nxt_int;
	return (char *)alloc_int;
}


/** 
 *  frees all memory above addr so that it can be reallocated.
 *
 *  @param addr[in] starting address of memory to be free. Must be in allocated memory
 *  @param ierr[out] return code
 *
 *  Note:  mem_alloc_free will fail if called with invalid input.  There is no scenario where there is a normal
 *  return and ierr != 0.  We may want to change this interface.
*/
void F77_NAME(mem_alloc_free,MEM_ALLOC_FREE)(char *addr, f_int *ierr){
	*ierr = 0;
	assert(base_ptr <= addr && (long long) addr <= top_int); //given address must be in allocated memory
	nxt_ptr = addr;
}

/**  
*  mem_alloc_free_all returns all managed memory to OS and makes memory management inactive.
*  After calling this subroutine, the mem_alloc_init may be called to reinitialize the memory
*  management system.
*/
void F77_NAME(mem_alloc_free_all,MEM_ALLOC_FREE_ALL)(){
	if (base_ptr != 0) free(base_ptr);
	base_ptr = NULL;
	nxt_ptr = NULL;
	total_bytes = 0;
}


/**  
 *  Frees all memory allocated by mem_alloc so that it can be reallocated.  The
 *  memory system remains active.  Freed memory is not reinitialized to 0.
*/
void F77_NAME(mem_alloc_reset,MEM_ALLOC_RESET)(){
	nxt_ptr = base_ptr;
}


/** 
 *  returns the amount of memory allocated
 *
 *     *param nbused[out] amount of memory already allocated by mem_alloc
 */
void F77_NAME(mem_alloc_query, MEM_ALLOC_QUERY)(long long *nbused){
	*nbused = (long long)nxt_ptr - (long long)base_ptr;
}


/** 
 *  returns the number of objects of a particular size that can be allocated.
 *
 *     *param size[in] size of the object 
 *     *param number[out] number of objects that can be allocated
 */
void F77_NAME(mem_query_free, MEM_QUERY_FREE)(f_int* size,f_int *number){
        *number = (f_int) (total_bytes)/(*size);
}


/** 
 * returns the base address of the managed memory
 */
long long F77_NAME(get_mem_base_addr,GET_MEM_BASE_ADDR)(){
	return (long long)base_ptr;
}

/** returns address of base(ix) where base is an array with 'size'-byte elements
 *
 *   @param[in] base first element of array
 *   @param[in] ix   index of element whose address is to be returned
 *   @param[in] size size of elements of array
 *
 *   @returns address of base(ix)
 */
long long F77_NAME(c_loc64, C_LOC64)(char *base, long long *ix, f_int *size)
{
   long long addr;
   addr = (long long) base + (*ix-1)*(*size);
   return addr;
}

long long F77_NAME(get_max_heap_usage,GET_MAX_HEAP_USAGE) () 
{
   return total_bytes;

}

/*for debugging*/
void print_mem_values(){
	printf("base_ptr = %p\n",  (void *)base_ptr);
	printf("nxt_ptr =  %p\n",  (void *) nxt_ptr);
	printf("top_int = %lld\n", top_int);
	printf("total_bytes = %lld\n", total_bytes);
}
