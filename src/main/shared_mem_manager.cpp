#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <vector>
#include <string>
#include <fstream>
#include<mpi.h>
#include "f77_name.h"
#include "f_types.h"
#include <sys/mman.h>

//Global Variables
const char* const HOSTS_FILENAME = "hosts.out";
int globalrank;
int localrank;
MPI_Comm localcomm;
MPI_Comm leadercomm;
//shared memory ID   
int shm_id;  
//shared memory key
int shm_key;
void *buffer;
int buffersize;
struct shmid_ds shmid_struct;
using namespace std;

extern "C" {

/**  Function to split the processes so that all processes
 *   on a single host form a group
 *
 *     *param name[in] hostname of the machine
 *     *param name_len[in] length of hostname
 */
void split_processes( char* name, int name_len)
{
  long myhash;
  MPI_Comm_rank( MPI_COMM_WORLD, &globalrank );
  locale loc;
  const collate<char>& coll = use_facet<collate<char> >(loc);
  myhash = coll.hash(name,name+name_len);
  myhash = myhash % 1000000000;
  int color= abs((int)myhash);
  //This call creates local MPI communicator for every machine
  MPI_Comm_split( MPI_COMM_WORLD, color, 0, &localcomm );
  MPI_Comm_rank(localcomm, &localrank);

  //All leaders will establish their own leadercomm
  if(localrank==0) {
	MPI_Comm_split( MPI_COMM_WORLD, 0, 0, &leadercomm );
  }
  else {
        //Others still need to call Mpi_split, because MPI_Split
        //expects that all PEs in the MPI WORLD call it
	MPI_Comm_split( MPI_COMM_WORLD, 1, 0, &leadercomm );
  }
}


/** Function to save all the hostnames along with the shared memory identifiers
 *     *param hostname[in]
 *     *param hostname_len[in]
 *     *param shm_key[in]
*/
void save_hostname_and_shm_key(char* hostname, int hostname_len,int shm_key)
{
  int leader_count;
  //get number of local leaders 
  MPI_Comm_size(leadercomm, &leader_count);
  //All PEs in leadercomm will send their hostnames to PE 0
  if(globalrank != 0) {
	MPI_Send(hostname, hostname_len, MPI_CHAR, 0, 0, leadercomm);
  } else {
	char pe_hostname[MPI_MAX_PROCESSOR_NAME];
	int pe_hostname_len;
        MPI_Status status[leader_count];
	ofstream hosts_file;
	hosts_file.open (HOSTS_FILENAME,ios_base::out);
	hosts_file <<shm_key<<endl;
        //PE 0 will put its own hostname
	hosts_file.write(hostname, hostname_len); hosts_file << endl;
	//PE 0 will receive hostnames from all the other leaders
	for (int i=1;i<leader_count;i++)
	{
		MPI_Recv( &pe_hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, i, 0,leadercomm, &status[i]);
		MPI_Get_count( &status[i],  MPI_CHAR, &pe_hostname_len);
		hosts_file.write(pe_hostname, pe_hostname_len); hosts_file << endl;
	}
	hosts_file.close();
  }
}

/** This function allocates shared memory between nodes on the same machine
 *
 *     *param buf[out] The start address of the buffer in FORTRAN
 *     *param size[in] The size of memory in bytes, to be allocated
 *     *param offset[out] Offset that will be used by FORTRAN code when accessing this memory 
 */
void F77_NAME(malloc_shared_mem, MALLOC_SHARED_mem) (void *buf,int *size,long long *offset)
{

  long long off;
  char name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  //Get processors name (hostname)
  MPI_Get_processor_name(name, &name_len);
  //Split the processes based on hostname
  MPI_Comm_rank(MPI_COMM_WORLD, &globalrank );
  split_processes(name,name_len);

  if(localrank==0)
      printf("PE with global rank %d is local leader for %s\n",globalrank,name);

  //PE 0 will randonly choose a shm_key
  if (globalrank == 0){
	srand (time(NULL));
	shm_key = rand()%100000;
	 printf("SHM key is %d\n",shm_key);
  }
  //The shm_key will be broadcasted to all the PEs
  MPI_Bcast (&shm_key, 1,MPI_INT, 0 , MPI_COMM_WORLD );

  //All local leaders will call this
  //function to save the hostnames and the shm_key to a file
  if(localrank==0)
  save_hostname_and_shm_key(name,name_len,shm_key);
  MPI_Barrier(MPI_COMM_WORLD);
  sleep(30);
  //Increase the size by 8 bytes
  *size +=8;
  buffersize = *size;
  //Get the shared segment
  shm_id = shmget(shm_key, *size, 0666 | IPC_CREAT);
  if (shm_id < 0)
      printf("shmget error!\n");

  // Get the virtual address
  buffer= shmat(shm_id, NULL, 0);

  char *ptr = (char *) buffer;
  off = (long long) ( ptr - (char *) buf);
  off= off/4;
  off= off +1;
  *offset = off;

//  printf("Mem size %d, ptr allocated at %ld. buf points at %ld and offset is %ld\n",*size,ptr, buf,*offset);
//  printf("Max addr shud be %x\n", ptr + *size );
  return;
}

/** Function to broadcast the shared memory contents to all nodes
*
*/
void F77_NAME(bcast_shared_mem, BCAST_SHARED_MEM) ()
{
   if(localrank==0){
       MPI_Bcast (buffer,buffersize/sizeof(int),MPI_INT, 0 , leadercomm );
   }
}

/** Function to free the allocated shared memory
*
*/ 
void F77_NAME(free_shared_mem, FREE_SHARED_MEM) ()
{
  int rc=0;
  printf("free_shared_mem got called.\n");
  
  if(localrank==0)
  {
      //Make sure that noone is attached to memory
      while(rc!=-1) 
         rc=shmdt(buffer);

      //Delete the shared segment
      rc = shmctl(shm_id, IPC_RMID, &shmid_struct);
      if (rc==-1)
          printf("shmctl failed!\n");

  }
  else {
	//Detach from memory
	rc=shmdt(buffer);
  }

}


void F77_NAME(protect_shared_mem, PROTECT_SHARED_MEM) ()
{
printf("mprotech output: %d\n",mprotect(buffer,buffersize, PROT_READ));

}

void F77_NAME(detect_mem, DETECT_MEM) (int *a)
{

printf("Got call with address %ld  and sees %d\n",a,*a);

}
};
