
#include <stdio.h>
#include <unistd.h>

void
#ifdef _UNICOS
     BWCHDIR
#else
#  ifdef C_SUFFIX
     bwchdir_
#  else
     bwchdir
#  endif /* C_SUFFIX */
#endif /* _UNICOS */
(long * iproc)
{
    /*first let's simply use path names from a given file*/
    FILE *f;
    char mypath[256];
    long i;

    if ((f=fopen("MPIDIRS","r"))==NULL)
    {
        fprintf(stderr,"Cannot open MPIDIRS, %d\n",*iproc);
        _exit(10);
    }
    for (i=*iproc; i>0; --i)
        fscanf(f,"%s",mypath);
    printf("Process %d WORKDIR %s\n",*iproc,mypath);
    i=chdir(mypath);
    if (i)
    {
        perror("Cannot chdir");
        _exit(10);
    }
    fclose(f);
}

