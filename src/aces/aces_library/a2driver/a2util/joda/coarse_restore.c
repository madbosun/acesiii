
/*
 * This routine restores the precious files for geometry optimizations and
 * finite difference calculations from szSaveDir/(OLD|CURRENT).
 *
 * It is called EVERY TIME joda runs since the first joda does not know
 * whether the program is a restart and the files are simply missing.
 */

#include <stdio.h>	/* for *printf */
#include <unistd.h>	/* for access */
#include <stdlib.h>	/* for system */
#include <string.h>	/* for strcmp */
#include <strings.h>	/* for strncpy */

#include "f77_name.h"   /* for F77_NAME */
#include "f_types.h"    /* for f_int */

void
F77_NAME(coarse_restore,COARSE_RESTORE)
(const char * szSaveDirIn, f_int * iErr)
{

#define CMDLEN 256
#define BAKLEN 128
    char szSaveDir[80], szCmd[CMDLEN], szBak[BAKLEN];
    if (!strncpy(szSaveDir,szSaveDirIn,80))
    {
        printf("\nERROR: unable to copy directory string\n");
        *iErr=1;
        return;
    }

 /* TEST IF szSaveDir/OLD EXISTS */
    snprintf(&szBak[0],BAKLEN,"%s/OLD",szSaveDir);
    if (access(szBak,F_OK))
    {
     /* SWITCH TO szSaveDir/CURRENT */
        snprintf(&szBak[0],BAKLEN,"%s/CURRENT",szSaveDir);
        if (access(szBak,F_OK))
        {
         /* no directories in szSaveDir -> firstrun */
            *iErr=0;
            return;
        }
    }

 /* TEST IF szBak IS USEABLE */
    if (access(szBak,R_OK|X_OK))
    {
        printf("ERROR: %s does not have the proper permissions\n",szBak);
        *iErr=1;
        return;
    }

 /* RESTORE PRECIOUS FILES */
    {
#include "coarse_precious.h" /* for coarse_precious[] */
#define BUFLEN 80
#define FILLEN 80
        int i=0; char szBuf[BUFLEN], szFile[FILLEN];
        FILE *fpFILES = fopen("FILES","r"); /* check FILES for overrides */

        printf("Files restored from %s:",szBak);
        while (*coarse_precious[i])
        {
            strncpy(&szFile[0],coarse_precious[i],FILLEN);
            if (fpFILES)
            {
             /* scan FILES for a matching line */
                while (fgets(&szBuf[0],BUFLEN,fpFILES))
                {
                    char *cz = &szBuf[0]; int i = 0;
                    while ((*cz != ' ') && (*cz != '\n')) cz++; *cz++ = NULL;
                    if (*cz)
                    { char *c = cz; while (*c != '\n') c++; *c = NULL; }
                    if (!strncmp(&szBuf[0],&szFile[0],FILLEN))
                    { strncpy(&szFile[0],cz,FILLEN); i = 1; }
                    if (i) break;
                }
                rewind(fpFILES);
            }
            snprintf(&szCmd[0],CMDLEN,"%s/%s",szBak,coarse_precious[i]);
            if (access(&szCmd[0],F_OK))
                printf(" !%s",coarse_precious[i]);
            else
            {
                snprintf(&szCmd[0],CMDLEN,"/bin/cp -f %s/%s %s",
                         szBak,coarse_precious[i],szFile);
                if (system(&szCmd[0]))
                {
                    printf("\nERROR: unable to restore %s/%s to %s\n",
                           szBak,coarse_precious[i],szFile);
                    *iErr=1;
                    return;
                }
                printf(" %s",coarse_precious[i]);
            }
            i++;
        }
        printf("\n");
        if (fpFILES) fclose(fpFILES);
    }

    *iErr=0;
    return;
}

