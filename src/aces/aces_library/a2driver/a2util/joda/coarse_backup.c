
/*
 * This routine saves the precious files from geometry optimizations and
 * finite difference calculations to szSaveDir/CURRENT.
 */

#include <stdio.h>	/* for *printf */
#include <stdlib.h>	/* for system */
#include <unistd.h>	/* for access */
#include <sys/types.h>	/* for mkdir */
#include <sys/stat.h>	/* for mkdir */
#include <string.h>	/* for strcmp */
#include <strings.h>	/* for strncpy */

#include "f77_name.h"	/* for F77_NAME */
#include "f_types.h"	/* for f_int */

void
F77_NAME(coarse_backup,COARSE_BACKUP)
(const char * szSaveDirIn, f_int * iErr)
{

#define CMDLEN 256
#define CURLEN 128
#define OLDLEN 128
    char szSaveDir[80], szCmd[CMDLEN], szCurrent[CURLEN], szOld[OLDLEN];
    if (!strncpy(szSaveDir,szSaveDirIn,80))
    {
        printf("\nERROR: unable to copy directory string\n");
        *iErr=1;
        return;
    }
    snprintf(&szCurrent[0],CURLEN,"%s/CURRENT",szSaveDir);
    snprintf(&szOld[0],    OLDLEN,"%s/OLD",    szSaveDir);

 /* CONDITION szSaveDir */
    if (access(szSaveDir,F_OK))
    {
        printf("Attempting to create %s for coarse grain restart . . . ",
               szSaveDir);
        if (mkdir(&szSaveDir[0],S_IRWXU))
        {
            printf("\nERROR: unable to create backup directory\n");
            *iErr=1;
            return;
        }
        printf("done\n");
    }
    else /* szSaveDir exists */
    {
        if (access(szSaveDir,R_OK|W_OK|X_OK))
        {
            printf("ERROR: %s does not have the proper permissions\n",
                   szSaveDir);
            *iErr=1;
            return;
        }
    }

 /* REMOVE szSaveDir/OLD */
    if (!access(&szOld[0],F_OK))
    {
        printf("WARNING: %s already exists.\n"
               "         The backup directory may have been contaminated.\n"
               "         Attempting to remove directory . . . ",szOld);
        snprintf(&szCmd[0],CMDLEN,"/bin/rm -rf %s",szOld);
        if (system(&szCmd[0]))
        {
            printf("\nERROR: unable to remove %s\n",szOld);
            *iErr=1;
            return;
        }
        printf("done\n");
    }

 /* MOVE CURRENT TO OLD */
    if (!access(&szCurrent[0],F_OK))
    {
        printf("Renaming %s to %s . . . ",szCurrent,szOld);
        if (rename(&szCurrent[0],&szOld[0]))
        {
            printf("\nERROR: unable to rename %s\n",szCurrent);
            *iErr=1;
            return;
        }
        printf("done\n");
    }

 /* MAKE CURRENT AND SAVE PRECIOUS FILES */
    if (mkdir(&szCurrent[0],S_IRWXU))
    {
        printf("ERROR: unable to create %s\n",szCurrent);
        *iErr=1;
        return;
    }
    else
    {
#include "coarse_precious.h" /* for coarse_precious[] */
#define BUFLEN 80
#define FILLEN 80
        int i=0; char szBuf[BUFLEN], szFile[FILLEN];
        FILE *fpFILES = fopen("FILES","r"); /* check FILES for overrides */

	printf("Files copied to %s:",szCurrent);
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
            if (access(&szFile[0],F_OK))
                printf(" !%s",coarse_precious[i]);
            else
            {
                snprintf(&szCmd[0],CMDLEN,"/bin/cp -f %s %s/%s",
                         szFile,szCurrent,coarse_precious[i]);
                if (system(&szCmd[0]))
                {
                    printf("\nERROR: unable to copy %s to %s/%s\n",
                           szFile,szCurrent,coarse_precious[i]);
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

 /* REMOVE szSaveDir/OLD */
    if (!access(&szOld[0],F_OK))
    {
        printf("Removing old back up directory . . . ");
        snprintf(&szCmd[0],CMDLEN,"/bin/rm -rf %s",szOld);
        if (system(&szCmd[0]))
        {
            printf("\nERROR: unable to remove %s\n",szOld);
            *iErr=1;
            return;
        }
        printf("done\n");
    }

    *iErr=0;
    return;
}

