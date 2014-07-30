#ifndef _ECP_COM_
#define _ECP_COM_
C
C This file contain all the ECP variables that need to be known 
C across multiple files.
C
#include "mxatms.par"
#include "baslims.par"

      Dimension Bcoefs(0:Maxang,0:Maxang),Dfact(-1:Maxang)
   
      common /ECPVAR/ Bcoefs, Dfact
      
#endif /* _ECP_COM_ */
