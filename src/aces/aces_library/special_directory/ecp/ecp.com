#ifndef _ECP_COM_
#define _ECP_COM_
C
C This file contain all the ECP variables that need to be known
C across multiple files.
C
      common /ECP_INT_VARS/Zlm(Lmnpwr), Lmnval(3,84),
     &                     Istart(0:Maxang),Iend(0:Maxang),
     &                     Ideg(0:Maxang),Lmf(Maxangpwr),
     &                     Lml(Maxangpwr),
     &                     Lmx(Lmnpwr),Lmy(Lmnpwr),Lmz(Lmnpwr),
     &                     Pi,Fpi,Sqpi2,Sqrt_Fpi,R_intcutoff

      common/ECP_POT_VARS/clp(Mxecpprim),zlp(Mxecpprim),
     &                    nlp(Mxecpprim),kfirst(Maxang,Max_centers),
     &                    klast(Maxang,Max_centers),
     &                    llmax(Max_centers)

      Common/ECP_INTGRD_VARS/Ideg_grd(0:Maxang),
     &                       Istart_grd(0:Maxang),Iend_grd(0:Maxang),
     &                       Lmnval_grd(7,Lmnmaxg)

      common /pseud / nelecp(Max_centers),ipseux(Max_centers),ipseud 

      common /nshel / expnt(ndi13),contr(ndi13,ndico),
     &                numcon(ndi13),katom(max_shells),
     &                ktype(max_shells),
     &                kprim(max_shells),kbfn(max_shells),
     &                kmini(max_shells),
     &                kmaxi(max_shells),nprims(max_shells),
     &                ndegen(max_shells),
     &                nshell,nbf

      Common /Qstore/Alpha,Beta,Xval
     
CSSS      Common /Ishuffle/Iorder(max_cbf),Inao(max_shells+1)

      Common /RadAng_sums/Rad_Sum(Maxang,Maxang), 
     &                    Ang_sum(Maxang,Maxang)
   
      Common /Fints/Fijk(0:4*Maxang,0:4*Maxang,0:4*Maxang)

      common /factorials/Fact(0:2*Maxang),Fac2(-1:4*Maxang),
     &                   Faco(0:2*Maxang),
     &                   Bcoefs(0:2*Maxang,0:2*Maxang),
     &                   Fprod(2*Maxang, 2*maxang)
  
#endif /* _ECP_COM_ */

