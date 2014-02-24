      integer iter, nroots, chunk_type,RL,aux_vec_num,
     & targetroot
      double precision DavRvec(800,800),DavLvec(800,800),
     & DavRoots(800),DavRvec_initguess(800,800),
     & DavLvec_initguess(800,800),DavOverRmax(20),
     & DavOverLmax(20)
      common /EOMDAVINFO/ iter,nroots,DavRvec,DavLvec,
     & DavRoots,chunk_type,RL,DavRvec_initguess,
     & DavLvec_initguess,aux_vec_num,DavOverRmax,
     & DavOverLmax,targetroot

