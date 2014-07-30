 
       integer maxfrag, maxapf, maxsegs 
       parameter (maxfrag = 1000)  
       parameter (maxapf  = 100)  
       parameter (maxsegs = 1000)  
        
       logical frag_calc, fdrop_flag 
       integer nfrags, natoms_frag, watom_frag
       common /FRAGDEF/frag_calc, fdrop_flag, nfrags, 
     *                natoms_frag(maxfrag),
     *                watom_frag(maxfrag,maxapf) 
 
       integer fragAO, n_occ, n_aocc, n_bocc, n_virt, n_avirt, n_bvirt
       common /FORBS/fragAO(maxfrag), 
     *              n_occ(maxfrag), 
     *              n_aocc(maxfrag), 
     *              n_bocc(maxfrag),
     *              n_virt(maxfrag),
     *              n_avirt(maxfrag),
     *              n_bvirt(maxfrag)
 
       integer bocc_frag, eocc_frag, baocc_frag, eaocc_frag,
     *        bbocc_frag, ebocc_frag,
c 
     *        bvirt_frag, evirt_frag, bavirt_frag, eavirt_frag,
     *        bbvirt_frag, ebvirt_frag
 
       common /SEGRANGE/bocc_frag(maxfrag,maxapf), 
     *                 eocc_frag(maxfrag,maxapf),
     *                 baocc_frag(maxfrag,maxapf), 
     *                 eaocc_frag(maxfrag,maxapf),
     *                 bbocc_frag(maxfrag,maxapf), 
     *                 ebocc_frag(maxfrag,maxapf),
c 
     *                 bvirt_frag(maxfrag,maxapf), 
     *                 evirt_frag(maxfrag,maxapf),
     *                 bavirt_frag(maxfrag,maxapf), 
     *                 eavirt_frag(maxfrag,maxapf),
     *                 bbvirt_frag(maxfrag,maxapf), 
     *                 ebvirt_frag(maxfrag,maxapf)

       integer w_frag_AO, w_frag_occ, w_frag_aocc, w_frag_bocc, 
     *                   w_frag_virt, w_frag_avirt, w_frag_bvirt 

       common /SEG_FLIST/ w_frag_AO(maxsegs), 
     *                   w_frag_occ(maxsegs), 
     *                   w_frag_aocc(maxsegs), 
     *                   w_frag_bocc(maxsegs), 
     *                   w_frag_virt(maxsegs), 
     *                   w_frag_avirt(maxsegs), 
     *                   w_frag_bvirt(maxsegs) 
 
      integer wfrag 
      double precision afill 

      common /FRAG_FILL/wfrag, afill(maxapf,maxapf) 

