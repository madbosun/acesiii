      subroutine testnuc(NUCLEI,XX)

      implicit none 
      
      integer NUCLEI,i
      integer XX(1:NUCLEI)

           do i=1,NUCLEI

           write(6,*) XX(i), i , '<<----'
      
           enddo

        return
        end
