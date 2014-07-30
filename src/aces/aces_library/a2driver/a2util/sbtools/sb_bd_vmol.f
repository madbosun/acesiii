
      blockdata sb_bd_vmol
      implicit none


c ilnbuf : the length of a buffer to read in chunks of integral files
c (There is no practical reason why this would be anything but 600.)

      integer           ilnbuf
      common /vmol_com/ ilnbuf
      save   /vmol_com/






      data ilnbuf /600/
      end

