      real*8 function getheff(array,index)
cjp because in calling routine array is declared as integer!
      real*8 array
      integer index
      dimension array(*)
      getheff = array(index)
      return
      end
