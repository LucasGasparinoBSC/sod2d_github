module mod_geom

      contains

         !!subroutine char_length_qua()

         !!   implicit none

         !!   e12 = (/ abs(elcod(1,2)-elcod(1,1)), abs(elcod(2,2)-elcod(2,1)) /)
         !!   e23 = (/ abs(elcod(1,3)-elcod(1,2)), abs(elcod(2,3)-elcod(2,2)) /)
         !!   e34 = (/ abs(elcod(1,4)-elcod(1,3)), abs(elcod(2,4)-elcod(2,3)) /)
         !!   e41 = (/ abs(elcod(1,1)-elcod(1,4)), abs(elcod(2,1)-elcod(2,4)) /)

         !!   l12 = sqrt(dot_product(e12,e12))
         !!   l23 = sqrt(dot_product(e23,e23))
         !!   l34 = sqrt(dot_product(e34,e34))
         !!   l41 = sqrt(dot_product(e41,e41))

         !!   h = minval((/l12,l23,l34,l41/))

         !!end subroutine char_length_qua

end module mod_geom
