
      FUNCTION myisnan (val)
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        myisnan = 0
        if(isnan(val).eqv..true.) myisnan =1
      END FUNCTION

      SUBROUTINE myisnan_d(val, vald, intout)
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
        vald = 0.0
        intout = 0
        if(isnan(val).eqv..true.) then
          intout= 1
        end if
      END SUBROUTINE
