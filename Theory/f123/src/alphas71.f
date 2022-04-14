C ================================================================
      function alQCD(Q,iset)
C========================================================================
C     *** THIS IS A DUMMY FUNCTION **********
C     This links into a dummy routine for the Cteq4 alphas values
C F. March 8th 2017. I modified that such that it calls a hook defined in the c++ code
C========================================================================
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

c      alQCD = ascteq4(Q,iset)
c we don't need any set in the nCTEQ++ code
       alQCD = AlQCDDISF123LFortranWrapperHook (Q)

      RETURN
      END


