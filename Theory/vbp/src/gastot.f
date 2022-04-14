      FUNCTION GASTOT(TA, TB, XA, XB, TAU)
C                                                   -=-=- gastot
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
cjfo  E866 xf mods
      common/E866xf/ixfx
      if(ixfx.eq.1)then
         gastot=(ta+tb)*(tau**2+(ta*tb)**2)
     >         /2./(ta*tb)**2/(ta+xb)/(tb+xa)
      else
         GASTOT = ( TAU + TA*TB ) * ( TAU**2 + (TA*TB)**2 )
     >          / (TA*TB)**2 / (TA +XA) / (TB + XB)
      endif

      RETURN
c                                            ****************************

      ENTRY HASTOT (TA, TB, XA, XB, TAU)
cjfo
      if(ixfx.eq.1)then
         hastot=-1./ta/tb/(ta+tb)
      else
         HASTOT = - 2. * TAU * ( TAU + TA * TB)
     >          / TA / TB / ( TA * XB + TB * XA )**2
      endif
      RETURN
C                        ****************************
      END
C                                                          =-=-= GluQrk
