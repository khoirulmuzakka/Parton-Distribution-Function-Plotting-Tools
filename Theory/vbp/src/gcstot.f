      FUNCTION GCSTOT(TA, TB, XA, XB, TAU)
C                                                   -=-=- gcstot

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
cjfo E866 xf mods
      common/E866xf/ixfx
      if(ixfx.eq.1)then
         gcstot=(tau**2+(tau-ta*tb)**2)/2./ta**3/tb**2/(tb+xa)
      else
         GCSTOT = XB * ( TAU + TA * TB ) * ( TAU**2 + 
     >             (TAU - TA * TB)**2 )
     >             / TA**3 / TB**2 / ( TA * XB + TB * XA ) / ( TB + XB )
      endif

      RETURN
C                                            ****************************
      ENTRY HCSTOT (TA, TB, XA, XB, TAU)
      if(ixfx.eq.1)then
         hcstot=(ta*(tb+xa)*(tb-xb)+2.*tau*(ta+tb))/2./ta**2/tb**2
     >         /(ta+tb)**2
      else
         HCSTOT =  TAU * (TAU + TA * TB )
     >               * ( TA * TB**2 * XA + TAU * 
     >               (TA * XB + 2.* TB * XA) )
     >               / ( TA * TB )**2 / (TA * XB + TB * XA )**3
      endif

      RETURN
C                        ****************************
      END
C                                                          =-=-= VbpAux
